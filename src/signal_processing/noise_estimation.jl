abstract type NoiseEstimation end
abstract type OptimFamily <: NoiseEstimation end

"""
    BayesianEst

Bayesian noise estimation

**Fields**
* `prior::Symbol`: Prior distribution of the noise parameters to estimate - Symbol
    * `:invgamma`: Inverse gamma distribution (default)
    * `:uniform`: Uniform distribution
"""
@show_data struct BayesianEst <: OptimFamily
    prior :: Symbol

    BayesianEst(prior = :invgamma) = new(prior)
end

"""
    GCVEst

Generalized Cross-Validation (GCV) noise estimation

This method has been proposed by Garcia in [1]

**Fields**
* `nothing`

**Reference**

[1] Garcia, D. (2010). Robust smoothing of gridded data in one and higher dimensions with missing values. Computational Statistics and Data Analysis, 54(5), 1167-1178
"""
struct GCVEst <: OptimFamily end

"""
    LCurveEst

L-curve noise estimation

This method is based on the method proposed by Hansen in [1]

**Fields**
* `nothing`

**Reference**

[1] Hansen, P. C. (1999). The L-curve and its use in the numerical treatment of inverse problems. Computational Inverse Problems in Electrocardiology, 119-142
"""
struct LCurveEst <: OptimFamily end

"""
    DerricoEst

D'Errico noise estimation

This method has been proposed by John D'Errico in [1]

**Fields**
* `nothing`

**Reference**

[1] John D'Errico (2023). Estimatenoise, MATLAB Central File Exchange. Retrieved December 7, 2023. (https://www.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise)
"""
struct DerricoEst <: NoiseEstimation end

"""
    varest(x, method::NoiseEstimation; batch_size = 0, summary = mean)
    varest(x; batch_size = 0, summary = mean)

Estimates the noise variance of a signal `x` using a given method

**Inputs**
* `x`: Signal
* `method`: Noise estimation method
    * `BayesianEst`: Bayesian noise estimation (To be implemented)
    * `GCVEst`: Generalized Cross-Validation (GCV) noise estimation
    * `LCurveEst`: L-curve noise estimation
    * `DerricoEst`: D'Errico noise estimation (default)
* `batch_size::Int`: Batch size for batch processing (default = 0)
* `summary`: Summary function for batch processing (default = mean)

**Output**
* `noisevar`: Noise variance - Vector{Real}


**Note**
* varest(x; ...) = varest(x, DerricoEst(); ...)
"""
function varest(x, method::NoiseEstimation; batch_size::Int = 0, summary = mean)
    varestfun = let
        if method isa OptimFamily
            x -> varest_optim(x, method)
        elseif method isa DerricoEst
            varest_derrico
        end
    end

    if batch_size == 0
        return varestfun(x)
    else
        # Check type of x
        x isa Vector ? x = reshape(x, 1, :) : nothing

        nr, nc = size(x)

        # Create batches
        batches_candidate = Vector{typeof(x)}[x[:, i:min(i + batch_size - 1, end)] for i in 1:batch_size:nc]

        # Valid batches
        len_b = length.(batches_candidate)
        valid_batches = findall(len_b .≥ batch_size)

        batches = batches_candidate[valid_batches]
        nb = length(batches)

        noisevar = similar(real(x), nr, nb)
        for (b, batch) in enumerate(batches)
            noisevar[:, b] .= varestfun(batch)
        end

        return vec(summary(noisevar, dims = 2))
    end
end

# Default method
varest(x; batch_size::Int = 0, summary = mean) = varest(x, DerricoEst(), batch_size = batch_size, summary = summary)

"""
    varest_bayesian(x, method::OptimFamily)

Estimates the noise variance of a signal `x` using Bayesian denoising.

**Inputs**
* `x`: Signal
* `method`: Noise estimation method
    * `BayesianEst`: Bayesian noise estimation
    * `GCVEst`: Generalized Cross-Validation (GCV) noise estimation
    * `LCurveEst`: L-curve noise estimation

**Output**
* `noisevar`: Noise variance - Vector{Float64}
"""
function varest_optim(x, method::OptimFamily)
    varestfun = let
        if method isa BayesianEst
            x -> varest_bayesian(x, method.prior)
        elseif method isa GCVEst
            varest_gcv
        elseif method isa LCurveEst
            varest_lcurve
        end
    end

    ndim = ndims(x)
    if ndim == 1
        noisevar = [varestfun(x)]
    elseif ndim == 2
        nx = size(x, 1)
        noisevar = similar(real(x), nx)
        @views for (idx, xi) in enumerate(eachrow(x))
             noisevar[idx] = varestfun(xi)
        end
    end

    return noisevar
end

"""
    varest_bayesian(x, prior)

Estimates the noise variance of a signal `x` using Bayesian regularization.

Note: This function is not intended to be used directly

# Input
* `x`: Signal - Vector{Float64}
* `prior::Symbol`: Prior distribution of the noise parameters to estimate - Symbol
    * `:invgamma`: Inverse gamma distribution (default)
    * `:uniform`: Uniform distribution

**Output**
* `noisevar`: Noise variance
"""
function varest_bayesian(x, prior = :invgamma)
    # Eigenvalues of the second difference matrix
    n = length(x)
    s = similar(real(x))
    z = similar(x)

    @. s = 2(1. - cos((0:n-1)π/n))
    @. s[s == 0.] = 1e-8
    s² = s.^2

    # Calulation of DCT-2
    z .= dct(x)

    lb = -5.
    ub = 5.

    optimfunc = Λ -> bayesfun!(Λ, z, s², prior)
    res = optimize(optimfunc, lb, ub)
    λ = 10^only(Optim.minimizer(res))

    fₖ = @. (1. + λ*s²)/s²
    if prior == :uniform
        vₖ = mean(@. abs2(z)/fₖ)
    elseif prior == :invgamma
        # Define α depending on the type of x
        eltype(x) == Complex{Float64} ? α = 1. : α = 0.5

        γₐ = 1.
        γₛ = 1.
        βₐ = 1e-10
        βₛ = 1e-10
        vₖ = (α*sum(@. abs2(z)/fₖ) + βₐ + βₛ/λ)/(2. + γₛ + γₐ +  α*n)
    end

    return λ*vₖ
end

"""
    bayesfun!(Λ, z, s²)

Function to be optimized in `varest_bayesian`

Note: This function is not intended to be used directly

**Inputs**
* `Λ`: Parameter to be optimized
* `z`: Signal
* `s²`: Eigenvalues of the smoothing matrix

**Output**
* `J`: Functional to be minimized
"""
function bayesfun!(Λ, z, s², prior)
    n = length(z)

    fₖ = @. (1. + 10. ^Λ*s²)/s²

    # Define α depending on the type of x
    eltype(z) == Complex{Float64} ? α = 1. : α = 0.5

    if prior == :uniform
        vₖ = mean(@. abs2(z)/fₖ)

        J = α*sum(log, fₖ) + (α*n - 2)*log(vₖ)
    elseif prior == :invgamma
        γₐ = 1.
        γₛ = 1.
        βₐ = 1e-10
        βₛ = 1e-10
        M = 2. + γₛ + γₐ +  α*n
        vₖ = (α*sum(@. abs2(z)/fₖ) + βₐ .+ βₛ ./(10. .^Λ))/M

        J = α*sum(log, fₖ) .+ (M - 2)*log.(vₖ) .+ (1. + γₛ)*log.(10. .^Λ)
    end

    return only(J)
end

"""
    varest_gcv(x)

Estimates the noise variance of a signal `x` using the Generalized Cross-Validation (GCV) method as proposed by Garcia in [1]

# Input
* `x`: Signal

**Output**
* `noisevar`: Noise variance

**Reference**

[1] Garcia, D. (2010). Robust smoothing of gridded data in one and higher dimensions with missing values. Computational Statistics and Data Analysis, 54(5), 1167-1178.
"""
function varest_gcv(x)
    n = length(x)
    s = similar(real(x))
    z = similar(x)

    @. s = 2(1. - cos((0:n-1)π/n))
    @. s[s == 0.] = 1e-8
    s² = s.^2

    z .= dct(x)

    hmin = 1e-6
    hmax = 0.99
    lb = (((1 + sqrt(1+8*hmax.^(2)))/4hmax^2)^2 - 1)/16
    ub = (((1 + sqrt(1+8*hmin.^(2)))/4hmin^2)^2 - 1)/16

    optimfunc = Λ -> gcvfun!(Λ, z, s²)
    res = optimize(optimfunc, log10(lb), log10(ub))
    λ = 10^only(Optim.minimizer(res))

    fₖ = @. 1/(1 + λ*s²) - 1.

    return mean(abs2, z.*fₖ)
end

"""
    gcvfun!(Λ, z, s²)

Function to be optimized in `varest_gcv`

Note: This function is not intended to be used directly

**Inputs**
* `Λ`: Parameter to be optimized
* `z`: Signal
* `s²`: Eigenvalues of the smoothing matrix

**Output**
* `J`: Functional to be minimized
"""
function gcvfun!(Λ, z, s²)
    fₖ = @. 1/(1 + 10^Λ*s²) - 1.
    noisevar = mean(abs2, z.*fₖ)
    J = noisevar/mean(fₖ)^2

    return only(J)
end

"""
    varest_lcurve(x)

Estimates the noise variance of a signal `x` using the L-curve method as proposed by Hansen in [1]

# Input
* `x`: Signal

**Output**
* `noisevar`: Noise variance

**Reference**

[1] Hansen, P. C. (1999). The L-curve and its use in the numerical treatment of inverse problems. Computational Inverse Problems in Electrocardiology, 119-142
"""
function varest_lcurve(x)
    n = length(x)
    s = similar(real(x))
    z = similar(x)

    @. s = 2(cos((0:n-1)π/n) - 1.)
    @. s[s == 0.] = 1e-8
    s² = s.^2

    z .= dct(x)

    npoints = 500
    lb = -22.
    ub = 22.
    Λ = 10 .^LinRange(lb, ub, npoints)

    ρ = similar(s, npoints)
    η = similar(s, npoints)

    for (i, λ) in enumerate(Λ)
        fₖ = @. 1/(1 + λ*s²)
        ρ[i] = mean(abs2, @. z*(fₖ - 1.))
        η[i] = mean(abs2, @. z*s*fₖ)
    end

    k = curvature(log.(ρ), log.(η), Λ, method = :linear)
    notnan = findall(@. !isnan(k))
    pos_max = argmax(k[notnan])

    return ρ[notnan][pos_max]
end

"""
    varest_derrico(x)

Estimates the noise variance of a signal `x` using the method proposed by John D'Errico in [1].

# Input
* `x`: Signal

**Output**
* `noisevar`: Noise variance

**Reference**

[1] John D'Errico (2023). Estimatenoise (https://www.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise), MATLAB Central File Exchange. Retrieved December 7, 2023.
"""
function varest_derrico(x)

    noisevar = similar(real(x), size(x, 1))
    if !isreal(x)
        noisevar .= varest_derrico_real(real(x)) .+ varest_derrico_real(imag(x))
    else
        noisevar .= varest_derrico_real(x)
    end

    return noisevar
end

"""
    varest_derrico_real(x)

Estimates the noise variance of a real signal `x` using the method proposed by John D'Errico in [1].

# Input
* `x`: Real signal

**Output**
* `noisevar`: Noise variance

**Reference**

[1] John D'Errico (2023). Estimatenoise (https://www.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise), MATLAB Central File Exchange. Retrieved December 7, 2023.
"""
function varest_derrico_real(x)
    ndim = ndims(x)
    ndim == 1 ? x = reshape(x, 1, :) : nothing

    nd, ns = size(x)

    # NOTE: The comments are those of John D'Errico in the original code.
    # The idea here is to form a linear combination of successive elements
    # of the series. If the underlying form is locally nearly linear, then
    # a [1 -2 1] combination (for equally spaced data) will leave only
    # the noise remaining. Next, if we assume the measurement noise was
    # iid, N(0,σ²), then we can try to back out the noise variance.
    nfda = 6
    np = 14
    fda = Vector{eltype(x)}[similar(x, i+1) for i in 1:nfda]
    # Normalization to unit norm
    fda[1] .= [1., -1.]./√2.
    fda[2] .= [1., -2., 1.]./√6.
    fda[3] .= [1., -3., 3., -1.]./√20.
    fda[4] .= [1., -4., 6., -4., 1.]./√70.
    fda[5] .= [1., -5., 10., -10., 5., -1.]./√252.
    fda[6] .= [1., -6., 15., -20., 15., -6., 1.]./√924.

    # Compute an interquantile range, like the distance between the 25 and 75 points. This trims off the trash at each end, potentially corrupted if there are discontinuities in the curve. It also deals simply with a non-zero mean in this data. Actually do this for several different interquantile ranges, then take a median.
    # NOTE: While I could have used other methods for the final variance estimation, this method was chosen to avoid outlier issues when the curve may have isolated discontinuities in function value or a derivative. The following points correspond to the central 90, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, and 20 percent ranges.
    perc = [0.05; range(0.1, 0.4, step = 0.025)]
    z = @. √2*erfinv(1. - 2perc)

    noisevar = similar(x, nd)
    σe = fill(NaN, nd, nfda)
    Q = similar(x, nd, np)
    @views for (i, fdai) in enumerate(fda)
        posnd = (i+1):ns
        ntrim = ns - i
        noisedata = similar(x, nd, ntrim)
        p = similar(x, ntrim)
        for (j, xj) in enumerate(eachrow(x))
            noisedata[j, :] .= conv(xj, fdai)[posnd]
        end

        if ntrim ≥ 2
            # Sorting will provide the necessary percentiles after interpolation.
            sort!(noisedata, dims = 2,  alg = InsertionSort)
            p .= (0.5 .+ collect(Float64, 1:ntrim))./(ntrim + 0.5)

            for (k, nk) in enumerate(eachrow(noisedata))
                itp = LinearInterpolation(nk, p, extrapolation = ExtrapolationType.Extension)
                @. Q[k, :] = (itp(1 - perc) - itp(perc))/2z
            end

            # Trim off any nans first, since if the series was short enough, some of those percentiles were not present.
            notnan = findall(@. !isnan(Q[1, :]))

            # Our noise std estimate is given by the median of the interquantile range(s). This is an ad hoc, but hopefully effective, way of estimating the measurement noise present in the signal.
            σe[:, i] .= median(Q[:, notnan], dims = 2)
        end
    end

    # Drop those estimates which failed for lack of enough data
    notnan = findall(@. !isnan(@view σe[1, :]))

    # Use median of these estimates to get a noise estimate.
    noisevar .= median(σe[:, notnan], dims = 2).^2.

    # Use an adhoc correction to remove the bias in the noise estimate. This correction was determined by examination of a large number of random samples.
    noisevar ./= (1. + 15(ns + 1.225)^(-1.245))

    # return ndim == 1 ? only(noisevar) : noisevar

    return noisevar
end

"""
    estimated_SNR(x, var)

Estimates the SNR of a signal `x` with a given variance `var`.

**Inputs**
* `x`: Signal
* `var`: Variance

**Output**
* `SNR`: signal to noise ratio [dB] - Float64
"""
function estimated_SNR(x, var)
    En = vec(mean(abs2, x, dims = ndims(x)))

    SNR = En./var

    return 10log10.(SNR)
end
