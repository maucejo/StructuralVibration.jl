"""
    agwn(x, snr_dB, rst = true)

Adds a Gaussian White Noise (AWGN) to a signal `x` with a given SNR.

# Inputs
- `x`: signal - Matrix{ComplexF64}
- `snr_dB`: signal to noise ratio [dB] - Float64
- `rst`: reset the random number generator - Bool

# Output
- `y`: noisy signal - Matrix{ComplexF64}

# Example
```julia-repl
julia> y = agwn(x, 25.)
```
"""
function agwn(x, snr_dB, rst = true)
    # Reset the RNG if required
    if rst
        rng = MersenneTwister(1000)
    end

    N, L = size(x)                          # Dimensions des données
    SNR = 10^(snr_dB/10.)                   # SNR en échelle linéaire
    En = sum(abs2, x, dims = 2)/L           # Énergie du signal
    V = En/SNR                              # Variance du bruit

    σ = sqrt.(V)                           # Écart-type
    n = σ.*randn(rng, eltype(x), N, L)    # Bruit

    return x .+ n
end

"""
    estimated_SNR(x, var)

Estimates the SNR of a signal `x` with a given variance `var`.

# Inputs
- `x`: signal - Matrix{ComplexF64}
- `var`: variance - Float64

# Output
- `SNR`: signal to noise ratio [dB] - Float64

# Example
```julia-repl
julia> SNR = estimated_SNR(x, 1e-3)
```
"""
function estimated_SNR(x, var)
    L = size(x, 2)
    En = sum(abs2, x, dims = 2)/L

    SNR = En./var

    return 10log10.(SNR)
end

"""
    varest(x, method = :bayes)

Estimates the noise variance of a signal `x` using the method specified by `method`.

# Inputs
- `x`: signal - Matrix{ComplexF64}
- `method`: method to estimate the noise variance - Symbol
    - `:bayes`: Bayesian denoising
    - `:derrico`: method proposed by John D'Errico in [1]

# Output
- `noisevar`: noise variance - Vector{Float64}

# Example
```julia-repl
julia> noisevar = varest(x, :bayes)
```

# Reference
[1] John D'Errico (2023). Estimatenoise (https://www.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise), MATLAB Central File Exchange. Retrieved December 7, 2023.
"""
function varest(x, method = :bayes)
    if method == :bayes
        return noisevar1D(x)
    elseif method == :derrico
        return estimatenoise(x)
    else
        error("Method not implemented")
    end
end

"""
    noisevar1D(x)

Estimates the noise variance of a signal `x` using Bayesian denoising.

# Inputs
- `x`: signal - Matrix{ComplexF64}

# Output
- `noisevar`: noise variance - Vector{Float64}

# Example
```julia-repl
julia> noisevar = noisevar1D(x)
```
"""
function noisevar1D(x)
    # Initialisation
    ndim = ndims(x)

    if ndim == 1
        noisevar = noisevar1D_(x)
    elseif ndim == 2
        nx = size(x, 1)
        noisevar = Vector{Float64}(undef, nx)
        @inbounds @views for idx ∈ 1:nx
             noisevar[idx] = noisevar1D_(x[idx, :])
        end
    end

    return noisevar
end

"""
    noisevar1D_(x)

Estimates the noise variance of a signal `x` using Bayesian regularization.

Note: This function is not intended to be used directly

# Input
- `x`: signal - Vector{Float64}

# Output

- `noisevar`: noise variance - Float64
"""
function noisevar1D_(x)
    # Valeurs propres de la matrice de lissage d'ordre 1
    n = length(x)
    s = Vector{Float64}(undef, n)
    z = Vector{eltype(x)}(undef, n)

    @. s = 2(cos((0:n-1)π/n) - 1.)
    @. s[s == 0.] = 1e-8
    s² = s.^2

    # Calul de la DCT-2
    z .= dct(x)

    lb = -5.
    ub = 5.
    # u0 = 0.

    optimfunc = L -> func!(L, z, s²)
    res = optimize(optimfunc, lb, ub)
    # res = optimize(optimfunc, [0.], LBFGS(), autodiff = :forward)
    λ = 10^only(Optim.minimizer(res))

    fₖ = @. (1. + λ*s²)/s²
    gₖ = mean(@. abs2(z)/fₖ)

    # γₐ = 1e-3
    # γₛ = 1e-3
    # βₐ = 1e-3
    # βₛ = 1e-3
    # gₖ = (sum(@. abs2(z)/fₖ) + βₐ + βₛ/λ)/(n + γₐ - γₛ)

    return λ*gₖ
end

"""
    func!(L, z, s²)

Function to be optimized in `noisevar1D_`.

Note: This function is not intended to be used directly

# Inputs

- `L`: parameter to be optimized - Float64
- `z`: signal - Vector{Float64}
- `s²`: eigenvalues of the smoothing matrix - Vector{Float64}

# Output
- `f`: function to be optimized - Float64
"""
function func!(L, z, s²)
    n = length(z)
    fₖ = @. (1. + 10. ^L*s²)/s²
    gₖ = mean(@. abs2(z)/fₖ)

    return sum(log, fₖ) + (n - 2)*log(gₖ)

    # γₐ = 1e-3
    # γₛ = 1e-3
    # βₐ = 1e-3
    # βₛ = 1e-3
    # gₖ = (sum(@. abs2(z)/fₖ) + βₐ .+ βₛ ./(10. .^L))/(n + γₐ - γₛ)

    # return only(sum(log, fₖ) .+ (n + γₐ - γₛ - 2)*log.(gₖ) .+ (1. + γₛ)*log.(10. .^L))
end

"""
    estimatenoise(x)

Estimates the noise variance of a signal `x` using the method proposed by John D'Errico in [1].

# Input
- `x`: signal - Matrix{ComplexF64}

# Output
- `noisevar`: noise variance - Vector{Float64}

# Example
```julia-repl
julia> noisevar = estimatenoise(x)
```

# Reference
[1] John D'Errico (2023). Estimatenoise (https://www.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise), MATLAB Central File Exchange. Retrieved December 7, 2023.
"""
function estimatenoise(x)
    if !isreal(x)
        return estimatenoise_(real(x)) + estimatenoise_(imag(x))
    end

    return estimatenoise_(x)
end

"""
    estimatenoise_(x)

Estimates the noise variance of a signal `x` using the method proposed by John D'Errico in [1].

Note: This function is not intended to be used directly

# Input
- `x`: signal - Vector{Float64}

# Output
- `noisevar`: noise variance - Float64

# Example
```julia-repl
julia> noisevar = estimatenoise_(x)
```
"""
function estimatenoise_(x)
    nd, ns = size(x)
    # The idea here is to form a linear combination of successive elements
    # of the series. If the underlying form is locally nearly linear, then
    # a [1 -2 1] combination (for equally spaced data) will leave only
    # the noise remaining. Next, if we assume the measurement noise was
    # iid, N(0,σ²), then we can try to back out the noise variance.
    nfda = 6
    np = 14
    fda = Vector{Vector{Float64}}(undef, nfda)
    perc = Vector{Float64}(undef, np)
    # Normalization to unit norm
    fda[1] = [1., -1.]./√(2.)
    fda[2] = [1., -2., 1.]./√(6.)
    fda[3] = [1., -3., 3., -1.]./√(20.)
    fda[4] = [1., -4., 6., -4., 1.]./√(70.)
    fda[5] = [1., -5., 10., -10., 5., -1.]./√(252.)
    fda[6] = [1., -6., 15., -20., 15., -6., 1.]./√(924.)

    # Compute an interquantile range, like the distance between the 25 and 75 points. This trims off the trash at each end, potentially corrupted if there are discontinuities in the curve. It also deals simply with a non-zero mean in this data. Actually do this for several different interquantile ranges, then take a median.
    # NOTE: While I could have used other methods for the final variance estimation, this method was chosen to avoid outlier issues when the curve may have isolated discontinuities in function value or a derivative. The following points correspond to the central 90, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, and 20 percent ranges.
    perc .= [0.05; range(0.1, 0.4, step = 0.025)]
    z = @. √(2)*erfinv(1. - 2perc)

    noisevar = Vector{Float64}(undef, nd)
    σₑ = fill(NaN, nd, nfda)
    Q = Matrix{Float64}(undef, nd, np)
    @inbounds @views for i ∈ 1:nfda
        fdai = fda[i]
        posnd = (i+1):ns
        ntrim = ns - i
        noisedata = Matrix{Float64}(undef, nd, ntrim)
        p = Vector{Float64}(undef, ntrim)
        @inbounds for j ∈ 1:nd
            xj = x[j, :]
            noisedata[j, :] .= conv(xj, fdai)[posnd]
        end
        # noisedata = conv([1.], fdai, x)[:, posnd] # Less time efficient

        if ntrim ≥ 2
            # Sorting will provide the necessary percentiles after interpolation.
            sort!(noisedata, dims = 2,  alg = InsertionSort)
            p .= (0.5 .+ collect(Float64, 1:ntrim))./(ntrim + 0.5)

            @inbounds for k ∈ 1:nd
                itp = LinearInterpolation(noisedata[k, :], p)
                @. Q[k, :] = (itp(1 - perc) - itp(perc))/2z
            end

            # Trim off any nans first, since if the series was short enough, some of those percentiles were not present.
            notnan = findall(@. !isnan(Q[1, :]))

            # Our noise std estimate is given by the median of the interquantile range(s). This is an ad hoc, but hopefully effective, way of estimating the measurement noise present in the signal.
            @inbounds for j ∈ 1:nd
                σₑ[j, i] = median(Q[j, notnan])
            end
        end
    end

    # Drop those estimates which failed for lack of enough data
    notnan = findall(@. !isnan(@view σₑ[1, :]))

    # Use median of these estimates to get a noise estimate.
    @inbounds @views for j ∈ 1:nd
        noisevar[j] = median(σₑ[j, notnan])^2.
    end

    # Use an adhoc correction to remove the bias in the noise estimate. This correction was determined by examination of a large number of random samples.
    noisevar ./= (1. + 15(ns + 1.225)^(-1.245))

    return noisevar
end