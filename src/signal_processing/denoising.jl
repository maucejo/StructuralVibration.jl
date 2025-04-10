abstract type DenoisingMethod end
abstract type RegDenoising <: DenoisingMethod end

"""
    BayesDenoising(prior = :invgamma)

Regularization-based denoising using Bayesian inference for estimating the regularization parameter

**Fields**
* `prior::Symbol`: Prior distribution over the hyperparameters used for computing the regularization parameter
    * `:invgamma`: Inverse Gamma distribution
    * `:uniform`: Uniform distribution
"""
@show_data struct BayesDenoising <: RegDenoising
    prior::Symbol

    BayesDenoising(prior = :invgamma) = new(prior)
end

"""
    GCVDenoising()

Regularization-based denoising method using the GCV method for estimating the regularization parameter

**Fields**
* `nothing`

**Reference**

[1] Garcia, D. (2010). Robust smoothing of gridded data in one and higher dimensions with missing values. Computational Statistics and Data Analysis, 54(5), 1167-1178
"""
struct GCVDenoising <: RegDenoising end

"""
    LCurveDenoising()

Regularization-based denoising method using the L-curve method for estimating the regularization parameter

**Fields**
* `nothing`

**Reference**

[1] Hansen, P. C. (1999). The L-curve and its use in the numerical treatment of inverse problems. Computational Inverse Problems in Electrocardiology, 119-142
"""
struct LCurveDenoising <: RegDenoising end

"""
    KalmanDenoising(; rts = false)

Kalman filter denoising method

**Fields**
* `rts`: Flag to enable the Rauch-Tung-Striebel smoother
"""
@show_data struct KalmanDenoising <: DenoisingMethod
    rts::Bool

    KalmanDenoising(; rts = false) = new(rts)
end

"""
    denoising(y::AbstractArray, alg)

Denoises a signal `y`

**Inputs**
* `y`: Noisy signal
* `alg`: Denoising method
    * `BayesDenoising`: Bayesian Regularization denoising method (to be implemented)
    * `GCVDenoising`: GCV denoising method
    * `LCurveDenoising`: L-curve denoising method
    * `KalmanDenoising`: Kalman filter denoising method

**Output**
* `x`: Denoised signal
"""
function denoising(y::T, alg::RegDenoising) where {T <: AbstractArray}
    reg_denoise = let
        if alg isa BayesDenoising
            (x, s², U) -> bayes_denoise(x, s², U, alg.prior)
        elseif alg isa GCVDenoising
            gcv_denoise
        elseif alg isa LCurveDenoising
            lcurve_denoise
        end
    end


    ndim = ndims(y)
    n = size(y, ndim)

    # Initialization
    s = similar(real(y), n)
    z = similar(y)
    U = similar(y, n, n)

    # Eigenvalues of the second difference matrix
    @. s = 2(1. - cos((0:n-1)π/n))
    @. s[s == 0.] = 1e-8
    s² = s.^2

    # Corresponding eigenvectors
    scale = similar(s)
    scale[1] = sqrt(1/n)
    scale[2:end] .= sqrt(2/n)
    U = [cos((k - 1)*(2j - 1)π/2n) for j = 1:n, k = 1:n].*scale'

    # Calculation of DCT-2
    z .= dct(y, ndim)

    x = similar(y)
    if ndim == 1
        x .= reg_denoise(z, s², U)
    else
        for (idz, zi) in enumerate(eachrow(z))
            x[idz, :] .= reg_denoise(zi, s², U)
        end
    end

    return x
end

function bayes_denoise(z, s², U, prior)
    lb = -5.
    ub = 5.

    optimfunc = Λ -> bayesfun!(Λ, z, s², prior)
    res = optimize(optimfunc, lb, ub)
    λ = 10^only(Optim.minimizer(res))

    return U*(z./(1 .+ λ*s²))
end

function gcv_denoise(z, s², U)
    hmin = 1e-6
    hmax = 0.99
    lb = (((1 + sqrt(1+8*hmax.^(2)))/4hmax^2)^2 - 1)/16
    ub = (((1 + sqrt(1+8*hmin.^(2)))/4hmin^2)^2 - 1)/16

    optimfunc = Λ -> gcvfun!(Λ, z, s²)
    res = optimize(optimfunc, log10(lb), log10(ub))
    λ = 10^only(Optim.minimizer(res))

    return U*(z./(1 .+ λ*s²))
end

function lcurve_denoise(z, s², U)
    npoints = 500
    lb = -22.
    ub = 22.
    Λ = 10 .^LinRange(lb, ub, npoints)

    ρ = similar(Λ)
    η = similar(Λ)

    for (i, λ) in enumerate(Λ)
        fₖ = @. 1/(1 + λ*s²)
        ρ[i] = mean(abs2, @. z*(fₖ - 1.))
        η[i] = mean(abs2, @. z*√(s²)*fₖ)
    end

    k = curvature(log.(ρ), log.(η), Λ, method = :linear)
    notnan = findall(@. !isnan(k))
    pos_max = argmax(k[notnan])

    return U*(z./(1 .+ Λ[notnan][pos_max]*s²))
end

function denoising(y, alg::KalmanDenoising)

    vest = varest(y)

    R = Diagonal(vest)
    ny = size(vest, 1)

    # Optimization pass
    lb = -10.
    ub = mean(vest)
    objfun = λ -> kf_objfun!(λ, y, R)
    res = optimize(objfun, lb, ub)
    Q = 10^only(Optim.minimizer(res))*I(ny)

    return kalman_denoise(y, Q, R, alg)
end

function kalman_denoise(y, Q, R, alg)

    nd = ndims(y)
    if nd == 1
        y = reshape(y, 1, length(y))
    end

    # Initialization
    nx, nt = size(y)
    xest = similar(y)
    Ppred = similar(y, nx, nx)
    Pest = [similar(Ppred) for _ in 1:nt]

    xest[:, 1] .= y[:, 1]
    Pest[1] .= R

    # Other initializations
    xpred = similar(y, nx)
    ik = similar(xpred)
    Kk = similar(Ppred)
    Sk = similar(Ppred)

    for k = 2:nt
        # Prediction
        xpred .= xest[:, k - 1]
        @. Ppred = Pest[k - 1] + Q

        # Update
        @. ik = y[:, k] - xpred
        @. Sk = Ppred + R
        Kk .= Ppred/Sk
        xest[:, k] .= xpred .+ Kk*ik
        Pest[k] .= (I - Kk)*Ppred*(I - Kk)' + Kk*R*Kk'
    end

    if alg.rts
        # RTS pass
        x = similar(y)
        P = similar(Ppred)

        x[:, end] .= xest[:, end]
        P .= Pest[end]
        for k in reverse(1:nt-1)
            xpred .= xest[:, k]
            Ppred .= Pest[k] + Q

            Kk .= Pest[k]/Ppred
            x[:, k] .= xpred .+ Kk*(x[:, k + 1] - xpred)
            P .= Pest[k] .+ Kk*(P .- Pest[k + 1])*Kk'
        end

        xest .= x
    end

    if nd == 1
        return vec(xest)
    end

    return xest
end

function kf_objfun!(λ, y, R)
    if ndims(y) == 1
        y = reshape(y, 1, length(y))
    end

    # Initialisation
    nx, nt = size(y)
    x = y[:, 1]

    P = similar(y, nx, nx)
    P .= R

    xpred = similar(y, nx)
    Ppred = similar(P)
    ik = similar(xpred)
    Sk = similar(P)
    Kk = similar(P)

    # Application of the initial conditions
    Q = (10 .^λ).*I(nx)

    # Filtering loop
    J = zero(λ)
    for k = 2:nt
        # Prediction
        xpred .= x
        @. Ppred = P + Q

        # Update
        @. ik = y[:, k] - xpred
        @. Sk = Ppred + R
        Kk .= Ppred/Sk
        x .= xpred .+ Kk*ik
        P .= (I - Kk)*Ppred*(I - Kk)' + Kk*R*Kk'

        J += logdet(Sk) + ik'*(Sk\ik)
    end

    return only(J)
end