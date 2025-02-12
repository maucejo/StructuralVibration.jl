abstract type DenoisingMethod end

"""
    RegDenoising(prior = :invgamma)

Bayesian Regularization denoising method

# Fields
* `prior`: Prior distribution over the hyperparameters used for computing the regularization parameter
    * `:invgamma`: Inverse Gamma distribution
    * `:uniform`: Uniform distribution
"""
struct RegDenoising <: DenoisingMethod
    prior::Symbol

    RegDenoising(prior = :invgamma) = new(prior)
end

"""
    KalmanDenoising(; rts = false)

Kalman filter denoising method

# Fields
* `rts`: Flag to enable the Rauch-Tung-Striebel smoother
"""
struct KalmanDenoising <: DenoisingMethod
    rts::Bool

    KalmanDenoising(; rts = false) = new(rts)
end

"""
    denoising(y::AbstractArray, alg::RegDenoising)

Denoises a signal `y`

# Inputs
* `y`: Noisy signal
* `alg`: Denoising method
    * `RegDenoising`: Bayesian Regularization denoising method
    * `KalmanDenoising`: Kalman filter denoising method

# Output
* `x`: Denoised signal
"""
function denoising(y::AbstractArray, alg::RegDenoising)
   ndim = ndims(y)
   sy = size(y)
   n = sy[ndim]

   # Initialization
   s = undefs(n)
   z = undefs(eltype(y), sy)
   U = undefs(n, n)

   # Eigenvalues of the second difference matrix
   @. s = 2(1. - cos((0:n-1)π/n))
   @. s[s == 0.] = 1e-8
   s² = s.^2

   # Corresponding eigenvectors
   scale = undefs(n)
   scale[1] = sqrt(1/n)
   scale[2:end] .= sqrt(2/n)
   U = [cos((k - 1)*(2j - 1)π/2n) for j = 1:n, k = 1:n].*scale'

   # Calulation of DCT-2
   z .= dct(y, ndim)

   x = undefs(sy)

    if ndim == 1
        x .= reg_denoise(z, s², U, alg.prior)
    else
        for (idz, zi) in enumerate(eachrow(z))
            x[idz, :] .= reg_denoise(zi, s², U, alg.prior)
        end
    end

    return x
end

function reg_denoise(z, s², U, prior)
    lb = -5.
    ub = 5.

    optimfunc = L -> bayesfun!(L, z, s², prior)
    res = optimize(optimfunc, lb, ub)
    λ = 10^only(Optim.minimizer(res))

    return U*(z./(1 .+ λ*s²))
end

function denoising(y::AbstractArray, alg::KalmanDenoising)

    # Conversion
    y isa Vector ? y = transpose(y) : nothing

    v̂ = varest(y)

    R = Diagonal(v̂)
    ny = size(y, 1)

    # Optimization pass
    lb = -10.
    ub = mean(v̂)
    objfun = λ -> kf_objfun!(λ, y, R)
    res = optimize(objfun, lb, ub)
    Q = 10^only(Optim.minimizer(res))*I(ny)

    return kalman_denoise(y, Q, R, alg)
end

function kalman_denoise(y, Q, R, alg)

    # Initialization
    nx, nt = size(y)
    xest = undefs(eltype(y), nx, nt)
    Pest = [undefs(eltype(y), nx, nx) for _ in 1:nt]

    xest[:, 1] .= y[:, 1]
    Pest[1] .= R

    # Other initializations
    xpred = undefs(eltype(y), nx)
    Ppred = undefs(eltype(y), nx, nx)
    ik = undefs(eltype(y), nx)
    Kk = undefs(eltype(y), nx, nx)
    Sk = undefs(eltype(y), nx, nx)

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
        x = undefs(eltype(y), nx, nt)
        P = undefs(eltype(y), nx, nx)

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

    return nx == 1 ? vec(xest) : xest
end

function kf_objfun!(λ, y, R)
    # Initialisation
    nx, nt = size(y)
    x = y[:, 1]

    P = undefs(eltype(y), nx, nx)
    P .= R

    xpred = undefs(eltype(y), nx)
    Ppred = undefs(eltype(y), nx, nx)
    ik = undefs(eltype(y), nx)
    Sk = undefs(eltype(y), nx, nx)
    Kk = undefs(eltype(y), nx, nx)

    # Application of the initial conditions
    Q = (10 .^λ).*I(nx)

    # Filtering loop
    J = 0.
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