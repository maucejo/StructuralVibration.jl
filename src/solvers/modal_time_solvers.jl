"""
    FreeModalTimeProblem(K, M, ξₙ, Nₘ = size(K, 1); ismodal = false)

Structure containing data for the modal time solver

# Fields
* `K`: Stiffness matrix
* `M`: Mass matrix
* `ξₙ`: Damping ratios
* `Nₘ`: Number of modes to retain inf the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data
"""
struct FreeModalTimeProblem
    K :: VecOrMat{Float64}
    M :: Matrix{Float64}
    ξₙ :: Vector{Float64}
    u0 :: Tuple{Vector{Float64}, Vector{Float64}}
    t
    Nₘ :: Int
    ismodal :: Bool

    function FreeModalTimeProblem(K, M, ξₙ, u0, t, Nₘ = size(K, 1), ismodal = false)
        if !isa(ξₙ, Array)
            ξₙ = fill(ξₙ, Nₘ)
        elseif length(ξₙ) != Nₘ
            error("The number of damping ratios must be equal to Nₘ")
        end

        new(K, M, ξₙ, u0, t, Nₘ, ismodal)
    end
end

"""
    HarmonicModalTimeProblem(K, M, ξₙ, u0, t, F, ω = 0., Nₘ = size(K, 1); ismodal = false)

Structure containing data for the modal time solver for computing the forced response due to an harmonic excitation

# Constructor
* `K`: Stiffness matrix (or modal stiffness matrix)
* `M`: Mass matrix (or mass-normalized mode shapes)
* `ξₙ`: Damping ratios
* `u0`: Initial conditions
    * `x₀`: Initial displacement (or modal displacement)
    * `v₀`: Initial velocity (or modal velocity)
* `t`: Time points at which to evaluate the response
* `F`: Amplitude vector (or modal force amplitude vector)
* `freq`: Excitation frequency
* `Nₘ`: Number of modes to retain in the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data

# Fields
* `K`: Stiffness matrix (or modal stiffness matrix)
* `M`: Mass matrix (or mass-normalized mode shapes)
* `ξₙ`: Damping ratios
* `u0`: Initial conditions
    * `x₀`: Initial displacement (or modal displacement)
    * `v₀`: Initial velocity (or modal velocity)
* `t`: Time points at which to evaluate the response
* `F`: Amplitude vector (or modal force amplitude vector)
* `ω`: Excitation angular frequency
* `Nₘ`: Number of modes to retain in the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data
"""
struct HarmonicModalTimeProblem
    K :: Matrix{Float64}
    M :: Matrix{Float64}
    ξₙ :: Vector{Float64}
    u0 :: Tuple{Vector{Float64}, Vector{Float64}}
    t
    F :: Vector{Float64}
    ω :: Float64
    Nₘ :: Int
    ismodal :: Bool

    function HarmonicModalTimeProblem(K, M, ξₙ, u0, t, F, freq, Nₘ = size(K, 1); ismodal = false)
        if !isa(ξₙ, Array)
            ξₙ = fill(ξₙ, Nₘ)
        elseif length(ξₙ) != Nₘ
            error("The number of damping ratios must be equal to Nₘ")
        end

        new(K, M, ξₙ, u0, t, F, 2π*freq, Nₘ, ismodal)
    end
end

"""
    ForcedModalTimeProblem(K, M, ξₙ, u0, t, F, Nₘ = size(K, 1); ismodal = false)

Structure containing data for modal time solver for computing the forced response due to an arbitrary excitation

# Fields
* `K`: Stiffness matrix (or modal stiffness matrix)
* `M`: Mass matrix (or mass-normalized mode shapes)
* `ξₙ`: Damping ratios
* `u0`: Initial conditions
    * `u0[1]`: Initial displacement (or modal displacement)
    * `u0[2]`: Initial velocity (or modal velocity)
* `t`: Time points at which to evaluate the response
* `F`: External force matrix (or modal force matrix)
* `Nₘ`: Number of modes to retain in the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data
"""
struct ForcedModalTimeProblem
    K
    M
    ξₙ :: Vector{Float64}
    u0 :: Tuple{Vector{Float64}, Vector{Float64}}
    t
    F
    Nₘ :: Int
    ismodal :: Bool

    function ForcedModalTimeProblem(K, M, ξₙ, u0, t, F, Nₘ = size(K, 1); ismodal = false)
        if !isa(ξₙ, Array)
            ξₙ = fill(ξₙ, Nₘ)
        elseif length(ξₙ) != Nₘ
            error("The number of damping ratios must be equal to Nₘ")
        end

        new(K, M, ξₙ, u0, t, F, Nₘ, ismodal)
    end
end

"""
    ModalTimeSolution(u, du, ddu)

Structure containing problem solutions

# Fields
* `u`: Displacement matrix or vector
* `du`: Velocity matrix or vector
* `ddu`: Acceleration matrix or vector
"""
struct ModalTimeSolution
    u
    du
    ddu
end

"""
    ModalImpulseSolution(u)

Structure containing the impulse response of a modal system

# Field
* `u`: Impulse response of a Mdof system
"""
struct ModalImpulseSolution
    u :: Union{Matrix{Float64}, Vector{Matrix{Float64}}}
end


"""
    solve(prob::FreeModalTimeProblem)

Compute the free response of a multi-degrees of freedom (Mdof) system using the modal approach.

# Inputs
* `prob`: Structure containing the parameters of the Mdof problem

# Output
* `sol`: ModalTimeSolution structure containing the response of the system at the given time points
"""
function solve(prob::FreeModalTimeProblem)
    (; K, M, ξₙ, u0, t, Nₘ, ismodal) = prob
    x₀, v₀ = u0
    nt = length(t)

    # Modal analysis
    if !ismodal
        λ, Φ = eigen(K, M)
        λₘ = λ[1:Nₘ]
        Φₘ = Φ[:, 1:Nₘ]
        ωₙ = .√λₘ;
        # Note: The mode shapes are mass-normalized, so Mₙ = I

        # Modal initial conditions
        qₓ = Φₘ'*M*x₀
        qᵥ = Φₘ'*M*v₀
    else
        Φₘ = M[:, 1:Nₘ]
        if K isa Vector
            ωₙ = .√K[1:Nₘ]
        else
            ωₙ = .√diag(K)[1:Nₘ]
        end
        qₓ = x₀[1:Nₘ]
        qᵥ = v₀[1:Nₘ]
    end

    # Modal coordinate calculation
    q = undefs(nt, Nₘ)
    dq = undefs(nt, Nₘ)
    ddq = undefs(nt, Nₘ)
    for (m, (ωₘ, ξₘ, Aₘ, qₘ))  in enumerate(zip(ωₙ, ξₙ, qₓ, qᵥ))
        if ωₘ == 0.
            @. q[:, m] = Aₘ + qₘ*t
            @. dq[:, m] = qₘ
            @. ddq[:, m] = 0.
        elseif ξₘ < 1.
            Ωₘ = ωₘ*sqrt(1 - ξₘ^2)
            Ξₘ = ξₘ*ωₘ

            # Constant Bₙ
            Bₘ = (qₘ + Ξₘ*Aₘ)/Ωₘ

            @. q[:, m] = (Aₘ*cos(Ωₘ*t) + Bₘ*sin(Ωₘ*t))*exp(-Ξₘ*t)

            @. dq[:, m] = Ωₘ*(-Aₘ*sin(Ωₘ*t) + Bₘ*cos(Ωₘ*t))*exp(-ξ*ωₘ*t) - Ξₘ*q[:, m]

            @. ddq[:, m] = -2Ξₘ*dq[:, m] - ωₘ^2*q[:, m]

        elseif ξₘ == 1.
            Bₘ = qₘ + ω*Aₘ - ω*ρₘ*sin(ϕₘ)

            @. q[:, m] = (Aₘ + Bₘ*t)*exp(-ωₘ*t)

            @. dq[:, m] = Bₘ*exp(-ωₘ*t) - ωₘ*q[:, m]

            @. ddq[:, m] = -2ωₘ*dq[:, m] - ωₘ^2*q[:, m]

        else
            βₘ = ωₘ*sqrt(ξₘ^2 - 1)
            Bₘ = (qₘ + Ξₘ*Aₘ + ω*ρₘ*sin(ϕₘ))/βₘ

            @. q[:, m] = (Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t))*exp(-Ξₘ*t)

            @. dq[:, m] = βₘ*(Aₘ*sinh(βₘ*t) + Bₘ*cosh(βₘ*t))*exp(-Ξₘ*t) - Ξₘ*q[:, m]

            @. ddq[:, m] = -2Ξₘ*dq[:, m] - ωₘ^2*q[:, m]
        end
    end

    # Computation of the displacement
    u = Φₘ*q';
    du = Φₘ*dq';
    ddu = Φₘ*ddq';

    return ModalTimeSolution(u, du, ddu)
end

"""
    solve(prob::ForcedModalTimeProblem)

Compute the forced response of a multi-degrees of freedom (Mdof) system due to an harmonic excitation using the modal approach.

# Inputs
* `prob`: Structure containing the parameters of the Mdof problem

# Output
* `sol`: ModalTimeSolution structure containing the response of the system at the given time points
"""
function solve(prob::HarmonicModalTimeProblem)
    (; K, M, ξₙ, u0, t, F, ω, Nₘ, ismodal) = prob
    x₀, v₀ = u0
    nt = length(t)

    if size(F, 2) ≠ 1
        error("The external force amplitude must be a vector")
    end

    if !ismodal
        # Modal analysis
        λ, Φ = eigen(K, M)
        λₘ = λ[1:Nₘ]
        Φₘ = Φ[:, 1:Nₘ]
        ωₙ = .√λₘ;
        # Note: The mode shapes are mass-normalized, so Mₙ = I

        # Modal initial conditions
        qₓ = Φₘ'*M*x₀
        qᵥ = Φₘ'*M*v₀
    else
        Φₘ = M[:, 1:Nₘ]
        if K isa Vector
            ωₙ = .√K[1:Nₘ]
        else
            ωₙ = .√diag(K)[1:Nₘ]
        end
        qₓ = x₀[1:Nₘ]
        qᵥ = v₀[1:Nₘ]
    end

    # Modal viscous Damping vector
    Ξ = ξₙ.*ωₙ

    # Modal participation factor
    Lₙ = Φₘ'*F

    # Particular solution
    Qₚ = @. Lₙ/(ωₙ^2 - ω^2 + 2im*Ξ*ω)
    ρ = abs.(Qₚ)
    ϕ = angle.(Qₚ)

    A = @. qₓ - ρ*cos(ϕ)

    # Modal coordinate calculation
    q = undefs(nt, Nₘ)
    dq = undefs(nt, Nₘ)
    ddq = undefs(nt, Nₘ)
    qh = undefs(nt)
    dqh = undefs(nt)
    ddqh = undefs(nt)

    for (m, (ωₘ, ξₘ, Aₘ, qₘ, ρₘ, ϕₘ))  in enumerate(zip(ωₙ, ξₙ, A, qᵥ, ρ, ϕ))
        if ξₘ == 0. && ωₘ == ω
            Aₘ = qₓ[m]
            Bₘ = qₘ/ω
            ρₘ = Lₙ[m]/2ω

            @. q[:, m] = Aₘ*cos(ω*t) + Bₘ*sin(ω*t) + ρₘ*t*sin(ω*t)
            @. dq[:, m] = -Aₘ*ω*sin(ω*t) + Bₘ*ω*cos(ω*t) + ρₘ*(sin(ω*t) + ω*t*cos(ω*t))
            @. ddq[:, m] = -Aₘ*ω^2*cos(ω*t) - Bₘ*ω^2*sin(ω*t) + ρₘ*(2ω*cos(ω*t) - ω^2*t*sin(ω*t))
        else
            if ωₘ == 0.
                @. qh = Aₘ + qₘ*t
                @. dqh = qₘ
                @. ddqh = 0.
            elseif ξₘ < 1.
                Ωₘ = ωₘ*sqrt(1 - ξₘ^2)
                Ξₘ = ξₘ*ωₘ

                # Constant Bₙ
                Bₘ = (qₘ + Ξₘ*Aₘ + ω*ρₘ*sin(ϕₘ))/Ωₘ

                @. qh = (Aₘ*cos(Ωₘ*t) + Bₘ*sin(Ωₘ*t))*exp(-Ξₘ*t)
                @. dqh = Ωₘ*(-Aₘ*sin(Ωₘ*t) + Bₘ*cos(Ωₘ*t))*exp(-Ξₘ*t) - Ξₘ*qh
                @. ddqh = -2Ξₘ*dqh - ωₘ^2*qh
            elseif ξₘ == 1.
                Bₘ = qₘ + ωₘ*Aₘ + ρₘ*ω*sin(ϕₘ)

                @. qh = (Aₘ + Bₘ*t)*exp(-ωₘ*t)
                @. dqh = Bₘ*exp(-ωₘ*t) - ωₘ*qh
                @. ddqh = -2ωₘ*dqh - ωₘ^2*qh
            else
                βₘ = ωₘ*sqrt(ξₘ^2 - 1)
                Bₘ = (qₘ + Ξₘ*Aₘ + ω*ρₘ*sin(ϕₘ))/βₘ

                @. qh = (Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t))*exp(-Ξₘ*t)
                @. dqh = βₘ*(Aₘ*sinh(βₘ*t) + Bₘ*cosh(βₘ*t))*exp(-Ξₘ*t) - Ξₘ*qh
                @. ddqh = -2Ξₘ*dqh - ωₘ^2*qh
            end

            @. q[:, m] = qh + ρₘ*cos(ω*t + ϕₘ)
            @. dq[:, m] = dqh - ρₘ*ω*sin(ω*t + ϕₘ)
            @. ddq[:, m] = ddqh - ρₘ*ω^2*cos(ω*t + ϕₘ)
        end
    end

    u = Φₘ*q';
    du = Φₘ*dq';
    ddu = Φₘ*ddq';

    return ModalTimeSolution(u, du, ddu)
end


"""
    solve(prob::ForcedModalTimeProblem)

Compute the forced response of a multi-degrees of freedom (Mdof) system due to an arbitrary excitation using the modal approach.

# Inputs
* `prob`: Structure containing the parameters of the Mdof problem
* `method`: Method to compute the Duhamel's integral
    * :interp: Interpolation + Gaussian quadrature (default)
    * :conv: Convolution

# Output
* `sol`: ModalTimeSolution structure containing the response of the system at the given time points
"""
function solve(prob::ForcedModalTimeProblem, method = :interp)
    (; K, M, ξₙ, u0, t, F, Nₘ, ismodal) = prob
    x₀, v₀ = u0
    nt = length(t)
    Δt = t[2] - t[1]

    if !ismodal
        # Modal analysis
        λ, Φ = eigen(K, M)
        λₘ = λ[1:Nₘ]
        Φₘ = Φ[:, 1:Nₘ]
        ωₙ = .√λₘ;
        # Note: The mode shapes are mass-normalized, so Mₙ = I

        # Modal initial conditions
        qₓ = Φₘ'*M*x₀
        qᵥ = Φₘ'*M*v₀

        # Modal participation factor
        Lₙ = Φₘ'*F
    else
        # Mode shapes and natural frequencies are already provided
        Φₘ = M[:, 1:Nₘ]
        if K isa Vector
            ωₙ = .√K[1:Nₘ]
        else
            ωₙ = .√diag(K)[1:Nₘ]
        end

        # Modal initial condition
        qₓ = x₀[1:Nₘ]
        qᵥ = v₀[1:Nₘ]

         # Modal participation factor
        Lₙ = F[1:Nₘ]
    end

    # Modal coordinate calculation
    q = undefs(nt, Nₘ)
    qh =undefs(nt)
    h = undefs(nt)
    for (m, (ωₘ, ξₘ, qxₘ, qvₘ)) in enumerate(zip(ωₙ, ξₙ, qₓ, qᵥ))
        if ωₘ == 0.
            @. qh = qxₘ + qvₘ*t
            @. h = t
        elseif ξₘ < 1.
            Ωₘ = ωₘ*sqrt(1 - ξₘ^2)
            Ξₘ = ξₘ*ωₘ
            Aₘ = qxₘ
            Bₘ = (qvₘ + Ξₘ*qxₘ)/Ωₘ
            @. qh = (Aₘ*cos(Ωₘ*t) + Bₘ*sin(Ωₘ*t))*exp(-Ξₘ*t)
            @. h = exp(-Ξₘ*t)*sin(Ωₘ*t)/Ωₘ
        elseif ξₘ == 1.
            Aₘ = qxₘ
            Bₘ = qvₘ + ωₘ*qxₘ
            @. qh = (Aₘ + Bₘ*t)*exp(-ωₘ*t)
            @. h = t*exp(-ωₘ*t)
        else
            βₘ = ωₘ*sqrt(ξₘ^2 - 1)
            Ξₘ = ξₘ*ωₘ
            Aₘ = qxₘ
            Bₘ = (qvₘ + Ξₘ*qxₘ)/βₘ
            @. qh = (Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t))*exp(-Ξₘ*t)
            @. h = exp(-Ξₘ*t)*sinh(βₘ*t)/βₘ
        end

        if method == :interp
            q[:, m] .= qh .+ duhamel_integral(Lₙ[m, :], h, t)
        else
            q[:, m] .= qh .+ Δt*conv(Lₙ[m, :], h)[1:nt]
        end
    end

    u = Φₘ*q';
    du = gradient(u, t);
    ddu = gradient(du, t);

    return ModalTimeSolution(u, du, ddu)
end

"""
    impulse_response(K::Matrix{Float64}, M::Matrix{Float64}, ξₙ, t, Nₘ = size(K, 1); ismat = false)

Compute the impulse response of a multi-degrees of freedom (Mdof) system using the modal approach

# Inputs
* `K`: Stiffness matrix
* `M`: Mass matrix
* `ξₙ`: Damping ratios
* `t`: Time points at which to evaluate the response
* `Nₘ`: Number of modes to retain in the modal basis
* `ismat`: Flag to indicate if the output should be a matrix

# Output
* `sol`: ModalImpulseSolution
"""
function impulse_response(K::Matrix{Float64}, M::AbstractMatrix{Float64}, ξₙ, t, Nₘ = size(K, 1); ismat = false)
    λ, Φ = eigen(K, M)
    fₙ = .√λ[1:Nₘ]/2π
    Φₙ = Φ[:, 1:Nₘ]
    ndofs = size(K, 1)
    nt = length(t)

    if !isa(ξₙ, Array)
        ξₙ = fill(ξₙ, Nₘ)
    elseif length(ξₙ) != Nₘ
        error("The number of damping ratios must be equal to Nₘ")
    end

    h = [zeros(ndofs, ndofs) for _ in 1:nt]
    hsdof = undefs(nt)

    for (m, (fₘ, ξₘ, Φₘ)) in enumerate(zip(fₙ, ξₙ, eachcol(Φₙ)))
        # The modes are mass-normalized
        sdof = Sdof(1., fₘ, ξₘ)
        prob = SdofFreeTimeProblem(sdof, [0., 1.], t)
        hsdof .= solve(prob).u

        for (i, hi) in enumerate(hsdof)
            h[i] .+= Φₘ*hi*Φₘ'
        end
    end

    if ismat
        return ModalImpulseSolution(reshape(reduce(hcat, h), ndofs, ndofs, :))
    end

    return ModalImpulseSolution(h)
end