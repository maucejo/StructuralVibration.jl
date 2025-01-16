"""
    ModalTimeProblem(K, M, ξₙ, Nₘ = size(K, 1); ismodal = false)

Structure containing data for the modal time solver

# Fields
* `K`: Stiffness matrix
* `M`: Mass matrix
* `ξₙ`: Damping ratios
* `Nₘ`: Number of modes to retain inf the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data
"""
@with_kw struct FreeModalTimeProblem
    K :: Matrix{Float64}
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

# Fields
 * `K`: Stiffness matrix (or modal stiffness matrix)
    * `M`: Mass matrix (or mass-normalized mode shapes)
    * `ξₙ`: Damping ratios
    * `u0`: Initial conditions
        * `x₀`: Initial displacement (or modal displacement)
        * `v₀`: Initial velocity (or modal velocity)
* `t`: Time points at which to evaluate the response
* `F`: Amplitude vector (or modal force amplitude vector)
* `ω`: Excitation frequency (for harmonic excitation)
* `Nₘ`: Number of modes to retain in the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data
"""
@with_kw struct HarmonicModalTimeProblem
    K :: Matrix{Float64}
    M :: Matrix{Float64}
    ξₙ :: Vector{Float64}
    u0 :: Tuple{Vector{Float64}, Vector{Float64}}
    t
    F :: Vector{Float64}
    ω :: Float64
    Nₘ :: Int
    ismodal :: Bool

    function HarmonicModalTimeProblem(K, M, ξₙ, u0, t, F, ω, Nₘ = size(K, 1); ismodal = false)
        if !isa(ξₙ, Array)
            ξₙ = fill(ξₙ, Nₘ)
        elseif length(ξₙ) != Nₘ
            error("The number of damping ratios must be equal to Nₘ")
        end

        new(K, M, ξₙ, u0, t, F, ω, Nₘ, ismodal)
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
@with_kw struct ForcedModalTimeProblem
    K :: Matrix{Float64}
    M :: Matrix{Float64}
    ξₙ :: Vector{Float64}
    u0 :: Tuple{Vector{Float64}, Vector{Float64}}
    t
    F :: Matrix{Float64}
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
    ModalTimeSolution(D, V, A)

Structure containing problem solutions

# Fields
* `D`: Displacement matrix
* `V`: Velocity matrix
* `A`: Acceleration matrix
"""
@with_kw struct ModalTimeSolution
    D :: Matrix{Float64}
    V :: Matrix{Float64}
    A :: Matrix{Float64}
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
        Φₘ = M
        ωₙ = .√(diag(K))
        qₓ = x₀
        qᵥ = v₀
    end

    # Modal coordinate calculation
    q = typeof(M)(undef, nt, Nₘ)
    dq = typeof(M)(undef, nt, Nₘ)
    ddq = typeof(M)(undef, nt, Nₘ)
    for (m, (ωₘ, ξₘ, Aₘ, qₘ))  in enumerate(zip(ωₙ, ξₙ, qₓ, qᵥ))
        if ξₘ < 1.
            Ωₘ = ωₘ*sqrt(1 - ξₘ^2)
            Ξₘ = ξₘ*ωₘ

            # Constant Bₙ
            Bₘ = (qₘ + Ξₘ*Aₘ)/Ωₘ

            @. q[:, m] = (Aₘ*cos(Ωₘ*t) + Bₘ*sin(Ωₘ*t))*exp(-Ξₘ*t)

            @. dq[:, m] = (-Aₘ*Ξₘ*cos(Ωₘ*t) - Bₘ*Ξₘ*sin(Ωₘ*t) - Aₘ*Ωₘ*sin(Ωₘ*t) + Bₘ*Ωₘ*cos(Ωₘ*t))*exp(-Ξₘ*t)

            @. ddq[:, m] = (Aₘ*Ξₘ^2*cos(Ωₘ*t) + Bₘ*Ξₘ^2*sin(Ωₘ*t) + 2*Aₘ*Ξₘ*Ωₘ*sin(Ωₘ*t) - 2*Bₘ*Ξₘ*Ωₘ*cos(Ωₘ*t) - Aₘ*Ωₘ^2*cos(Ωₘ*t) - Bₘ*Ωₘ^2*sin(Ωₘ*t))*exp(-Ξₘ*t)

        elseif ξₘ == 1.
            Bₘ = qₘ + ω*Aₘ - ω*ρₘ*sin(ϕₘ)

            @. q[:, m] = (Aₘ + Bₘ*t)*exp(-ωₘ*t)

            @. dq[:, m] = (Bₘ - ωₘ*(Aₘ + Bₘ*t))*exp(-ωₘ*t)

            @. ddq[:, m] = (ωₘ^2*(Aₘ + Bₘ*t) - 2Bₘ*ωₘ)*exp(-ωₘ*t)

        else
            βₘ = ωₘ*sqrt(ξₘ^2 - 1)
            Bₘ = (qₘ + Ξₘ*Aₘ + ω*ρₘ*sin(ϕₘ))/βₘ

            @. q[:, m] = (Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t))*exp(-Ξₘ*t)

            @. dq[:, m] = (Aₘ*βₘ*sinh(βₘ*t) + Bₘ*βₘ*cosh(βₘ*t) - Ξₘ*(Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t)))*exp(-Ξₘ*t)

            @. ddq[:, m] = (Aₘ*βₘ^2*cosh(βₘ*t) + Bₘ*βₘ^2*sinh(βₘ*t) - 2*Ξₘ*(Aₘ*βₘ*sinh(βₘ*t) + Bₘ*βₘ*cosh(βₘ*t)) + Ξₘ^2*(Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t)))*exp(-Ξₘ*t)
        end
    end

    # Computation of the displacement
    D = Φₘ*q';
    V = Φₘ*dq';
    A = Φₘ*ddq';

    return ModalTimeSolution(D, V, A)
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
        Φₘ = M
        ωₙ = .√(diag(K))
        qₓ = x₀
        qᵥ = v₀
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
    q = typeof(M)(undef, nt, Nₘ)
    dq = typeof(M)(undef, nt, Nₘ)
    ddq = typeof(M)(undef, nt, Nₘ)
    for (m, (ωₘ, ξₘ, Aₘ, qₘ, ρₘ, ϕₘ))  in enumerate(zip(ωₙ, ξₙ, A, qᵥ, ρ, ϕ))
        if ξₘ < 1.
            Ωₘ = ωₘ*sqrt(1 - ξₘ^2)
            Ξₘ = ξₘ*ωₘ

            # Constant Bₙ
            Bₘ = (qₘ + Ξₘ*Aₘ + ω*ρₘ*sin(ϕₘ))/Ωₘ

            @. q[:, m] = (Aₘ*cos(Ωₘ*t) + Bₘ*sin(Ωₘ*t))*exp(-Ξₘ*t) + ρₘ*cos(ω*t + ϕₘ)

            @. dq[:, m] = (-Aₘ*Ξₘ*cos(Ωₘ*t) - Bₘ*Ξₘ*sin(Ωₘ*t) - Aₘ*Ωₘ*sin(Ωₘ*t) + Bₘ*Ωₘ*cos(Ωₘ*t))*exp(-Ξₘ*t) - ρₘ*ω*sin(ω*t + ϕₘ)

            @. ddq[:, m] = (Aₘ*Ξₘ^2*cos(Ωₘ*t) + Bₘ*Ξₘ^2*sin(Ωₘ*t) + 2*Aₘ*Ξₘ*Ωₘ*sin(Ωₘ*t) - 2*Bₘ*Ξₘ*Ωₘ*cos(Ωₘ*t) - Aₘ*Ωₘ^2*cos(Ωₘ*t) - Bₘ*Ωₘ^2*sin(Ωₘ*t))*exp(-Ξₘ*t) - ρₘ*ω^2*cos(ω*t + ϕₘ)

        elseif ξₘ == 1.
            Bₘ = qₘ + ω*Aₘ - ω*ρₘ*sin(ϕₘ)

            @. q[:, m] = (Aₘ + Bₘ*t)*exp(-ωₘ*t) + ρₘ*cos(ω*t + ϕₘ)

            @. dq[:, m] = (Bₘ - ωₘ*(Aₘ + Bₘ*t))*exp(-ωₘ*t) - ρₘ*ω*sin(ω*t + ϕₘ)

            @. ddq[:, m] = (ωₘ^2*(Aₘ + Bₘ*t) - 2Bₘ*ωₘ)*exp(-ωₘ*t) - ρₘ*ω^2*cos(ω*t + ϕₘ)

        else
            βₘ = ωₘ*sqrt(ξₘ^2 - 1)
            Bₘ = (qₘ + Ξₘ*Aₘ + ω*ρₘ*sin(ϕₘ))/βₘ

            @. q[:, m] = (Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t))*exp(-Ξₘ*t) + ρₘ*cos(ω*t + ϕₘ)

            @. dq[:, m] = (Aₘ*βₘ*sinh(βₘ*t) + Bₘ*βₘ*cosh(βₘ*t) - Ξₘ*(Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t)))*exp(-Ξₘ*t) - ρₘ*ω*sin(ω*t + ϕₘ)

            @. ddq[:, m] = (Aₘ*βₘ^2*cosh(βₘ*t) + Bₘ*βₘ^2*sinh(βₘ*t) - 2*Ξₘ*(Aₘ*βₘ*sinh(βₘ*t) + Bₘ*βₘ*cosh(βₘ*t)) + Ξₘ^2*(Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t)))*exp(-Ξₘ*t) - ρₘ*ω^2*cos(ω*t + ϕₘ)
        end
    end

    D = Φₘ*q';
    V = Φₘ*dq';
    A = Φₘ*ddq';

    return ModalTimeSolution(D, V, A)
end


"""
    solve(prob::ModalTimeProblem)

Compute the forced response of a multi-degrees of freedom (Mdof) system due to an arbitrary excitation using the modal approach.

# Inputs
* `prob`: Structure containing the parameters of the Mdof problem

# Output
* `sol`: ModalTimeSolution structure containing the response of the system at the given time points
"""
function solve(prob::ForcedModalTimeProblem)
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
        Φₘ = M
        ωₙ = .√(diag(K))

        # Initial modal displacement and velocity
        qₓ = x₀
        qᵥ = v₀

        # Modal participation factor
        Lₙ = F
    end

    # Modal coordinate calculation
    q = zeros(nt, Nₘ)
    qₕ = zeros(nt)
    h = zeros(nt)
    for (m, (ωₘ, ξₘ, Aₘ, Bₘ)) in enumerate(zip(ωₙ, ξₙ, Aₙ, Bₙ))
        if ξₘ < 1.
            Ωₘ = ωₘ*sqrt(1 - ξₘ^2)
            Ξₘ = ξₘ*ωₘ
            Aₘ = qₓ
            Bₘ = (qᵥ + Ξₘ*qₓ)/Ωₘ
            @. qₕ = (Aₘ*cos(Ωₘ*t) + Bₘ*sin(Ωₘ*t))*exp(-Ξₘ*t)
            @. h = exp(-Ξₘ*t)*sin(Ωₘ*t)/Ωₘ
        elseif ξₘ == 1.
            Aₘ = qₓ
            Bₘ = qᵥ + ωₘ*qₓ
            @. qₕ = (Aₘ + Bₘ*t)*exp(-ωₘ*t)
            @. h = t*exp(-ωₘ*t)
        else
            βₘ = ωₘ*sqrt(ξₘ^2 - 1)
            Ξₘ = ξₘ*ωₘ
            Aₘ = qₓ
            Bₘ = (qᵥ + Ξₘ*qₓ)/βₘ
            @. qₕ = (Aₘ*cosh(βₘ*t) + Bₘ*sinh(βₘ*t))*exp(-Ξₘ*t)
            @. h = exp(-Ξₘ*t)*sinh(βₘ*t)/βₘ
        end

        q[:, m] .= qₕ .+ Δt*conv(Lₙ[n, :], h)[1:nt]
    end

    D = Φₘ*q';
    V = gradient(D, t);
    A = gradient(V, t);

    return ModalTimeSolution(D, V, A)
end