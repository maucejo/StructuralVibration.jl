"""
    DiscreteTimeProblem(K::Matrix{Float64}, M::Matrix{Float64}, C::Matrix{Float64}, u0::Tuple{Vector{Float64}, Vector{Float64}}, h::Float64, F::Matrix{Float64})

Structure containing data for the time solver

# Fields
* K: Stiffness matrix
* M: Mass matrix
* C: Damping matrix
* u0: Initial conditions
    * u0[1]: Initial displacement
    * u0[2]: Initial velocity
* h: Time step
* F: External force matrix
"""
@with_kw struct DiscreteTimeProblem
    K :: Matrix{Float64}
    M :: Matrix{Float64}
    C :: Matrix{Float64}
    u0 :: Tuple{Vector{Float64}, Vector{Float64}}
    h :: Float64
    F :: Matrix{Float64}
end

"""
    DiscreteTimeSolution(D, V, A)

Structure containing problem solutions

# Fields
* D: Displacement matrix
* V: Velocity matrix
* A: Acceleration matrix
"""
@with_kw struct DiscreteTimeSolution
    D :: Matrix{Float64}
    V :: Matrix{Float64}
    A :: Matrix{Float64}
end

# Solvers
"""
    CentralDiff()

Central difference time solver
"""
struct CentralDiff end

"""
    RK4()

Fourth-order Runge-Kutta time solver
"""
struct RK4 end

# Newmark family
abstract type NewmarkFamily end

"""
    FoxGoodwin()

Fox-Goodwin time solver
"""
struct FoxGoodwin <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    FoxGoodwin() = new(0., 0., 0.5, 1/12, "Fox-Goodwin...")
end

"""
    LinearAcceleration()

Linear acceleration time solver
"""
struct LinearAcceleration <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    LinearAcceleration() = new(0., 0., 0.5, 1/6, "Linear acceleration...")
end

"""
    Newmark(; γ₀ = 0.5, β₀ = 0.25)

Newmark time solver
"""
struct Newmark <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    Newmark(; γ₀ = 0.5, β₀ = 0.25) = new(0., 0., γ₀, β₀, "Newmark...")
end

"""
    HHT(; γ₀ = 0.5, β₀ = 0.25, ρ = 1., αf = Inf)

Hilber-Hughes-Taylor time solver
"""
struct HHT <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    function HHT(; γ₀ = 0.5, β₀ = 0.25, ρ = 1., αf = Inf)
        if αf ≠ Inf
            (0. < αf < 1/3) ? error("αf must be in [0, 1/3[") : nothing
            return new(αf, 0., γ₀, β₀, "HHT...")
        else
            ρ < 0.5 ? error("ρ must be in [0.5, 1]") : nothing
            return new((1. - ρ)/(1. + ρ), 0., γ₀, β₀, "HHT...")
        end
    end
end

"""
    WBZ(; γ₀ = 0.5, β₀ = 0.25, ρ = 1., αₘ = Inf)

Wood-Bossak-Zienkiewicz time solver
"""
struct WBZ <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    function WBZ(; γ₀ = 0.5, β₀ = 0.25, ρ = 1., αₘ = Inf)
        if αₘ ≠ Inf
            (αₘ ≤ 0.5) ? error("αₘ must be ≤ 0.5") : nothing
            return new(0., αₘ, γ₀, β₀, "WBZ...")
        else
            (ρ > 1.) ? error("ρ must be in [0, 1]") : nothing
            return new(0., (ρ - 1.)/(ρ + 1.), γ₀, β₀, "WBZ...")
        end
    end
end

"""
    GeneralizedAlpha(; γ₀ = 0.5, β₀ = 0.25, ρ = 1., αf = Inf, αₘ = Inf)

Generalized-α time solver
"""
struct GeneralizedAlpha <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    function GeneralizedAlpha(; γ₀ = 0.5, β₀ = 0.25, ρ = 1., αf = Inf, αₘ = Inf)
        if αf ≠ Inf && αₘ ≠ Inf
            (αₘ ≤ αf ≤ 0.5) ? error("αₘ ≤ αf ≤ 0.5") : nothing
        else
            (ρ > 1.) ? error("ρ must be in [0, 1]") : nothing

            αf = ρ/(ρ + 1.)
            αₘ = 3αf - 1.
        end

        return new(αf, αₘ, γ₀, β₀, "Generalized-α...")
    end
end

"""
    MidPoint()

Mid-point rule time solver
"""
struct MidPoint <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    MidPoint() = new(0.5, 0.5, 0.5, 0.25, "Mid-point rule...")
end

## Algorithms
# Central-difference
function solve(prob::DiscreteTimeProblem, alg::CentralDiff)
    (; K, M, C, u0, h, F) = prob

    Nddl, nt = size(F)
    # Initialization of the result matrices
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Intermediate vectors
    rhs = zeros(Nddl)

    # Computation of the initial acceleration
    D[:, 1] .= u0[1]
    V[:, 1] .= u0[2]

    rhs .= F[:, 1] - C*V[:, 1] - K*D[:, 1]
    lu!(M)
    A[:, 1] = M\rhs

    D_1 = D[:, 1] - h.*V[:, 1] + (h^2 .*A[:, 1]./4)

    p = Progress(nt - 1, desc = "Central difference...", color = :black, showspeed = true)
    @views @inbounds for n in 1:nt-1
        next!(p)
        if n == 1
            D[:, n+1] = 2D[:, n] - D_1 + (h^2*A[:, n])
        else
            D[:, n+1] = 2D[:, n] - D[:, n-1] + (h^2*A[:, n])
        end

        V[:, n+1] = (D[:, n+1] - D[:, n])/h

        rhs .= F[:, n+1] - C*V[:, n+1] - K*D[:, n+1]
        A[:, n+1] = M\rhs
    end

    return DiscreteTimeSolution(D, V, A)
end

# Fourth-order Runge-Kutta
function solve(prob::DiscreteTimeProblem, alg::RK4)
    (; K, M, C, u0, h, F) = prob

    Nddl, nt = size(F)
    # Initialization of the result matrices
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Initialization of the intermediate vectors
    k₁ = zeros(Nddl)
    k₂ = zeros(Nddl)
    k₃ = zeros(Nddl)
    k₄ = zeros(Nddl)
    KD = zeros(Nddl)
    CK = zeros(Nddl)
    Fn_2 = zeros(Nddl)

    # Computation of the initial acceleration
    D[:, 1] .= u0[1]
    V[:, 1] .= u0[2]

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    lu!(M)
    A[:, 1] = M\rhs0

    p = Progress(nt - 1, desc = "RK4...", color=:black, showspeed=true)
    @views @inbounds for n = 1:nt-1
        next!(p)
        Fn_2 .= (F[:, n+1] + F[:, n])/2
        CK .= (C + h.*K./2)*V[:, n]
        KD .= K*D[:, n]

        k₁ .= A[:, n]
        k₂ .= Fn_2 - CK - h*(C*k₁)./2 - KD
        k₃ .= Fn_2 - CK - h*(C*k₂)./2 - (h^2*(K*k₁)/4) - KD
        k₄ .= F[:, n+1] - (C + h*K)*V[:, n] - h*(C*k₃) - (h^2*(K*k₂)/2) - KD

        D[:, n+1] .= D[:, n] + h*V[:, n] + (h^2*(k₁ + (M\(k₂ + k₃)))/6)
        V[:, n+1] .= V[:, n] + h*(k₁ + (M\(2k₂ + 2k₃ + k₄)))/6
        A[:, n+1] .= M\(F[:, n+1] - C*V[:, n+1] - K*D[:, n+1])
    end

    return DiscreteTimeSolution(D, V, A)
end

#Newmark family
function solve(prob::DiscreteTimeProblem, alg::NewmarkFamily)
    (; K, M, C, u0, h, F) = prob
    (; αf, αₘ, γ₀, β₀, name) = alg

    Nddl, nt = size(F)
    # Initialization of the result matrices
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Intermediate vectors
    rhs = zeros(Nddl)

    # Computation of the initial acceleration
    D[:, 1] .= u0[1]
    V[:, 1] .= u0[2]

    rhs .= F[:, 1] - C*V[:, 1] - K*D[:, 1]
    A[:, 1] = M\rhs

    # Newmark scheme parameters
    γ = γ₀*(1. + 2αf - 2αₘ)
    β = β₀*(1. + αf - αₘ)^2

    a₁ = (1. - γ)*h
    a₂ = (0.5 - β)*h^2
    a₃ = γ*h
    a₄ = β*h^2

    b₈ = 1. - αf
    b₁ = b₈*a₁
    b₂ = b₈*a₂
    b₃ = b₈*a₃
    b₄ = b₈*a₄
    b₅ = b₈*h
    b₆ = 1. - αₘ
    b₇ = αₘ
    b₉ = αf

    # Computation of the effective stiffness matrix
    S = @. b₆*M + b₃*C + b₄*K
    lu!(S)

    p = Progress(nt - 1, desc = name, color = :black, showspeed = true)
    @views @inbounds for n = 1:nt-1
        next!(p)
        rhs .= b₈*F[:, n+1] + b₉*F[:, n] - C*(b₁*A[:, n] + V[:, n]) - K*(b₂*A[:, n] + b₅*V[:, n] + D[:, n]) - b₇*M*A[:, n]

        A[:, n+1] = S\rhs
        V[:, n+1] = @. V[:, n] + a₁*A[:, n] + a₃*A[:, n+1]
        D[:, n+1] = @. D[:, n] + h*V[:, n] + a₂*A[:, n] + a₄*A[:, n+1]
    end

    return DiscreteTimeSolution(D, V, A)
end

# Default solver
solve(prob:: DiscreteTimeProblem) = solve(prob, GeneralizedAlpha())