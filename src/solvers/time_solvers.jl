"""
Structure containing data for the time solver

# Parameters
* K: Stiffness matrix
* M: Mass matrix
* C: Damping matrix
* t: Collection of calculation time steps
"""
@with_kw struct LinearTimeProblem
    K :: Matrix{Float64}
    M :: Matrix{Float64}
    C :: Matrix{Float64}
    F :: Matrix{Float64}
    h

    LinearTimeProblem(K, M, C, F, t) = new(K, M, C, F, t[2] - t[1])
end

"""
Structure containing problem solutions

# Parameters
* D: Displacement matrix
* V: Velocity matrix
* A: Acceleration matrix
"""
@with_kw struct TimeSolution
    D :: Matrix{Float64}
    V :: Matrix{Float64}
    A :: Matrix{Float64}
end

# Solvers
# Central difference scheme
struct CentralDiff end

# Fourth-order Runge-Kutta
struct RK4 end

# Newmark-Family abstract type
abstract type NewmarkFamily end

# Fox-Goodwin
struct FoxGoodwin <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    FoxGoodwin() = new(0., 0., 0.5, 1/12, "Fox-Goodwin...")
end

# Linear Acceleration
struct LinearAcceleration <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    LinearAcceleration() = new(0., 0., 0.5, 1/6, "Linear acceleration...")
end

# Newmark
struct Newmark <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    Newmark(; γ₀ = 0.5, β₀ = 0.25) = new(0., 0., γ₀, β₀, "Newmark...")
end

# Hilber-Hughes-Taylor
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

# Wood-Bossak-Zienkiewicz
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

# Generalized-α
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

# Mid-Point rule
struct MidPoint <: NewmarkFamily
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64
    name::String

    MidPoint() = new(0.5, 0.5, 0.5, 0.25, "Mid-point rule...")
end

# Algorithms

# Central-difference
function solve(prob::LinearTimeProblem, u0, alg::CentralDiff)

    (; K, M, C, F, h) = prob

    nt = size(F, 2)
    Nddl = size(K, 1)

    # Initialization of the result matrices
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Computation of the initial acceleration
    D[:, 1] = u0.D₀
    V[:, 1] = u0.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    LU = lu(M)
    A[:, 1] = LU\rhs0

    D_1 = D[:, 1] - h.*V[:, 1] + (h^2 .*A[:, 1]./4)

    p = Progress(nt - 1, desc = "Central difference...", color = :black, showspeed = true)
    @views @inbounds for n in 1:nt-1
        next!(p)
        if n == 1
            D[:, n+1] = 2D[:, n] - D_1 + (h^2 .*A[:, n])
        else
            D[:, n+1] = 2D[:, n] - D[:, n-1] + (h^2 .*A[:, n])
        end

        V[:, n+1] = (D[:, n+1] - D[:, n])./h

        rhs = F[:, n+1] - C*V[:, n+1] - K*D[:, n+1]
        A[:, n+1] = LU\rhs
    end

    return TimeSolution(D, V, A)
end

# Fourth-order Runge-Kutta
function solve(prob::LinearTimeProblem, u0, alg::RK4)

    (; K, M, C, F, h) = prob

    nt = size(F, 2)
    Nddl = size(K, 1)

    # Initialization of the result matrices
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Computation of the initial acceleration
    D[:, 1] = u0.D₀
    V[:, 1] = u0.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    LU = lu(M)
    A[:, 1] = LU\rhs0

    p = Progress(nt - 1, desc = "RK4...", color=:black, showspeed=true)
    @views @inbounds for n = 1:nt-1
        next!(p)
        Fn_2 = (F[:, n+1] + F[:, n])./2

        k₁ = A[:, n]

        CK = (C + h.*K./2)*V[:, n]
        KD = K*D[:, n]
        k₂ = Fn_2 - CK - h.*(C*k₁)./2 - KD
        k₃ = Fn_2 - CK - h.*(C*k₂)./2 - (h^2 .*(K*k₁)./4) - KD
        k₄ = F[:, n+1] - (C + h.*K)*V[:, n] - h.*(C*k₃) - (h^2 .*(K*k₂)./2) - KD

        D[:, n+1] = D[:, n] + h.*V[:, n] + (h^2 .*(k₁ + (LU\(k₂ + k₃)))./6)
        V[:, n+1] = V[:, n] + h.*(k₁ + (LU\(2k₂ + 2k₃ + k₄)))./6
        A[:, n+1] = LU\(F[:, n+1] - C*V[:, n+1] - K*D[:, n+1])
    end

    return TimeSolution(D, V, A)
end

#Newmark family
function solve(prob::LinearTimeProblem, u0, alg::NewmarkFamily)

    (; K, M, C, F, h) = prob
    (; αf, αₘ, γ₀, β₀, name) = alg

    nt = size(F, 2)
    Nddl = size(K, 1)

    # Initialization of the result matrices
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Computation of the initial acceleration
    D[:, 1] = u0.D₀
    V[:, 1] = u0.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    A[:, 1] = M\rhs0

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
    LU = lu(S)

    p = Progress(nt - 1, desc = name, color = :black, showspeed = true)
    @views @inbounds for n = 1:nt-1
        next!(p)
        rhs = b₈.*F[:, n+1] + b₉.*F[:, n] - C*(b₁.*A[:, n] + V[:, n]) - K*(b₂.*A[:, n] + b₅.*V[:, n] + D[:, n]) - b₇.*(M*A[:, n])

        A[:, n+1] = LU\rhs
        V[:, n+1] = @. V[:, n] + a₁*A[:, n] + a₃*A[:, n+1]
        D[:, n+1] = @. D[:, n] + h*V[:, n] + a₂*A[:, n] + a₄*A[:, n+1]
    end

    return TimeSolution(D, V, A)
end

# Default solver
solve(prob::LinearTimeProblem, u0) = solve(prob, u0, GeneralizedAlpha())