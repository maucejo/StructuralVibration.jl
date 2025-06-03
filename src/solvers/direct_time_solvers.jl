"""
    DirectTimeProblem(K, M, C, F, u0, t)

Structure containing data for the time solver

**Constructor**
* `K::AbstractMatrix`: Stiffness matrix
* `M::AbstractMatrix`: Mass matrix
* `C:AbstractMatrix`: Damping matrix
* `F::AbstractMatrix`: External force matrix
* `u0::Tuple`: Initial conditions
    * `u0[1]`: Initial displacement (or modal displacement)
    * `u0[2]`: Initial velocity (or modal velocity)
* `t::AbstractVector`: Time points at which to evaluate the response

**Fields**
* `K`: Stiffness matrix
* `M`: Mass matrix
* `C`: Damping matrix
* `F`: External force matrix
* `u0`: Initial conditions
    * `u0[1]`: Initial displacement
    * `u0[2]`: Initial velocity
* `h::Real`: Time step
"""
@show_data struct DirectTimeProblem{Tk <: AbstractMatrix, Tm <: AbstractMatrix, Tc <: AbstractMatrix, Tf <: AbstractMatrix, Tu <: AbstractVector, Th <: Real}
    K::Tk
    M::Tm
    C::Tc
    F::Tf
    u0::Tuple{Tu, Tu}
    h::Th

    function DirectTimeProblem(K::Tk, M::Tm, C::Tc, F::Tf, u0::Tuple{Tu, Tu}, t::AbstractRange) where {Tk, Tm, Tc, Tf, Tu}
        h = t[2] - t[1]
        return new{Tk, Tm, Tc, Tf, Tu, typeof(h)}(K, M, C, F, u0, h)
    end
end

"""
    DirectTimeSolution(u, du, ddu)

Structure containing problem solutions

**Fields**
* `u`: Displacement matrix or vector
* `du`: Velocity matrix or vector
* `ddu`: Acceleration matrix or vector
"""
@show_data struct DirectTimeSolution{T <: Real}
    u::Matrix{T}
    du::Matrix{T}
    ddu::Matrix{T}
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
@show_data struct FoxGoodwin{T <: Real} <: NewmarkFamily
    αf::T
    αm::T
    γ0::T
    β0::T
    name::String

    FoxGoodwin(αf::T = 0., αm::T = 0., γ0::T = 0.5, β0::T = 1/12) where T = new{T}(αf, αm, γ0, β0, "Direct Time Problem - Fox-Goodwin...")
end

"""
    LinearAcceleration()

Linear acceleration time solver
"""
@show_data struct LinearAcceleration{T <: Real} <: NewmarkFamily
    αf::T
    αm::T
    γ0::T
    β0::T
    name::String

    LinearAcceleration(αf::T = 0., αm::T = 0., γ0::T = 0.5, β0::T = 1/6) where T = new{T}(αf, αm, γ0, β0, "Direct Time Problem - Linear acceleration...")
end

"""
    Newmark(; γ0 = 0.5, β0 = 0.25)

Newmark time solver
"""
@show_data struct Newmark{T <: Real} <: NewmarkFamily
    αf::T
    αm::T
    γ0::T
    β0::T
    name::String

    Newmark(; γ0::T = 0.5, β0::T = 0.25) where T = new{T}(0., 0., γ0, β0, "Direct Time Problem - Newmark...")
end

"""
    HHT(; γ0 = 0.5, β0 = 0.25, ρ = 1., αf = Inf)

Hilber-Hughes-Taylor time solver
"""
@show_data struct HHT{T <: Real} <: NewmarkFamily
    αf::T
    αm::T
    γ0::T
    β0::T
    name::String

    function HHT(; γ0::T = 0.5, β0::T = 0.25, ρ::T = 1., αf::T = Inf) where T
        if αf ≠ Inf
            (0. < αf < 1/3) ? error("αf must be in [0, 1/3[") : nothing
            return new{T}(αf, 0., γ0, β0, "Direct Time Problem - HHT...")
        else
            ρ < 0.5 ? error("ρ must be in [0.5, 1]") : nothing
            return new{T}((1. - ρ)/(1. + ρ), 0., γ0, β0, "Direct Time Problem - HHT...")
        end
    end
end

"""
    WBZ(; γ0 = 0.5, β0 = 0.25, ρ = 1., αm = Inf)

Wood-Bossak-Zienkiewicz time solver
"""
@show_data struct WBZ{T <: Real} <: NewmarkFamily
    αf::T
    αm::T
    γ0::T
    β0::T
    name::String

    function WBZ(; γ0::T = 0.5, β0::T = 0.25, ρ::T = 1., αm::T = Inf) where T
        if αm ≠ Inf
            (αm ≤ 0.5) ? error("αm must be ≤ 0.5") : nothing
            return new{T}(0., αm, γ0, β0, "Direct Time Problem - WBZ...")
        else
            (ρ > 1.) ? error("ρ must be in [0, 1]") : nothing
            return new{T}(0., (ρ - 1.)/(ρ + 1.), γ0, β0, "Direct Time Problem - WBZ...")
        end
    end
end

"""
    GeneralizedAlpha(; γ0 = 0.5, β0 = 0.25, ρ = 1., αf = Inf, αm = Inf)

Generalized-α time solver
"""
@show_data struct GeneralizedAlpha{T <: Real} <: NewmarkFamily
    αf::T
    αm::T
    γ0::T
    β0::T
    name::String

    function GeneralizedAlpha(; γ0::T = 0.5, β0::T = 0.25, ρ::T = 1., αf::T = Inf, αm = Inf) where T
        if αf ≠ Inf && αm ≠ Inf
            (αm ≤ αf ≤ 0.5) ? error("αm ≤ αf ≤ 0.5") : nothing
        else
            (ρ > 1.) ? error("ρ must be in [0, 1]") : nothing

            αf = ρ/(ρ + 1.)
            αm = 3αf - 1.
        end

        return new{T}(αf, αm, γ0, β0, "Direct Time Problem - Generalized-α...")
    end
end

"""
    MidPoint()

Mid-point rule time solver
"""
@show_data struct MidPoint{T <: Real} <: NewmarkFamily
    αf::T
    αm::T
    γ0::T
    β0::T
    name::String

    MidPoint(αf::T = 0.5, αm::T = 0.5, γ0::T = 0.5, β0::T = 0.25) where T = new{T}(αf, αm, γ0, β0, "Direct Time Problem - Mid-point rule...")
end

## Algorithms
"""
    solve(prob::DirectTimeProblem, alg; progress = true)

Direct time integration solver

**Inputs**
* `prob`: DirectTimeProblem structure
* `alg`: Numerical integration algorithm
    * `CentralDiff()`: Central difference
    * `RK4()`: Fourth-order Runge-Kutta
    * `Newmark()`: Newmark
    * `HHT()`: Hilber-Hughes-Taylor
    * `WBZ()`: Wood-Bossak-Zienkiewicz
    * `GeneralizedAlpha()`: Generalized-α (default)
    * `MidPoint()`: Mid-point rule
* `progress`: Show progress bar

**Output**
* `sol`: Solution structure containing the response of the system at the given time points
    * `u`: Displacement
    * `du`: Velocity
    * `ddu`: Acceleration
"""
function solve(prob::DirectTimeProblem, alg::CentralDiff; progress = true)
    (; K, M, C, F, u0, h) = prob

    Nddl, nt = size(F)
    # Initialization of the result matrices
    D = similar(K, Nddl, nt)
    V = similar(D)
    A = similar(D)

    # Intermediate vectors
    rhs = similar(F, Nddl)

    # Computation of the initial acceleration
    D[:, 1] .= u0[1]
    V[:, 1] .= u0[2]

    rhs .= F[:, 1] - C*V[:, 1] - K*D[:, 1]
    lu!(M)
    A[:, 1] = M\rhs

    D_1 = D[:, 1] - h.*V[:, 1] + (h^2 .*A[:, 1]./4)

    p = Progress(nt - 1, desc = "Direct Time Problem - Central difference...", showspeed = true)
    @views @inbounds for n in 1:nt-1
        progress ? next!(p) : nothing

        if n == 1
            D[:, n+1] = 2D[:, n] - D_1 + (h^2*A[:, n])
        else
            D[:, n+1] = 2D[:, n] - D[:, n-1] + (h^2*A[:, n])
        end

        V[:, n+1] = (D[:, n+1] - D[:, n])/h

        rhs .= F[:, n+1] - C*V[:, n+1] - K*D[:, n+1]
        A[:, n+1] = M\rhs
    end

    return DirectTimeSolution(D, V, A)
end

# Fourth-order Runge-Kutta
function solve(prob::DirectTimeProblem, alg::RK4; progress = true)
    (; K, M, C, F, u0, h) = prob

    Nddl, nt = size(F)
    # Initialization of the result matrices
    D = similar(K, Nddl, nt)
    V = similar(D)
    A = similar(D)

    # Initialization of the intermediate vectors
    k₁ = similar(F, Nddl)
    k₂ = similar(k₁)
    k₃ = similar(k₁)
    k₄ = similar(k₁)
    KD = similar(k₁)
    CK = similar(k₁)
    Fn_2 = similar(k₁)

    # Computation of the initial acceleration
    D[:, 1] .= u0[1]
    V[:, 1] .= u0[2]

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    lu!(M)
    A[:, 1] = M\rhs0

    p = Progress(nt - 1, desc = "Direct Time Problem - RK4...", showspeed = true)
    @views @inbounds for n = 1:nt-1
        progress ? next!(p) : nothing

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

    return DirectTimeSolution(D, V, A)
end

#Newmark family
function solve(prob::DirectTimeProblem, alg::NewmarkFamily; progress = true)
    (; K, M, C, F, u0, h) = prob
    (; αf, αm, γ0, β0, name) = alg

    Nddl, nt = size(F)
    # Initialization of the result matrices
    D = similar(K, Nddl, nt)
    V = similar(D)
    A = similar(D)

    # Intermediate vectors
    rhs = similar(F, Nddl)

    # Computation of the initial acceleration
    D[:, 1] .= u0[1]
    V[:, 1] .= u0[2]

    rhs .= F[:, 1] - C*V[:, 1] - K*D[:, 1]
    A[:, 1] = M\rhs

    # Newmark scheme parameters
    γ = γ0*(1. + 2αf - 2αm)
    β = β0*(1. + αf - αm)^2

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
    b₆ = 1. - αm
    b₇ = αm
    b₉ = αf

    # Computation of the effective stiffness matrix
    S = @. b₆*M + b₃*C + b₄*K
    lu!(S)

    p = Progress(nt - 1, desc = name, showspeed = true)
    @views @inbounds for n = 1:nt-1
        progress ? next!(p) : nothing

        rhs .= b₈*F[:, n+1] + b₉*F[:, n] - C*(b₁*A[:, n] + V[:, n]) - K*(b₂*A[:, n] + b₅*V[:, n] + D[:, n]) - b₇*M*A[:, n]

        A[:, n+1] = S\rhs
        V[:, n+1] = @. V[:, n] + a₁*A[:, n] + a₃*A[:, n+1]
        D[:, n+1] = @. D[:, n] + h*V[:, n] + a₂*A[:, n] + a₄*A[:, n+1]
    end

    return DirectTimeSolution(D, V, A)
end

# Default solver
solve(prob:: DirectTimeProblem; progress = true) = solve(prob, GeneralizedAlpha(), progress = progress)