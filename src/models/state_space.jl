"""
    ContinuousStateSpace(Ac::Matrix{Float64}, Bc::Matrix{Float64})

Continuous-time state-space model

# Fields
* `Ac`: Continuous-time state matrix A
* `Bc`: Continuous-time input matrix B
"""
struct ContinuousStateSpace
    Ac::Matrix{Float64}
    Bc::Matrix{Float64}
end

"""
    DiscreteStateSpace(Ad::Matrix{Float64}, Bd::Matrix{Float64}, Bdp::Matrix{Float64})

Discrete-time state-space model

# Fields
* `Ad`: Discrete-time state matrix A
* `Bd`: Discrete-time input matrix B
* `Bdp`: Discrete-time input matrix Bp (only for `:foh` method)
"""
struct DiscreteStateSpace
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    Bdp::Matrix{Float64}

    function DiscreteStateSpace(Ad, Bd, Bdp = Matrix{Float64}[])
        return new(Ad, Bd, Bdp)
    end
end

"""
    StateSpaceProblem(K::Matrix{Float64}, M::Matrix{Float64}, C::Matrix{Float64})

Structure containing data for the state-space model

# Constructor parameters
* `K`: Stiffness matrix
* `M`: Mass matrix
* `C`: Damping matrix
* `u0`: Initial conditions
* `h`: Time step
* `F`: External force matrix

# Fields
* css: ContinuousStateSpace
* u0: Initial conditions
* h: Time step
* F: External force matrix
"""
struct StateSpaceProblem
    css :: ContinuousStateSpace
    u0 :: Vector{Float64}
    h :: Float64
    F :: Matrix{Float64}

    StateSpaceProblem(K, M, C, u0, h, F) = new(ss_model(K, M, C), u0, h, F)
end

struct StateSpaceSolution
    D::Matrix{Float64}
    V::Matrix{Float64}
    A::Matrix{Float64}
end

"""
    c2d(Ac::Matrix{Float64}, Bc::Matrix{Float64}, h::Float64, method::Symbol)

Converts a continuous-time state-space model to a discrete-time state-space model.

# Inputs
* `Ac`: Continuous-time state-space matrix A.
* `Bc::Matrix{Float64}`: Continuous-time state-space matrix B.
* `h`: Sampling time.
* `method::Symbol`: Discretization method
    * `:zoh`: Zero-order Hold method
    * `:foh`: First-order Hold method
    * `:blh`: Band-limited Hold method

# Output
* `DiscreteStateSpace`: Discrete-time state-space model
    * `Ad`: Discrete-time state-space matrix A
    * `Bd`: Discrete-time state-space matrix B
    * `Bdp`: Discrete-time state-space matrix Bp (only for `:foh` method)
"""
function c2d(css::ContinuousStateSpace, h, method)
    (; Ac, Bc) = css
    Ad = exp(Ac*h)
    if method == :zoh
        Bd = (A .- I)*(A\Β)
        return DiscreteStateSpace(Ad, Bd)
    elseif method == :foh
        Bd = (A .- I)*(A\Β)
        Bg = (Ad .- Ac*h .- I)*(Ac^2\Bc)/h
        return DiscreteStateSpace(Ad, Bd .- Bg, Bg)
    elseif method == :blh
        Bd = A*Bc*h
        return DiscreteStateSpace(Ad, Bd)
    end
end

"""
    ss_model(K::Matrix{Float64}, M::Matrix{Float64}, C::Matrix{Float64})

Generates a continuous-time state-space model from the mass, damping, and stiffness matrices

# Inputs
* `K`: Stiffness matrix
* `M`: Mass matrix
* `C`: Damping matrix

# Output
`css`: ContinuousStateSpace
"""
function ss_model(K::Matrix{Float64}, M::Matrix{Float64}, C::Matrix{Float64})
    n = size(K, 1)
    lu!(M)
    Ac = [zeros(n) I; -M\K -M\C]
    Bc = [zeros(n); M\I]

    return ContinuousStateSpace(Ac, Bc)
end

"""
    eigenmode(css::ContinuousStateSpace, Nₘ)

Computes the eigenmodes of a continuous-time state-space model

# Inputs
* `css`: Continuous-time state-space model
* `Nₘ`: Number of eigenmodes to keep in the modal basis

# Outputs
* `λ`: Eigenvalues
* `Ψ`: Eigenvectors
"""

function eigenmode(css::ContinuousStateSpace, Nₘ)
    (; Ac) = css
    λ, Ψ = eigen(Ac)

    return λ[1:Nₘ], Ψ[:, 1:Nₘ]
end

"""
    solve(prob::DiscreteTimeProblem, u0::Vector{Float64}, h::Float64, F::Matrix{Float64}, method = :zoh)

Solves a discrete-time problem using the state-space model

# Inputs
* `prob`: Discrete-time problem
* `u0`: Vector of initial conditions
* `h`: Sampling time
* `F`: Matrix of external forces
* `method::Symbol`: Discretization method
    * `:zoh`: Zero-order Hold method
    * `:foh`: First-order Hold method
    * `:blh`: Band-limited Hold method

# Output
* `StateSpaceSolution`: Solution of the state-space model
"""
function solve(prob::StateSpaceProblem, method = :zoh)
    (; css) = prob
    dss = c2d(css, h, method)

    n, nt = size(F)
    m = Int(n/2)
    x = zeros(n, nt)
    A = zeros(m, nt)

    x[:, 1] .= u0[:]
    A[:, 1] .= css.Bd[m+1:end]*F[:, 1] .+ css.Ac[m+1:end, :]*x[:, 1]

    @views @inbounds for k in 1:nt-1
        if method == :zoh && method == :blh
            x[:, k+1] .= dss.Ad*x[:, k] .+ dss.Bd*F[:, k]
        elseif method == :foh
            x[:, k+1] .= dss.Ad*x[:, k] .+ dss.Bd*F[:, k] .+ dss.Bdp*F[:, k+1]
        end

        A[:, k+1] .= css.Bc[m+1:end]*F[:, k+1] .+ css.Ac[m+1:end, :]*x[:, k+1]
    end

    return StateSpaceSolution(x[1:m, :], x[m+1:end, :], A)
end


"""
    solve(prob::StateSpaceProblem, alg = RK4())

Solves a continuous-time problem using the state-space model

# Inputs
* `prob`: Continuous-time problem
* `alg`: Time integration algorithm

# Output
* `StateSpaceSolution`: Solution of the state-space model
"""
function solve(prob::StateSpaceProblem, alg = RK4())
    (; css) = prob
    (; Ac, Bc) = css

    n, nt = size(F)
    m = Int(n/2)
    x = zeros(n, nt)
    A = zeros(m, nt)

    # Intermediate vectors
    k₁ = zeros(n)
    k₂ = zeros(n)
    k₃ = zeros(n)
    k₄ = zeros(n)
    Fn_2 = zeros(n)

    x[:, 1] .= u0[:]

    @views @inbounds for k in 1:nt-1
        Fn_2 .= (F[:, k] .+ F[:, k+1])/2

        k₁ .= Ac*x[:, k] .+ Bc*F[:, k]
        k₂ .= Ac*(x[:, k] .+ h*k₁/2) .+ Bc*Fn_2
        k₃ .= Ac*(x[:, k] .+ h*k₂/2) .+ Bc*Fn_2
        k₄ .= Ac*(x[:, k] .+ h*k₃) .+ Bc*F[:, k+1]

        @. x[:, k+1] = x[:, k] + h*(k₁ + 2k₂ + 2k₃ + k₄)/6
        A[:, k+1] .= css.Bc[m+1:end]*F[:, k+1] .+ css.Ac[m+1:end, :]*x[:, k+1]
    end

    return StateSpaceSolution(x[1:m, :], x[m+1:end, :], A)
end