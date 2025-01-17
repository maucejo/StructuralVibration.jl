"""
    ContinuousStateSpace(Ac::Matrix{Float64}, Bc::Matrix{Float64})

Continuous-time state-space model

# Fields
* `Ac`: Continuous-time state matrix A
* `Bc`: Continuous-time input matrix B
"""
@with_kw struct ContinuousStateSpace
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
@with_kw struct DiscreteStateSpace
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    Bdp::Matrix{Float64}

    function DiscreteStateSpace(Ad, Bd, Bdp = Matrix{Float64}[])
        return new(Ad, Bd, Bdp)
    end
end

"""
    c2d(css::ContinuouStateSpace, h::Float64, method::Symbol)

Converts a continuous-time state-space model to a discrete-time state-space model.

# Inputs
* `css`: Continuous-time state-space model
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

    if method == :zoh
        Ad = exp(Ac*h)
        Bd = (Ad .- I)*(Ac\Βc)
        return DiscreteStateSpace(Ad, Bd)
    elseif method == :foh
        Ad = exp(Ac*h)
        Bd = (Ad .- I)*(Ac\Βc)
        Bg = (Ad .- Ac*h .- I)*(Ac^2\Bc)/h
        return DiscreteStateSpace(Ad, Bd .- Bg, Bg)
    elseif method == :blh
        Ad = exp(Ac*h)
        Bd = Ad*Bc*h
        return DiscreteStateSpace(Ad, Bd)
    elseif method == :rk4
        Ad = [24(I + Ac*h) + 12(Ac*h)^2 + 4(Ac*h)^3 + (Ac*h)^4]/24
        Bf = h*[12I + 8Ac*h + 3(Ac*h)^2 + (Ac*h)^3]*Bc/24
        Bg = h*[12I + 4Ac*h + (Ac*h)^2]*Bc/24
        return DiscreteStateSpace(Ad, Bf, Bg)
    end
end

"""
    ss_model(K, M, C)

Generates a continuous-time state-space model from the mass, damping, and stiffness matrices

# Inputs
* `K`: Stiffness matrix
* `M`: Mass matrix
* `C`: Damping matrix

# Output
`css`: ContinuousStateSpace
"""
function ss_model(K, M, C)
    n = size(K, 1)
    lu!(M)
    Ac = [zeros(n, n) I; -M\K -M\C]
    Bc = [zeros(n, n); M\I]

    return ContinuousStateSpace(Ac, Bc)
end

"""
    eigenmode(Ac, Nₘ)

Computes the eigenmodes of a continuous-time state-space model

# Inputs
* `Ac`: Continuous-time state-space matrix
* `Nₘ`: Number of eigenmodes to keep in the modal basis (default: empty)

# Outputs
* `λ`: Eigenvalues
* `Ψ`: Eigenvectors
"""

function eigenmode(Ac::Matrix{Float64}, Nₘ = Int[])
    λ, Ψ = eigen(Ac)

    if length(Nₘ) > 0
        return λ[1:Nₘ], Ψ[:, 1:Nₘ]
    end

    return λ, Ψ
end