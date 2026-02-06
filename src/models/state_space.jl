"""
    ContinuousStateSpace(Ac, Bc)

Continuous-time state-space model

**Fields**
* `Ac`: Continuous-time state matrix A
* `Bc`: Continuous-time input matrix B
"""
@show_data struct ContinuousStateSpace{T <: Real}
    Ac::Matrix{T}
    Bc::Matrix{T}
end

"""
    DiscreteStateSpace(Ad, Bd, Bdp)

Discrete-time state-space model

**Fields**
* `Ad`: Discrete-time state matrix A
* `Bd`: Discrete-time input matrix B
* `Bdp`: Discrete-time input matrix Bp (only for `:foh` method)
"""
@show_data struct DiscreteStateSpace{T <: Real}
    Ad::Matrix{T}
    Bd::Matrix{T}
    Bdp::Matrix{T}

    DiscreteStateSpace(Ad::Matrix{T}, Bd::Matrix{T}, Bdp::Matrix{T} = similar(Bd)) where T = new{T}(Ad, Bd, Bdp)
end

"""
    c2d(css::ContinuouStateSpace, h::Float64, method::Symbol)

Converts a continuous-time state-space model to a discrete-time state-space model.

**Inputs**
* `css`: Continuous-time state-space model
* `h`: Sampling time.
* `method`: Discretization method
    * `:zoh`: Zero-order Hold method
    * `:foh`: First-order Hold method
    * `:blh`: Band-limited Hold method
    * `:rk4`: 4th order Runge-Kutta method

**Output**
* `DiscreteStateSpace`: Discrete-time state-space model
    * `Ad`: Discrete-time state-space matrix A
    * `Bd`: Discrete-time state-space matrix B
    * `Bdp`: Discrete-time state-space matrix Bp (only for `:foh` and `:rk4` methods)
"""
function c2d(css::ContinuousStateSpace, h, method = :zoh)::DiscreteStateSpace
    (; Ac, Bc) = css

    if method == :zoh
        Ad = exp(Ac*h)
        Bd = (Ad - I)*(Ac\Bc)
        return DiscreteStateSpace(Ad, Bd)
    elseif method == :foh
        Ad = exp(Ac*h)
        Bd = (Ad - I)*(Ac\Bc)
        Bg = (Ad .- Ac*h - I)*(Ac^2\Bc)/h
        return DiscreteStateSpace(Ad, Bd .- Bg, Bg)
    elseif method == :blh
        Ad = exp(Ac*h)
        Bd = Ad*Bc*h
        return DiscreteStateSpace(Ad, Bd)
    elseif method == :rk4
        Ad = (24(I + Ac*h) + 12(Ac*h)^2 + 4(Ac*h)^3 + (Ac*h)^4)/24
        Bf = h*(12I + 8Ac*h + 3(Ac*h)^2 + (Ac*h)^3)*Bc/24
        Bg = h*(12I + 4Ac*h + (Ac*h)^2)*Bc/24
        return DiscreteStateSpace(Ad, Bf, Bg)
    end
end

"""
    ss_model(K, M, C)

Generates a continuous-time state-space model from the mass, damping, and stiffness matrices

**Inputs**
* `K`: Stiffness matrix
* `M`: Mass matrix
* `C`: Damping matrix

**Output**
* `css`: ContinuousStateSpace
"""
function ss_model(K, M, C)

    n = size(K, 1)
    M_fact = lu(M)
    Ac = [zeros(n, n) I; -M\K -M\C]
    Bc = [zeros(n, n); M\I]

    return ContinuousStateSpace(Ac, Bc)
end

"""
    ss_modal_model(ωn, ξn, ϕn)

Generates a continuous-time state-space model from the mass, damping, and stiffness matrices

**Inputs**
* `ωn`: Natural angular frequencies
* `ξn`: Damping ratios
* `ϕn`: Mass-normalized mode shapes

**Output**
* `css`: ContinuousStateSpace
"""
function ss_modal_model(ωn, ξn, ϕn)

    m, n = size(ϕn)
    Ac = [zeros(n, n) I; -Diagonal(ωn.^2) -Diagonal(2ξn.*ωn)]
    Bc = [zeros(n, m); ϕn']

    return ContinuousStateSpace(Ac, Bc)
end

"""
    eigenmode(Ac, n)

Computes the eigen of a continuous-time state-space model

**Inputs**
* `Ac`: Continuous-time state-space matrix
* `n`: Number of mode pairs to keep in the modal basis (default: size(Ac, 1) ÷ 2)

*Note: `n` is the number of mode pairs to keep in the basis*

**Outputs**
* `λ`: Eigenvalues
* `Ψ`: Eigenvectors
"""

function eigenmode(Ac, n::Int = size(Ac, 1) ÷ 2)

    # For mechanical systems with real-valued state-space matrices, the eigenvalues and eigenvectors come in complex conjugate pairs. To ensure type stability, we determine the appropriate type for the eigenvalues and eigenvectors based on whether `Ac` is real or complex.
    t = if isreal(Ac)
        Complex{eltype(Ac)}
    else
        eltype(Ac)
    end

    # Compute eigenmodes and ensure type stability
    res = eigen(Ac)
    λ, Ψ = t.(res.values), t.(res.vectors)
    sort_idx = sortperm(abs.(λ))

    return λ[sort_idx][1:2n], Ψ[:, sort_idx][:, 1:2n]
end

"""
    modal_parameters(λ)

Computes the natural angular frequencies and damping ratios from the complex eigenvalues

**Input**
* `λ`: Complex eigenvalues

**Outputs**
* `ωn`: Natural angular frequencies
* `ξn`: Damping ratios
"""
function modal_parameters(λ)
    λn = λ[2:2:end]
    ωn = abs.(λn)
    ξn = -real(λn)./ωn

    return ωn, ξn
end
