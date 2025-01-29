"""
    StateSpaceProblem(css::ContinuousStateSpace, u0::Vector{Float64}, h::Float64, F::Matrix{Float64})

Structure containing data for the state-space model

# Constructor parameters
* `css`: Continuous-time state space model
* `u0`: Initial conditions
* `h`: Time step
* `F`: External force matrix

# Fields
* css: ContinuousStateSpace
* u0: Initial conditions
* h: Time step
* F: External force matrix
"""
struct StateSpaceTimeProblem
    css :: ContinuousStateSpace
    u0 :: Vector{Float64}
    h :: Float64
    F :: Matrix{Float64}
end

"""
    StateSpaceTimeSolution(D, V, A)

Structure containing the solution of the state-space model

# Fields
* `D`: Displacement matrix or vector
* `V`: Velocity matrix or vector
* `A`: Acceleration matrix or vector
"""
struct StateSpaceTimeSolution
    D
    V
    A
end

"""
    StateSpaceFRFProblem(css::ContinuousStateSpace, freq, type::Symbol)

Structure containing the data feeding the direct solver for calculating an FRF

# Fields
* `css`: Continuous-time state space model
* `freq`: Frequencies of interest
* `Sₒ`: Selection matrix for observation points
* `Sₑ`: Selection matrix for excitation points

# Note: It is assumed that the output equation is of the form y = Sₒ*x
"""
struct StateSpaceFRFProblem
    css :: ContinuousStateSpace
    freq
    Sₒ
    Sₑ

    StateSpaceFRFProblem(css, freq, Sₒ = I(Int(size(css.Ac, 1)/2)), Sₑ = I(size(css.Bc, 2))) = new(css, freq, Sₒ, Sₑ)
end

"""
    StateSpaceFRFSolution(y)

Structure containing the solution of the state-space model

# Fields
* `y`: FRF matrix
"""
struct StateSpaceFRFSolution
    y :: Union{Matrix{ComplexF64}, Vector{Matrix{ComplexF64}}}
end

"""
    StateSpaceFreqProblem(css::ContinuousStateSpace, freq, F, Sₒ)

Structure containing the data feeding the modal solver for calculating the frequency reponse

# Fields
* `css`: Continuous-time state space model
* `freq`: Frequencies of interest
* `F`: External force matrix
* `Sₒ`: Selection matrix for observation points

# Output
* `StateSpaceFreqSolution`: Solution of the state-space model
"""
struct StateSpaceFreqProblem
    css :: ContinuousStateSpace
    freq
    F
    Sₒ

    StateSpaceFreqProblem(css, freq, F, Sₒ = I(Int(size(css.Ac, 1)/2))) = new(css, freq, F, Sₒ)
end

struct StateSpaceFreqSolution
    y
end

"""
    solve(prob::StateSpaceProblem, method = :zoh)

Solves a discrete-time problem using the state-space model

# Inputs
* `prob`: Discrete-time problem
* `method::Symbol`: Discretization method
    * `:zoh`: Zero-order Hold method
    * `:foh`: First-order Hold method
    * `:blh`: Band-limited Hold method
    * `:rk4`: Runge-Kutta 4th order method

# Output
* `StateSpaceSolution`: Solution of the state-space model
"""
function solve(prob::StateSpaceTimeProblem, method::Symbol = :zoh)
    (; css, u0, h, F) = prob
    dss = c2d(css, h, method)

    n, nt = size(F)
    m = Int(n/2)
    x = zeros(n, nt)
    A = zeros(m, nt)

    x[:, 1] .= u0[:]
    A[:, 1] .= css.Bd[m+1:end]*F[:, 1] .+ css.Ac[m+1:end, :]*x[:, 1]

    @views @inbounds for k in 1:nt-1
        if method == :zoh || method == :blh
            x[:, k+1] .= dss.Ad*x[:, k] .+ dss.Bd*F[:, k]
        elseif method == :foh || method == :rk4
            x[:, k+1] .= dss.Ad*x[:, k] .+ dss.Bd*F[:, k] .+ dss.Bdp*F[:, k+1]
        end

        A[:, k+1] .= css.Bc[m+1:end]*F[:, k+1] .+ css.Ac[m+1:end, :]*x[:, k+1]
    end

    return StateSpaceTimeSolution(x[1:m, :], x[m+1:end, :], A)
end


"""
    solve(prob::StateTimeSpaceProblem, alg = RK4())

Solves a continuous-time problem using the state-space model

# Inputs
* `prob`: Continuous-time problem
* `alg`: Time integration algorithm

# Output
* `StateSpaceSolution`: Solution of the state-space model
"""
function solve(prob::StateSpaceTimeProblem, alg::RK4)
    (; css, u0, h, F) = prob
    (; Ac, Bc) = css

    n, nt = size(F)
    m = Int(n/2)
    x = undefs(n, nt)
    A = undefs(m, nt)

    # Intermediate vectors
    k₁ = undefs(n)
    k₂ = undefs(n)
    k₃ = undefs(n)
    k₄ = undefs(n)
    Fn_2 = undefs(n)

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

    return StateSpaceTimeSolution(x[1:m, :], x[m+1:end, :], A)
end

"""
    solve(m::StateSpaceFRFProblem, type = :dis, ismat = false)

Computes the FRF matrix by direct method

# Parameter
* m: Structure containing the problem data
* type: Type of FRF to compute (:dis, :vel, :acc)
* ismat: Return the FRF matrix as a 3D array (default = false)

# Output
* sol: FRFSolution structure
"""
function solve(m::StateSpaceFRFProblem, type = :dis, ismat = false)

    # Initialisation
    (; css, freq, Sₒ, Sₑ) = m
    Ac, Bc = css
    Nₒ = size(Sₒ, 1)
    Nₑ = size(Sₑ, 1)
    Nstate, Nu = size(Bc)
    Ns = Int(Nstate/2)
    Nf = length(freq)

    FRF = [undefs(ComplexF64, Nₒ, Nₑ) for _ in 1:Nf]
    M = undefs(ComplexF64, Nstate, Nu)

    if type == :acc
        C = S₀*Ac[Ns+1:end, :]
        D = S₀*Bc[Ns+1:end, :]
    end

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "FRF calculation - Direct method...", showspeed = true)
    for (f, ω) in enumerate(ωf)
        next!(p)
        M .= (1im*ω*I - Ac)\Bc

        if type == :dis
            FRF[f] .= Sₒ*M[1:Ns, :]*Sₑ'
        elseif type == :vel
            FRF[f] .= Sₒ*M[Ns+1:end, :]*Sₑ'
        elseif type == :acc
            FRF[f] .= (C*M .+ D)*Sₑ'
        end
    end

    if ismat
        return reshape(reduce(hcat, FRF), Nₒ, Nₑ, :)
    end

    return StateSpaceFRFSolution(FRF)
end

"""
    solve(m::StateSpaceFRFProblem, Nₘ::Int, type = :dis, ismat = false)

Computes the FRF matrix by modal method

# Parameter
* m: Structure containing the problem data
* Nₘ: Number of eigenmodes to keep in the modal basis
* type: Type of FRF to compute (:dis, :vel, :acc)
* ismat: Return the FRF matrix as a 3D array (default = false)

# Output
* sol: StateSpaceFRFSolution
"""
function solve(m::StateSpaceFRFProblem, Nₘ::Int, type = :dis, ismat = false)

    # Initialisation
    (; css, freq, Sₒ, Sₑ) = m
    Nₒ = size(Sₒ, 1)
    Nₑ = size(Sₑ, 1)
    Nstate = size(css.Ac, 1)
    Ns = Int(Nstate/2)
    Nf = length(freq)

    # Compute modal information
    λ, Ψ = eigenmode(Ac, Nₘ)
    B = Ψ\css.Bc

    if type == :dis
        Ψₒ = Sₒ*Ψ[1:Ns, :]
    elseif type == :vel
        Ψₒ = Sₒ*Ψ[Ns+1:end, :]
    elseif type == :acc
        C = Sₒ*Ac[Ns+1:end, :]*Ψ
        D = Sₒ*Bc[Ns+1:end, :]*Sₑ'
    end
    Bₑ = B*Sₑ'
    Nm = length(λ)

    FRF = [undefs(ComplexF64, Nₒ, Nₑ) for _ in 1:Nf]
    M = undefs(Diagonal{ComplexF64}, Nm)
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "FRF calculation - Modal approach...", showspeed = true)
    for (f, ω) in enumerate(ωf)
        next!(p)
        @. M[indm] = 1/(1im*ω - λ)

        if type == :dis || type == :vel
            FRF[f] .= Ψₒ*M*Bₑ
        elseif type == :acc
            FRF[f] .*= C*M*Bₑ + D
        end
    end

    if ismat
        return reshape(reduce(hcat, FRF), Nₒ, Nₑ, :)
    end

    return StateSpaceFRFSolution(FRF)
end

"""
    solve(m::StateSpaceFreqProblem, type = :dis)

Computes the frequency response by direct method

# Inputs
* `m`: Structure containing the problem data
* `type::Symbol`: Type of FRF to compute
    * `:dis`: Displacement
    * `:vel`: Velocity
    * `:acc`: Acceleration

# Output
`sol`: StateSpaceFreqSolution
"""
function solve(m::StateSpaceFreqProblem, type = :dis)
    # Initialisation
    (; css, freq, F, Sₒ) = m
    Ac, Bc = css
    Nₒ = size(Sₒ, 1)
    Nstate, Nu = size(Bc)
    Ns = Int(Nstate/2)
    Nf = length(freq)

    y = undefs(ComplexF64, Nₒ, Nf)
    M = undefs(ComplexF64, Nstate, Nu)

    if type == :acc
        C = Sₒ*Ac[Ns+1:end, :]
        D = Sₒ*Bc[Ns+1:end, :]
    end

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "Frequency response - Direct method...", showspeed = true)
    for (f, (ω, Fₑ)) in enumerate(zip(ωf, eachcol(F)))
        next!(p)
        M .= (1im*ω*I - Ac)\Bc

        if type == :dis
            y[:, f] .= Sₒ*M[1:Ns, :]*Fₑ
        elseif type == :vel
            y[:, f] .= Sₒ*M[Ns+1:end, :]*Fₑ
        elseif type == :acc
            y[:, f] .= (C*M + D)*Fₑ
        end
    end

    return StateSpaceFreqSolution(y)
end

"""
    solve(m::StateSpaceFreqProblem, Nₘ = Int[], type = :dis)

Computes the frequency response by modal method

# Inputs
* `m`: Structure containing the problem data
* `Nₘ`: Number of eigenmodes to keep in the modal basis
* `type::Symbol`: Type of FRF to compute
    * `:dis`: Displacement
    * `:vel`: Velocity
    * `:acc`: Acceleration
* `Nₘ`: Number of eigenmodes to keep in the modal basis

# Output
`sol`: StateSpaceFreqSolution
"""
function solve(m::StateSpaceFreqProblem, Nₘ::Int, type = :dis)
    # Initialisation
    (; css, freq, F, Sₒ) = m
    Nₒ = size(Sₒ, 1)
    Nstate = size(css.Ac, 1)
    Ns = Int(Nstate/2)
    Nf = length(freq)

    # Compute modal information
    λ, Ψ = eigenmode(Ac, Nₘ)
    Bc = Ψ\css.Bc

    if type == :dis || type == :acc
        Ψₒ = Sₒ*Ψ[1:Ns, :]
    elseif type == :vel
        Ψₒ = Sₒ*Ψ[Ns+1:end, :]
    end
    u = Bc*F
    Nm = length(λ)

    y = undefs(ComplexF64, Nₒ, Nf)
    M = Diagonal{ComplexF64}(undef, Nm)
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "Frequency response - Modal approach...", showspeed = true)
    for (f, (ω, uₑ)) in enumerate(zip(ωf, eachcol(u)))
        next!(p)
        @. M[indm] = 1/(1im*ω - λ)

        y[:, f] .= Ψₒ*M*uₑ
        if type == :acc
            y[:, f] .*= -ω^2
        end
    end

    return StateSpaceFreqSolution(y)
end