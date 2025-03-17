"""
    StateSpaceProblem(css::ContinuousStateSpace, u0::Vector{Float64}, h::Float64, F::Matrix{Float64})

Structure containing data for the state-space model

# Constructor
* `css`: Continuous-time state space model
* `u0`: Initial conditions
* `t`: Time vector
* `F`: External force matrix

# Fields
* css: ContinuousStateSpace
* u0: Initial conditions
* h: Time step
* F: External force matrix
"""
@show_struct struct StateSpaceTimeProblem{T <: Real}
    css::ContinuousStateSpace
    u0::Vector{T}
    h::T
    F::Matrix{T}

    StateSpaceTimeProblem(css, u0::Vector{T}, t::AbstractVector, F::Matrix{T}) where T = new{T}(css, u0, t[2] - t[1], F)
end

"""
    StateSpaceTimeSolution(u, du, ddu)

Structure containing the solution of the state-space model

# Fields
* `u`: Displacement matrix or vector
* `du`: Velocity matrix or vector
* `ddu`: Acceleration matrix or vector
"""
@show_struct struct StateSpaceTimeSolution{T <: Real}
    u::Matrix{T}
    du::Matrix{T}
    ddu::Matrix{T}
end

"""
    StateSpaceFRFProblem(css::ContinuousStateSpace, freq, type, So, Se)

Structure containing the data feeding the direct solver for calculating an FRF

# Fields
* `css`: Continuous-time state space model
* `freq`: Frequencies of interest
* `So`: Selection matrix for observation points
* `Se`: Selection matrix for excitation points

# note: It is assumed that the output equation is of the form y = So*x
"""
@show_struct struct StateSpaceFRFProblem{T <: AbstractVector, S <: AbstractMatrix}
    css::ContinuousStateSpace
    freq::T
    So::S
    Se::S

    StateSpaceFRFProblem(css, freq::T, So::S = I(Int(size(css.Ac, 1)/2)), Se::S = I(size(css.Bc, 2))) where {T, S} = new{T, S}(css, freq, So, Se)
end

"""
    StateSpacemodalFRFProblem(css::ContinuousStateSpace, freq, type, So, Se, n)

Structure containing the data feeding the direct solver for calculating an FRF

# Fields
* `css`: Continuous-time state space model
* `freq`: Frequencies of interest
* `So`: Selection matrix for observation points
* `Se`: Selection matrix for excitation points
* `n`: Number of modes to keep in the modal basis

# note: It is assumed that the output equation is of the form y = So*x
"""
@show_struct struct StateSpaceModalFRFProblem{T <: AbstractVector, S <: AbstractMatrix}
    css::ContinuousStateSpace
    freq::T
    So::S
    Se::S
    n::Int

    StateSpaceModalFRFProblem(css, freq::T, So::S = I(Int(size(css.Ac, 1)/2)), Se::S = I(size(css.Bc, 2)), n = 0) where {T, S} = new{T, S}(css, freq, So, Se, n)
end

"""
    StateSpaceFRFSolution(u)

Structure containing the solution of the state-space model

# Field
* `u`: FRF matrix
"""
@show_struct struct StateSpaceFRFSolution{T <: Complex}
    u::Union{Matrix{T}, Vector{Matrix{T}}}
end

"""
    StateSpaceFreqProblem(css::ContinuousStateSpace, freq, F, So)

Structure containing the data feeding the modal solver for calculating the frequency reponse

# Fields
* `css`: Continuous-time state space model
* `freq`: Frequencies of interest
* `F`: External force matrix
* `So`: Selection matrix for observation points
"""
@show_struct struct StateSpaceFreqProblem{T <: AbstractVector, S <: AbstractMatrix, U <: AbstractMatrix}
    css::ContinuousStateSpace
    freq::T
    F::S
    So::U

    StateSpaceFreqProblem(css, freq::T, F::S, So::U = I(Int(size(css.Ac, 1)/2))) where {T, S, U} = new{T, S, U}(css, freq, F, So)
end

"""
    StateSpaceModalFreqProblem(css::ContinuousStateSpace, freq, F, So, n)

Structure containing the data feeding the modal solver for calculating the frequency reponse

# Fields
* `css`: Continuous-time state space model
* `freq`: Frequencies of interest
* `F`: External force matrix
* `So`: Selection matrix for observation points
* `n`: Number of modes to keep in the modal basis
"""
@show_struct struct StateSpaceModalFreqProblem{T <: AbstractVector, S <: AbstractMatrix, U <: AbstractMatrix}
    css::ContinuousStateSpace
    freq::T
    F::S
    So::U
    n::Int

    StateSpaceModalFreqProblem(css, freq::T, F::S, So::U = I(Int(size(css.Ac, 1)/2)), n = 0) where {T, S, U} = new{T, S, U}(css, freq, F, So, n)
end


"""
    StateSpaceFreqSolution(u)

Structure containing the solution of the state-space model

# Field
* `u`: Frequency response matrix
"""
@show_struct struct StateSpaceFreqSolution{T <: Complex}
    u::Matrix{T}
end

"""
    solve(prob::StateSpaceTimeProblem, method = :zoh; progress = true)

Solves a discrete-time problem using the state-space model

# Inputs
* `prob`: Discrete-time problem
* `method::Symbol`: Discretization method
    * `:zoh`: Zero-order Hold method
    * `:foh`: First-order Hold method
    * `:blh`: Band-limited Hold method
    * `:rk4`: Runge-Kutta 4th order method
* `progress`: Show progress bar (default = true)

# Output
* `StateSpaceSolution`: Solution of the state-space model
"""
function solve(prob::StateSpaceTimeProblem, method::Symbol = :zoh; progress = true)

    (; css, u0, h, F) = prob
    dss = c2d(css, h, method)

    nx = size(css.Ac, 1)
    nt = size(F, 2)
    m = Int(nx/2)
    x = similar(F, nx, nt)
    A = similar(F, m, nt)

    x[:, 1] .= u0[:]
    A[:, 1] .= css.Bc[m+1:end, :]*F[:, 1] .+ css.Ac[m+1:end, :]*x[:, 1]

    if method == :zoh
        name = "State Space Time problem - ZOH..."
    elseif method == :foh
        name = "Time problem - FOH..."
    elseif method == :blh
        name = "Time problem - BLH..."
    elseif method == :rk4
        name = "Time problem - RK4..."
    end

    p = Progress(nt, color = :black, desc = name, showspeed = true)
    @views @inbounds for k in 1:nt-1
        progress ? next!(p) : nothing

        if method == :zoh || method == :blh
            x[:, k+1] .= dss.Ad*x[:, k] .+ dss.Bd*F[:, k]
        elseif method == :foh || method == :rk4
            x[:, k+1] .= dss.Ad*x[:, k] .+ dss.Bd*F[:, k] .+ dss.Bdp*F[:, k+1]
        end

        A[:, k+1] .= css.Bc[m+1:end, :]*F[:, k+1] .+ css.Ac[m+1:end, :]*x[:, k+1]
    end

    return StateSpaceTimeSolution(x[1:m, :], x[m+1:end, :], A)
end


"""
    solve(prob::StateTimeSpaceProblem, alg = RK4(); progress = true)

Solves a continuous-time problem using the state-space model

# Inputs
* `prob`: Continuous-time problem
* `alg`: Time integration algorithm
* `progress`: Show progress bar (default = true)

# Output
* `StateSpaceSolution`: Solution of the state-space model
"""
function solve(prob::StateSpaceTimeProblem, alg::RK4; progress = true)
    (; css, u0, h, F) = prob
    (; Ac, Bc) = css

    nx = size(css.Ac, 1)
    nt = size(F, 2)
    m = Int(nx/2)
    x = similar(F, nx, nt)
    A = similar(F, m, nt)

    # Intermediate vectors
    k₁ = similar(F, nx)
    k₂ = similar(k₁)
    k₃ = similar(k₁)
    k₄ = similar(k₁)
    Fn_2 = similar(F, m)

    x[:, 1] .= u0[:]

    p = Progress(nt, color = :black, desc = "State Space Time Problem - RK4...", showspeed = true)
    @views @inbounds for k in 1:nt-1
        progress ? next!(p) : nothing

        Fn_2 .= (F[:, k] .+ F[:, k+1])/2

        k₁ .= Ac*x[:, k] .+ Bc*F[:, k]
        k₂ .= Ac*(x[:, k] .+ h*k₁/2) .+ Bc*Fn_2
        k₃ .= Ac*(x[:, k] .+ h*k₂/2) .+ Bc*Fn_2
        k₄ .= Ac*(x[:, k] .+ h*k₃) .+ Bc*F[:, k+1]

        @. x[:, k+1] = x[:, k] + h*(k₁ + 2k₂ + 2k₃ + k₄)/6
        A[:, k+1] .= css.Bc[m+1:end, :]*F[:, k+1] .+ css.Ac[m+1:end, :]*x[:, k+1]
    end

    return StateSpaceTimeSolution(x[1:m, :], x[m+1:end, :], A)
end

"""
    solve(m::StateSpaceFRFProblem, type = :dis; ismat = false, progress = true)

Computes the FRF matrix by direct method

# Parameter
* `m`: Structure containing the problem data
* `type`: Type of FRF to compute (:dis, :vel, :acc)
* `ismat`: Return the FRF matrix as a 3D array (default = false)
* `progress`: Show progress bar (default = true)

# Output
* `sol`: FRFSolution structure
"""
function solve(m::StateSpaceFRFProblem, type::Symbol = :dis; ismat = false, progress = true)

    # Initialisation
    (; css, freq, So, Se) = m
    (; Ac, Bc) = css
    no = size(So, 1)
    Ne = size(Se, 1)
    nstate, Nu = size(Bc)
    ns = Int(nstate/2)
    nf = length(freq)

    FRF = [similar([], Complex{eltype(Ac)}, no, Ne) for _ in 1:nf]
    M = similar([], Complex{eltype(Ac)}, nstate, Nu)

    if type == :acc
        C = S₀*Ac[ns+1:end, :]
        D = S₀*Bc[ns+1:end, :]
    end

    ωf = 2π*freq
    p = Progress(nf, color = :black, desc = "State Space FRF Problem - Direct method...", showspeed = true)
    for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        M .= (1im*ω*I - Ac)\Bc

        if type == :dis
            FRF[f] .= So*M[1:ns, :]*Se'
        elseif type == :vel
            FRF[f] .= So*M[ns+1:end, :]*Se'
        elseif type == :acc
            FRF[f] .= (C*M .+ D)*Se'
        end
    end

    if ismat
        return reshape(reduce(hcat, FRF), no, Ne, :)
    end

    return StateSpaceFRFSolution(FRF)
end

"""
    solve(m::StateSpaceModalFRFProblem, type = :dis; ismat = false, progress = true)

Computes the FRF matrix by modal method

# Parameter
* `m`: Structure containing the problem data
* `type`: Type of FRF to compute (:dis, :vel, :acc)
* `ismat`: Return the FRF matrix as a 3D array (default = false)
* `progress`: Show progress bar

# Output
* `sol`: StateSpaceFRFSolution
"""
function solve(m::StateSpaceModalFRFProblem, type = :dis; ismat = false, progress = true)

    # Initialisation
    (; css, freq, So, Se, n) = m
    (; Ac, Bc) = css
    no = size(So, 1)
    Ne = size(Se, 1)
    nstate = size(css.Ac, 1)
    ns = Int(nstate/2)
    nf = length(freq)

    # Compute modal information
    λ, Ψ = eigenmode(Ac, n)
    B = Ψ\Bc

    if type == :dis
        Ψo = So*Ψ[1:ns, :]
    elseif type == :vel
        Ψo = So*Ψ[ns+1:end, :]
    elseif type == :acc
        C = So*Ac[ns+1:end, :]*Ψ
        D = So*Bc[ns+1:end, :]*Se'
    end
    Be = B*Se'
    n = length(λ)

    FRF = [similar([], Complex{eltype(Ac)}, no, Ne) for _ in 1:nf]
    M = Diagonal(similar([], Complex{eltype(Ac)}, n))
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(nf, color = :black, desc = "State Space FRF Problem - Modal approach...", showspeed = true)
    for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        @. M[indm] = 1/(1im*ω - λ)

        if type == :dis || type == :vel
            FRF[f] .= Ψo*M*Be
        elseif type == :acc
            FRF[f] .*= C*M*Be + D
        end
    end

    if ismat
        return reshape(reduce(hcat, FRF), no, Ne, :)
    end

    return StateSpaceFRFSolution(FRF)
end

"""
    solve(m::StateSpaceFreqProblem, type = :dis; progress = true)

Computes the frequency response by direct method

# Inputs
* `m`: Structure containing the problem data
* `type::Symbol`: Type of FRF to compute
    * `:dis`: Displacement
    * `:vel`: Velocity
    * `:acc`: Acceleration
* `progress`: Show progress bar

# Output
`sol`: StateSpaceFreqSolution
"""
function solve(m::StateSpaceFreqProblem, type = :dis; progress = true)
    # Initialisation
    (; css, freq, F, So) = m
    (; Ac, Bc) = css
    no = size(So, 1)
    nstate, Nu = size(Bc)
    ns = Int(nstate/2)
    nf = length(freq)

    y = similar([], Complex{eltype(Ac)}, no, nf)
    M = similar([], Complex{eltype(Ac)}, nstate, Nu)

    if type == :acc
        C = So*Ac[ns+1:end, :]
        D = So*Bc[ns+1:end, :]
    end

    ωf = 2π*freq
    p = Progress(nf, color = :black, desc = "State Space Frequency Problem - Direct method...", showspeed = true)
    for (f, (ω, Fe)) in enumerate(zip(ωf, eachcol(F)))
        progress ? next!(p) : nothing

        M .= (1im*ω*I - Ac)\Bc

        if type == :dis
            y[:, f] .= So*M[1:ns, :]*Fe
        elseif type == :vel
            y[:, f] .= So*M[ns+1:end, :]*Fe
        elseif type == :acc
            y[:, f] .= (C*M + D)*Fe
        end
    end

    return StateSpaceFreqSolution(y)
end

"""
    solve(m::StateSpaceFreqProblem, type::Symbol = :dis; progress = true)

Computes the frequency response by modal method

# Inputs
* `m`: Structure containing the problem data
* `n`: Number of eigenodes to keep in the modal basis
* `type::Symbol`: Type of FRF to compute
    * `:dis`: Displacement
    * `:vel`: Velocity
    * `:acc`: Acceleration
* `progress`: Show progress bar

# Output
`sol`: StateSpaceFreqSolution
"""
function solve(m::StateSpaceModalFreqProblem, type = :dis; progress = true)
    # Initialisation
    (; css, freq, F, So, n) = m
    (; Ac, Bc) = css
    no = size(So, 1)
    nstate = size(Ac, 1)
    ns = Int(nstate/2)
    nf = length(freq)

    # Compute modal information
    λ, Ψ = eigenmode(Ac, n)
    Bc = Ψ\css.Bc

    if type == :dis || type == :acc
        Ψo = So*Ψ[1:ns, :]
    elseif type == :vel
        Ψo = So*Ψ[ns+1:end, :]
    end
    u = Bc*F
    n = length(λ)

    y = similar([], Complex{eltype(Ac)}, no, nf)
    M = Diagonal(similar([], Complex{eltype(Ac)}, n))
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(nf, color = :black, desc = "State Space Frequency Problem - Modal approach...", showspeed = true)
    for (f, (ω, ue)) in enumerate(zip(ωf, eachcol(u)))
        progress ? next!(p) : nothing

        @. M[indm] = 1/(1im*ω - λ)

        y[:, f] .= Ψo*M*ue
        if type == :acc
            y[:, f] .*= -ω^2
        end
    end

    return StateSpaceFreqSolution(y)
end