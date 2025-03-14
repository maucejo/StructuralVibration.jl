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
@with_kw struct StateSpaceTimeProblem
    css :: ContinuousStateSpace
    u0 :: Vector{Float64}
    h :: Float64
    F :: Matrix{Float64}

    StateSpaceTimeProblem(css, u0, t, F) = new(css, u0, t[2] - t[1], F)
end

"""
    StateSpaceTimeSolution(u, du, ddu)

Structure containing the solution of the state-space model

# Fields
* `u`: Displacement matrix or vector
* `du`: Velocity matrix or vector
* `ddu`: Acceleration matrix or vector
"""
@with_kw struct StateSpaceTimeSolution
    u
    du
    ddu
end

"""
    StateSpaceFRFProblem(css::ContinuousStateSpace, freq, type, So, Se)

Structure containing the data feeding the direct solver for calculating an FRF

# Fields
* `css`: Continuous-time state space model
* `freq`: Frequencies of interest
* `So`: Selection matrix for observation points
* `Se`: Selection matrix for excitation points

# Note: It is assumed that the output equation is of the form y = So*x
"""
@with_kw struct StateSpaceFRFProblem
    css :: ContinuousStateSpace
    freq
    So
    Se

    StateSpaceFRFProblem(css, freq, So = I(Int(size(css.Ac, 1)/2)), Se = I(size(css.Bc, 2))) = new(css, freq, So, Se)
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

# Note: It is assumed that the output equation is of the form y = So*x
"""
@with_kw struct StateSpaceModalFRFProblem
    css :: ContinuousStateSpace
    freq
    So
    Se
    n

    StateSpaceModalFRFProblem(css, freq, So = I(Int(size(css.Ac, 1)/2)), Se = I(size(css.Bc, 2)), n = 0) = new(css, freq, So, Se, n)
end

"""
    StateSpaceFRFSolution(u)

Structure containing the solution of the state-space model

# Field
* `u`: FRF matrix
"""
@with_kw struct StateSpaceFRFSolution
    u :: Union{Matrix{ComplexF64}, Vector{Matrix{ComplexF64}}}
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
@with_kw struct StateSpaceFreqProblem
    css :: ContinuousStateSpace
    freq
    F
    So

    StateSpaceFreqProblem(css, freq, F, So = I(Int(size(css.Ac, 1)/2))) = new(css, freq, F, So)
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
@with_kw struct StateSpaceModalFreqProblem
    css :: ContinuousStateSpace
    freq
    F
    So
    n

    StateSpaceModalFreqProblem(css, freq, F, So = I(Int(size(css.Ac, 1)/2)), n = 0) = new(css, freq, F, So, n)
end


"""
    StateSpaceFreqSolution(u)

Structure containing the solution of the state-space model

# Field
* `u`: Frequency response matrix
"""
@with_kw struct StateSpaceFreqSolution
    u
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
    x = undefs(nx, nt)
    A = undefs(m, nt)

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
    x = undefs(nx, nt)
    A = undefs(m, nt)

    # Intermediate vectors
    k₁ = undefs(nx)
    k₂ = undefs(nx)
    k₃ = undefs(nx)
    k₄ = undefs(nx)
    Fn_2 = undefs(m)

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
    Nₒ = size(So, 1)
    Ne = size(Se, 1)
    Nstate, Nu = size(Bc)
    Ns = Int(Nstate/2)
    Nf = length(freq)

    FRF = [undefs(ComplexF64, Nₒ, Ne) for _ in 1:Nf]
    M = undefs(ComplexF64, Nstate, Nu)

    if type == :acc
        C = S₀*Ac[Ns+1:end, :]
        D = S₀*Bc[Ns+1:end, :]
    end

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "State Space FRF Problem - Direct method...", showspeed = true)
    for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        M .= (1im*ω*I - Ac)\Bc

        if type == :dis
            FRF[f] .= So*M[1:Ns, :]*Se'
        elseif type == :vel
            FRF[f] .= So*M[Ns+1:end, :]*Se'
        elseif type == :acc
            FRF[f] .= (C*M .+ D)*Se'
        end
    end

    if ismat
        return reshape(reduce(hcat, FRF), Nₒ, Ne, :)
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
    Nₒ = size(So, 1)
    Ne = size(Se, 1)
    Nstate = size(css.Ac, 1)
    Ns = Int(Nstate/2)
    Nf = length(freq)

    # Compute modal information
    λ, Ψ = eigenode(css.Ac, n)
    B = Ψ\css.Bc

    if type == :dis
        Ψₒ = So*Ψ[1:Ns, :]
    elseif type == :vel
        Ψₒ = So*Ψ[Ns+1:end, :]
    elseif type == :acc
        C = So*Ac[Ns+1:end, :]*Ψ
        D = So*Bc[Ns+1:end, :]*Se'
    end
    Be = B*Se'
    n = length(λ)

    FRF = [undefs(ComplexF64, Nₒ, Ne) for _ in 1:Nf]
    M = Diagonal(undefs(ComplexF64, n))
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "State Space FRF Problem - Modal approach...", showspeed = true)
    for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        @. M[indm] = 1/(1im*ω - λ)

        if type == :dis || type == :vel
            FRF[f] .= Ψₒ*M*Be
        elseif type == :acc
            FRF[f] .*= C*M*Be + D
        end
    end

    if ismat
        return reshape(reduce(hcat, FRF), Nₒ, Ne, :)
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
    Nₒ = size(So, 1)
    Nstate, Nu = size(Bc)
    Ns = Int(Nstate/2)
    Nf = length(freq)

    y = undefs(ComplexF64, Nₒ, Nf)
    M = undefs(ComplexF64, Nstate, Nu)

    if type == :acc
        C = So*Ac[Ns+1:end, :]
        D = So*Bc[Ns+1:end, :]
    end

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "State Space Frequency Problem - Direct method...", showspeed = true)
    for (f, (ω, Fe)) in enumerate(zip(ωf, eachcol(F)))
        progress ? next!(p) : nothing

        M .= (1im*ω*I - Ac)\Bc

        if type == :dis
            y[:, f] .= So*M[1:Ns, :]*Fe
        elseif type == :vel
            y[:, f] .= So*M[Ns+1:end, :]*Fe
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
    Nₒ = size(So, 1)
    Nstate = size(css.Ac, 1)
    Ns = Int(Nstate/2)
    Nf = length(freq)

    # Compute modal information
    λ, Ψ = eigenode(css.Ac, n)
    Bc = Ψ\css.Bc

    if type == :dis || type == :acc
        Ψₒ = So*Ψ[1:Ns, :]
    elseif type == :vel
        Ψₒ = So*Ψ[Ns+1:end, :]
    end
    u = Bc*F
    n = length(λ)

    y = undefs(ComplexF64, Nₒ, Nf)
    M = Diagonal(undefs(ComplexF64, n))
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "State Space Frequency Problem - Modal approach...", showspeed = true)
    for (f, (ω, ue)) in enumerate(zip(ωf, eachcol(u)))
        progress ? next!(p) : nothing

        @. M[indm] = 1/(1im*ω - λ)

        y[:, f] .= Ψₒ*M*ue
        if type == :acc
            y[:, f] .*= -ω^2
        end
    end

    return StateSpaceFreqSolution(y)
end