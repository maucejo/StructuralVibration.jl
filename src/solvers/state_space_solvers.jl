"""
    StateSpaceProblem(css::ContinuousStateSpace, u0::Vector{Float64}, h::Float64, F::Matrix{Float64})

Structure containing data for the state-space model

**Constructor**
* `css`: Continuous-time state space model
* `u0`: Initial conditions
* `F`: External force matrix
* `t`: Time vector

**Fields**
* `css`: ContinuousStateSpace
* `F`: External force matrix
* `u0`: Initial conditions
* `h`: Time step
"""
@show_data struct StateSpaceTimeProblem{Tf <: AbstractMatrix, Tu <: AbstractVector, T <: Real}
    css::ContinuousStateSpace
    F::Tf
    u0::Tu
    h::T

    function StateSpaceTimeProblem(css, F::Tf, u0::Tu, t::AbstractRange) where {Tf, Tu}

        h = step(t)
        return new{Tf, Tu, typeof(h)}(css, F, u0, h)
    end
end

"""
    StateSpaceTimeSolution(u, du, ddu)

Structure containing the solution of the state-space model

**Fields**
* `u`: Displacement matrix or vector
* `du`: Velocity matrix or vector
* `ddu`: Acceleration matrix or vector
"""
@show_data struct StateSpaceTimeSolution{T <: Real}
    u::Matrix{T}
    du::Matrix{T}
    ddu::Matrix{T}
end

"""
    StateSpaceFRFProblem(css::ContinuousStateSpace, freq, type, So, Se)

Structure containing the data feeding the direct solver for calculating an FRF

**Fields**
* `css`: Continuous-time state space model
* `freq`: Frequencies of interest
* `So`: Selection matrix for observation points
* `Se`: Selection matrix for excitation points

**Note**
It is assumed that the output equation is of the form y = So*x
"""
@show_data struct StateSpaceFRFProblem{Tf <: AbstractRange, Ts <: AbstractMatrix}
    css::ContinuousStateSpace
    freq::Tf
    So::Ts
    Se::Ts

    StateSpaceFRFProblem(css, freq::Tf, So::Ts = I(Int(size(css.Ac, 1)/2)), Se::Ts = I(size(css.Bc, 2))) where {Tf, Ts} = new{Tf, Ts}(css, freq, So, Se)
end

"""
    StateSpacemodalFRFProblem(css::ContinuousStateSpace, freq, type, So, Se, n)

Structure containing the data feeding the direct solver for calculating an FRF

**Fields**
* `css`: Continuous-time state space model
* `freq`: Frequencies of interest
* `So`: Selection matrix for observation points
* `Se`: Selection matrix for excitation points
* `n`: Number of modes to keep in the modal basis
"""
@show_data struct StateSpaceModalFRFProblem{Tf <: AbstractRange, Ts <: AbstractMatrix}
    css::ContinuousStateSpace
    freq::Tf
    So::Ts
    Se::Ts
    n::Int

    function StateSpaceModalFRFProblem(css, freq::Tf, So::Ts = I(Int(size(css.Ac, 1)/2)), Se::Ts = I(size(css.Bc, 2)), n = size(css.Ac, 2)) where {Tf, Ts}

        isodd(n) || n == 0 ? throw(ArgumentError("n must be even and greater than 0")) : nothing

        return new{Tf, Ts}(css, freq, So, Se, n)
    end
end

"""
    StateSpaceFRFSolution(u)

Structure containing the solution of the state-space model

**Field**
* `u`: FRF matrix
"""
@show_data struct StateSpaceFRFSolution{T <: Complex}
    u::Union{Array{T, 3}, Vector{Matrix{T}}}
end

"""
    StateSpaceFreqProblem(css::ContinuousStateSpace, F, freq, So)

Structure containing the data feeding the modal solver for calculating the frequency response

**Fields**
* `css`: Continuous-time state space model
* `F`: External force matrix
* `freq`: Frequencies of interest
* `So`: Selection matrix for observation points
"""
@show_data struct StateSpaceFreqProblem{TF <: AbstractMatrix, Tf <: AbstractRange, Ts <: AbstractMatrix}
    css::ContinuousStateSpace
    F::TF
    freq::Tf
    So::Ts

    StateSpaceFreqProblem(css, F::TF, freq::Tf, So::Ts = I(Int(size(css.Ac, 1)/2))) where {TF, Tf, Ts} = new{TF, Tf, Ts}(css, F, freq, So)
end

"""
    StateSpaceModalFreqProblem(css::ContinuousStateSpace, F, freq, So, n)

Structure containing the data feeding the modal solver for calculating the frequency response

**Fields**
* `css`: Continuous-time state space model
* `F`: External force matrix
* `freq`: Frequencies of interest
* `So`: Selection matrix for observation points
* `n`: Number of modes to keep in the modal basis
"""
@show_data struct StateSpaceModalFreqProblem{TF <: AbstractMatrix, Tf <: AbstractRange, Ts <: AbstractMatrix}
    css::ContinuousStateSpace
    F::TF
    freq::Tf
    So::Ts
    n::Int

    function StateSpaceModalFreqProblem(css, F::TF, freq::Tf, So::Ts = I(Int(size(css.Ac, 1)/2)), n = size(css.Ac, 2)) where {TF, Tf, Ts}

        isodd(n) || n == 0 ? throw(ArgumentError("n must be even and greater than 0")) : nothing

        return new{TF, Tf, Ts}(css, F, freq, So, n)
end

"""
    StateSpaceFreqSolution(u)

Structure containing the solution of the state-space model

**Field**
* `u`: Frequency response matrix
"""
@show_data struct StateSpaceFreqSolution{T <: Complex}
    u::Matrix{T}
end

"""
    solve(prob::StateSpaceTimeProblem, method = :zoh; progress = true)

Solves a discrete-time problem using the state-space model

**Inputs**
* `prob`: Discrete-time problem
* `method`: Discretization method
    * `:zoh`: Zero-order Hold method
    * `:foh`: First-order Hold method
    * `:blh`: Band-limited Hold method
    * `:rk4`: Runge-Kutta 4th order method
* `progress`: Show progress bar (default = true)

**Output**
* `StateSpaceSolution`: Solution of the state-space model
"""
function solve(prob::StateSpaceTimeProblem, method = :zoh; progress = true)

    (; css, F, u0, h) = prob
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

    p = Progress(nt, desc = name, showspeed = true)
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

**Inputs**
* `prob`: Continuous-time problem
* `alg`: Time integration algorithm
* `progress`: Show progress bar (default = true)

**Output**
* `StateSpaceSolution`: Solution of the state-space model
"""
function solve(prob::StateSpaceTimeProblem, alg::RK4; progress = true)
    (; css, u0, F, h) = prob
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

    p = Progress(nt, desc = "State Space Time Problem - RK4...", showspeed = true)
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
    solve(prob::StateSpaceFRFProblem; type = :dis, ismat = false, progress = true)
    solve(prob::StateSpaceModalFRFProblem; type = :dis, ismat = false, progress = true)

Computes the FRF matrix by direct or modal method

**Inputs**
* `prob`: Structure containing the problem data
* `type`: Type of FRF to compute
    * `:dis`: Admittance
    * `:vel`: Mobility
    * `:acc`: Accelerance
* `ismat`: Return the FRF matrix as a 3D array (default = false)
* `progress`: Show progress bar (default = true)

**Output**
* `sol`: FRFSolution structure
    * `u`: FRF matrix
"""
function solve(prob::StateSpaceFRFProblem; type = :dis, ismat = false, progress = true)

    # Initialisation
    (; css, freq, So, Se) = prob
    (; Ac, Bc) = css
    no = size(So, 1)
    ne = size(Se, 1)
    nstate, nu = size(Bc)
    ns = Int(nstate/2)
    nf = length(freq)

    FRF = Matrix{Complex}[similar(Ac, Complex{eltype(Ac)}, no, ne) for _ in 1:nf]
    M = similar(Ac, Complex{eltype(Ac)}, nstate, nu)

    if type == :acc
        C = So*Ac[ns+1:end, :]
        D = So*Bc[ns+1:end, :]
    end

    ωf = 2π*freq
    p = Progress(nf, desc = "State Space FRF Problem - Direct method...", showspeed = true)
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
        return StateSpaceFRFSolution(reshape(reduce(hcat, FRF), no, ne, :))
    end

    return StateSpaceFRFSolution(FRF)
end

"""
    solve(prob::StateSpaceFRFProblem; type = :dis, ismat = false, progress = true)
    solve(prob::StateSpaceModalFRFProblem; type = :dis, ismat = false, progress = true)

Computes the FRF matrix by direct or modal method

**Inputs**
* `prob`: Structure containing the problem data
* `type`: Type of FRF to compute
    * `:dis`: Admittance
    * `:vel`: Mobility
    * `:acc`: Accelerance
* `ismat`: Return the FRF matrix as a 3D array (default = false)
* `progress`: Show progress bar (default = true)

**Output**
* `sol`: FRFSolution structure
    * `u`: FRF matrix
"""
function solve(prob::StateSpaceModalFRFProblem; type = :dis, ismat = false, progress = true)

    # Initialisation
    (; css, freq, So, Se, n) = prob
    (; Ac, Bc) = css
    no = size(So, 1)
    ne = size(Se, 1)
    nstate, nu = size(Bc)
    ns = Int(nstate/2)
    nf = length(freq)

    # Compute modal information
    λ, Ψ = eigenmode(Ac, n)

    B = similar(Ψ, n, nu)
    B .= Ψ\Bc

    if type == :dis
        Ψo = So*Ψ[1:ns, :]
    elseif type == :vel
        Ψo = So*Ψ[ns+1:end, :]
    elseif type == :acc
        C = similar(Ψ, no, n)
        D = similar(Bc, no, ne)

        C .= So*Ac[ns+1:end, :]*Ψ
        D .= So*Bc[ns+1:end, :]*Se'
    end
    Be = B*Se'

    # Signature for solving type instability
    FRF = Matrix{Complex}[similar(λ, no, ne) for _ in 1:nf]
    M = Diagonal(similar(λ, n))
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(nf, desc = "State Space FRF Problem - Modal approach...", showspeed = true)
    for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        @. M[indm] = 1/(1im*ω - λ)

        if type == :dis || type == :vel
            FRF[f] .= Ψo*M*Be
        elseif type == :acc
            FRF[f] .= C*M*Be .+ D
        end
    end

    if ismat
        return StateSpaceFRFSolution(reshape(reduce(hcat, FRF), no, ne, :))
    end

    return StateSpaceFRFSolution(FRF)
end

"""
    solve(prob::StateSpaceFreqProblem; type = :dis, progress = true)
    solve(prob::StateSpaceFreqProblem; type = :dis, progress = true)

Computes the frequency response by direct or modal method

**Inputs**
* `prob`: Structure containing the problem data
* `type::Symbol`: Type of FRF to compute
    * `:dis`: Displacement
    * `:vel`: Velocity
    * `:acc`: Acceleration
* `progress`: Show progress bar

**Output**
`sol`: Solution of the problem
    * `u`: Response spectrum matrix
"""
function solve(prob::StateSpaceFreqProblem; type = :dis, progress = true)
    # Initialisation
    (; css, freq, F, So) = prob
    (; Ac, Bc) = css
    no = size(So, 1)
    nstate, nu = size(Bc)
    ns = Int(nstate/2)
    nf = length(freq)

    y::Matrix{Complex} = similar(Ac, Complex{eltype(Ac)}, no, nf)
    M = similar(Ac, Complex{eltype(Ac)}, nstate, nu)

    if type == :acc
        C = So*Ac[ns+1:end, :]
        D = So*Bc[ns+1:end, :]
    end

    ωf = 2π*freq
    p = Progress(nf, desc = "State Space Frequency Problem - Direct method...", showspeed = true)
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
    solve(prob::StateSpaceFreqProblem; type = :dis, progress = true)
    solve(prob::StateSpaceFreqProblem; type = :dis, progress = true)

Computes the frequency response by direct or modal method

**Inputs**
* `prob`: Structure containing the problem data
* `type::Symbol`: Type of FRF to compute
    * `:dis`: Displacement
    * `:vel`: Velocity
    * `:acc`: Acceleration
* `progress`: Show progress bar

**Output**
`sol`: Solution of the problem
    * `u`: Response spectrum matrix
"""
function solve(prob::StateSpaceModalFreqProblem; type = :dis, progress = true)
    # Initialisation
    (; css, F, freq, So, n) = prob
    (; Ac, Bc) = css
    no = size(So, 1)
    nstate, nu = size(Bc)
    ns = Int(nstate/2)
    nf = length(freq)

    # Compute modal information
    λ, Ψ = eigenmode(Ac, n)

    B = similar(Ψ, n, nu)
    B .= Ψ\Bc

    if type == :dis
        Ψo = So*Ψ[1:ns, :]
    elseif type == :vel
        Ψo = So*Ψ[ns+1:end, :]
    elseif type == :acc
        C = similar(Ψ, no, n)
        D = similar(Bc, no, nu)

        C .= So*Ac[ns+1:end, :]*Ψ
        D .= So*Bc[ns+1:end, :]
    end
    n = length(λ)

    y = similar(λ, no, nf)
    M = Diagonal(similar(λ, n))
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(nf, desc = "State Space Frequency Problem - Modal approach...", showspeed = true)
    for (f, (ω, Fe)) in enumerate(zip(ωf, eachcol(F)))
        progress ? next!(p) : nothing

        @. M[indm] = 1/(1im*ω - λ)

        if type == :dis || type == :vel
            y[:, f] .= Ψo*M*B*Fe
        elseif type == :acc
            y[:, f] .= (C*M*B .+ D)*Fe
        end
    end

    return StateSpaceFreqSolution(y)
end