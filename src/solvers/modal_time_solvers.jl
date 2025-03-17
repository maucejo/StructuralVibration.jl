"""
    FreeModalTimeProblem(K, M, ξn, n = size(K, 1); ismodal = false)

Structure containing data for the modal time solver

# Fields
* `K`: Stiffness matrix
* `M`: Mass matrix
* `ξn`: Damping ratios
* `n`: Number of modes to retain inf the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data
"""
@show_struct struct FreeModalTimeProblem{T <: Real, S <: AbstractMatrix, U <: AbstractVector}
    K::VecOrMat{T}
    M::S
    ξn::Vector{T}
    u0::Tuple{Vector{T}, Vector{T}}
    t::V
    n::Int
    ismodal::Bool

    function FreeModalTimeProblem(K::VecOrMat{T}, M::S, ξn::Vector{T}, u0::Tuple{Vector{T}, Vector{T}}, t::U, n = size(K, 1), ismodal = false) where {T, S, U}
        if !isa(ξn, Array)
            ξn = fill(ξn, n)
        elseif length(ξn) != n
            error("The number of damping ratios must be equal to n")
        end

        new{T, S, U}(K, M, ξn, u0, t, n, ismodal)
    end
end

"""
    HarmonicModalTimeProblem(K, M, ξn, u0, t, F, ω = 0., n = size(K, 1); ismodal = false)

Structure containing data for the modal time solver for computing the forced response due to an harmonic excitation

# Constructor
* `K`: Stiffness matrix (or modal stiffness matrix)
* `M`: Mass matrix (or mass-normalized mode shapes)
* `ξn`: Damping ratios
* `u0`: Initial conditions
    * `x0`: Initial displacement (or modal displacement)
    * `v0`: Initial velocity (or modal velocity)
* `t`: Time points at which to evaluate the response
* `F`: Amplitude vector (or modal force amplitude vector)
* `freq`: Excitation frequency
* `n`: Number of modes to retain in the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data

# Fields
* `K`: Stiffness matrix (or modal stiffness matrix)
* `M`: Mass matrix (or mass-normalized mode shapes)
* `ξn`: Damping ratios
* `u0`: Initial conditions
    * `x0`: Initial displacement (or modal displacement)
    * `v0`: Initial velocity (or modal velocity)
* `t`: Time points at which to evaluate the response
* `F`: Amplitude vector (or modal force amplitude vector)
* `ω`: Excitation angular frequency
* `n`: Number of modes to retain in the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data
"""
@show_struc struct HarmonicModalTimeProblem{T <: Real, S <: AbstractMatrix, U <: AbstractVector}
    K::Matrix{T}
    M::S
    ξn::Vector{T}
    u0::Tuple{Vector{T}, Vector{T}}
    t::U
    F::Vector{T}
    ω::T
    n::Int
    ismodal::Bool

    function HarmonicModalTimeProblem(K::Matrix{T}, M::S, ξn::Vector{T}, u0::Tuple{Vector{T}, Vector{T}}, t::U, F::Vector{T}, freq::T, n = size(K, 1); ismodal = false) where {T, S, U}
        if !isa(ξn, Array)
            ξn = fill(ξn, n)
        elseif length(ξn) != n
            error("The number of damping ratios must be equal to n")
        end

        new{T, S, U}(K, M, ξn, u0, t, F, 2π*freq, n, ismodal)
    end
end

"""
    ForcedModalTimeProblem(K, M, ξn, u0, t, F, n = size(K, 1); ismodal = false)

Structure containing data for modal time solver for computing the forced response due to an arbitrary excitation

# Fields
* `K`: Stiffness matrix (or modal stiffness matrix)
* `M`: Mass matrix (or mass-normalized mode shapes)
* `ξn`: Damping ratios
* `u0`: Initial conditions
    * `u0[1]`: Initial displacement (or modal displacement)
    * `u0[2]`: Initial velocity (or modal velocity)
* `t`: Time points at which to evaluate the response
* `F`: External force matrix (or modal force matrix)
* `n`: Number of modes to retain in the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data
"""
@show_struct struct ForcedModalTimeProblem{T <: Real, S <: AbstractMatrix, U <: AbstractVector}
    K::VecOrMat{T}
    M::S
    ξn::Vector{T}
    u0::Tuple{Vector{T}, Vector{T}}
    t::U
    F::Matrix{T}
    n::Int
    ismodal::Bool

    function ForcedModalTimeProblem(K::VecOrMat{T}, M::S, ξn::Vector{T}, u0::Tuple{Vector{T}, Vector{T}}, t::U, F::Matrix{T}, n = size(K, 1); ismodal = false) where {T, S, U}
        if !isa(ξn, Array)
            ξn = fill(ξn, n)
        elseif length(ξn) != n
            error("The number of damping ratios must be equal to n")
        end

        new{T, S, U}(K, M, ξn, u0, t, F, n, ismodal)
    end
end

"""
    ModalTimeSolution(u, du, ddu)

Structure containing problem solutions

# Fields
* `u`: Displacement matrix or vector
* `du`: Velocity matrix or vector
* `ddu`: Acceleration matrix or vector
"""
@show_struct struct ModalTimeSolution{T <: Real}
    u::Matrix{T}
    du::Matrix{T}
    ddu::Matrix{T}
end

"""
    ModalImpulseSolution(u)

Structure containing the impulse response of a modal system

# Field
* `u`: Impulse response of a Mdof system
"""
@show_struct struct ModalImpulseSolution{T <: Real}
    u::Union{Matrix{T}, Vector{Matrix{T}}}
end


"""
    solve(prob::FreeModalTimeProblem)

Compute the free response of a multi-degrees of freedom (Mdof) system using the modal approach.

# Inputs
* `prob`: Structure containing the parameters of the Mdof problem

# Output
* `sol`: ModalTimeSolution structure containing the response of the system at the given time points
"""
function solve(prob::FreeModalTimeProblem)
    (; K, M, ξn, u0, t, n, ismodal) = prob
    x0, v0 = u0
    nt = length(t)

    # Modal analysis
    if !ismodal
        λ, Φ = eigen(K, M)
        λm = λ[1:n]
        Φm = Φ[:, 1:n]
        ωn = .√λm;
        # Note: The mode shapes are mass-normalized, so Mn = I

        # Modal initial conditions
        qx = Φm'*M*x0
        qv = Φm'*M*v0
    else
        Φm = M[:, 1:n]
        if K isa Vector
            ωn = .√K[1:n]
        else
            ωn = .√diag(K)[1:n]
        end
        qx = x0[1:n]
        qv = v0[1:n]
    end

    # Modal coordinate calculation
    q = similar(x0, nt, n)
    dq = similar(q)
    ddq = similar(q)
    for (m, (ωm, ξm, Am, qm))  in enumerate(zip(ωn, ξn, qx, qv))
        if ωm == 0.
            @. q[:, m] = Am + qm*t
            @. dq[:, m] = qm
            @. ddq[:, m] = 0.
        elseif ξm < 1.
            Ωm = ωm*sqrt(1 - ξm^2)
            Ξm = ξm*ωm

            # Constant Bn
            Bm = (qm + Ξm*Am)/Ωm

            @. q[:, m] = (Am*cos(Ωm*t) + Bm*sin(Ωm*t))*exp(-Ξm*t)

            @. dq[:, m] = Ωm*(-Am*sin(Ωm*t) + Bm*cos(Ωm*t))*exp(-ξ*ωm*t) - Ξm*q[:, m]

            @. ddq[:, m] = -2Ξm*dq[:, m] - ωm^2*q[:, m]

        elseif ξm == 1.
            Bm = qm + ω*Am - ω*ρm*sin(ϕm)

            @. q[:, m] = (Am + Bm*t)*exp(-ωm*t)

            @. dq[:, m] = Bm*exp(-ωm*t) - ωm*q[:, m]

            @. ddq[:, m] = -2ωm*dq[:, m] - ωm^2*q[:, m]

        else
            βm = ωm*sqrt(ξm^2 - 1)
            Bm = (qm + Ξm*Am + ω*ρm*sin(ϕm))/βm

            @. q[:, m] = (Am*cosh(βm*t) + Bm*sinh(βm*t))*exp(-Ξm*t)

            @. dq[:, m] = βm*(Am*sinh(βm*t) + Bm*cosh(βm*t))*exp(-Ξm*t) - Ξm*q[:, m]

            @. ddq[:, m] = -2Ξm*dq[:, m] - ωm^2*q[:, m]
        end
    end

    # Computation of the displacement
    u = Φm*q';
    du = Φm*dq';
    ddu = Φm*ddq';

    return ModalTimeSolution(u, du, ddu)
end

"""
    solve(prob::ForcedModalTimeProblem)

Compute the forced response of a multi-degrees of freedom (Mdof) system due to an harmonic excitation using the modal approach.

# Inputs
* `prob`: Structure containing the parameters of the Mdof problem

# Output
* `sol`: ModalTimeSolution structure containing the response of the system at the given time points
"""
function solve(prob::HarmonicModalTimeProblem)
    (; K, M, ξn, u0, t, F, ω, n, ismodal) = prob
    x0, v0 = u0
    nt = length(t)

    if size(F, 2) ≠ 1
        error("The external force amplitude must be a vector")
    end

    if !ismodal
        # Modal analysis
        λ, Φ = eigen(K, M)
        λm = λ[1:n]
        Φm = Φ[:, 1:n]
        ωn = .√λm;
        # Note: The mode shapes are mass-normalized, so Mn = I

        # Modal initial conditions
        qx = Φm'*M*x0
        qv = Φm'*M*v0
    else
        Φm = M[:, 1:n]
        if K isa Vector
            ωn = .√K[1:n]
        else
            ωn = .√diag(K)[1:n]
        end
        qx = x0[1:n]
        qv = v0[1:n]
    end

    # Modal viscous Damping vector
    Ξ = ξn.*ωn

    # Modal participation factor
    Ln = Φm'*F

    # Particular solution
    Qp = @. Ln/(ωn^2 - ω^2 + 2im*Ξ*ω)
    ρ = abs.(Qp)
    ϕ = angle.(Qp)

    A = @. qx - ρ*cos(ϕ)

    # Modal coordinate calculation
    q = similar(x0, nt, n)
    dq = similar(q)
    ddq = similar(q)
    qh = similar(t)
    dqh = similar(t)
    ddqh = similar(t)

    for (m, (ωm, ξm, Am, qm, ρm, ϕm))  in enumerate(zip(ωn, ξn, A, qv, ρ, ϕ))
        if ξm == 0. && ωm == ω
            Am = qx[m]
            Bm = qm/ω
            ρm = Ln[m]/2ω

            @. q[:, m] = Am*cos(ω*t) + Bm*sin(ω*t) + ρm*t*sin(ω*t)
            @. dq[:, m] = -Am*ω*sin(ω*t) + Bm*ω*cos(ω*t) + ρm*(sin(ω*t) + ω*t*cos(ω*t))
            @. ddq[:, m] = -Am*ω^2*cos(ω*t) - Bm*ω^2*sin(ω*t) + ρm*(2ω*cos(ω*t) - ω^2*t*sin(ω*t))
        else
            if ωm == 0.
                @. qh = Am + qm*t
                @. dqh = qm
                @. ddqh = 0.
            elseif ξm < 1.
                Ωm = ωm*sqrt(1 - ξm^2)
                Ξm = ξm*ωm

                # Constant Bn
                Bm = (qm + Ξm*Am + ω*ρm*sin(ϕm))/Ωm

                @. qh = (Am*cos(Ωm*t) + Bm*sin(Ωm*t))*exp(-Ξm*t)
                @. dqh = Ωm*(-Am*sin(Ωm*t) + Bm*cos(Ωm*t))*exp(-Ξm*t) - Ξm*qh
                @. ddqh = -2Ξm*dqh - ωm^2*qh
            elseif ξm == 1.
                Bm = qm + ωm*Am + ρm*ω*sin(ϕm)

                @. qh = (Am + Bm*t)*exp(-ωm*t)
                @. dqh = Bm*exp(-ωm*t) - ωm*qh
                @. ddqh = -2ωm*dqh - ωm^2*qh
            else
                βm = ωm*sqrt(ξm^2 - 1)
                Bm = (qm + Ξm*Am + ω*ρm*sin(ϕm))/βm

                @. qh = (Am*cosh(βm*t) + Bm*sinh(βm*t))*exp(-Ξm*t)
                @. dqh = βm*(Am*sinh(βm*t) + Bm*cosh(βm*t))*exp(-Ξm*t) - Ξm*qh
                @. ddqh = -2Ξm*dqh - ωm^2*qh
            end

            @. q[:, m] = qh + ρm*cos(ω*t + ϕm)
            @. dq[:, m] = dqh - ρm*ω*sin(ω*t + ϕm)
            @. ddq[:, m] = ddqh - ρm*ω^2*cos(ω*t + ϕm)
        end
    end

    u = Φm*q';
    du = Φm*dq';
    ddu = Φm*ddq';

    return ModalTimeSolution(u, du, ddu)
end


"""
    solve(prob::ForcedModalTimeProblem)

Compute the forced response of a multi-degrees of freedom (Mdof) system due to an arbitrary excitation using the modal approach.

# Inputs
* `prob`: Structure containing the parameters of the Mdof problem
* `method`: Method to compute the Duhamel's integral
    * `:filt`: Filtering using the Z-transform of the impulse response (default)
    * `:interp`: Interpolation + Gaussian quadrature
    * `:conv`: Convolution

# Output
* `sol`: ModalTimeSolution structure containing the response of the system at the given time points
"""
function solve(prob::ForcedModalTimeProblem; method = :interp)
    (; K, M, ξn, u0, t, F, n, ismodal) = prob
    x0, v0 = u0
    nt = length(t)
    Δt = t[2] - t[1]

    if !ismodal
        # Modal analysis
        λ, Φ = eigen(K, M)
        λm = λ[1:n]
        Φm = Φ[:, 1:n]
        ωn = .√λm;
        # Note: The mode shapes are mass-normalized, so Mn = I

        # Modal initial conditions
        qx = Φm'*M*x0
        qv = Φm'*M*v0

        # Modal participation factor
        Ln = Φm'*F
    else
        # Mode shapes and natural frequencies are already provided
        Φm = M[:, 1:n]
        if K isa Vector
            ωn = .√K[1:n]
        else
            ωn = .√diag(K)[1:n]
        end

        # Modal initial condition
        qx = x0[1:n]
        qv = v0[1:n]

         # Modal participation factor
        Ln = F[1:n, :]
    end

    # Modal coordinate calculation
    q = similar(u0, nt, n)
    qh = similar(t)
    if method == :interp || method == :conv
        h = similar(t)
    else
        num = undefs(3)
        denom = undefs(3)
    end

    for (m, (ωm, ξm, qxm, qvm)) in enumerate(zip(ωn, ξn, qx, qv))
        if ωm == 0.
            @. qh = qxm + qvm*t
            if method == :interp || method == :conv
                h = @. t
            else
                num .= [0., Δt, 0.]
                denom .= [1., -2., 1.]
            end
        elseif ξm < 1.
            Ωm = ωm*sqrt(1 - ξm^2)
            Ξm = ξm*ωm
            Am = qxm
            Bm = (qvm + Ξm*qxm)/Ωm
            @. qh = (Am*cos(Ωm*t) + Bm*sin(Ωm*t))*exp(-Ξm*t)
            if method == :interp || method == :conv
                @. h = exp(-Ξm*t)*sin(Ωm*t)/Ωm
            else
                α = exp(-Ξm*Δt)
                β = Ωm*Δt
                # Transfer function in the z-domain
                num .= [0., α*sin(β)/Ωm, 0.]
                denom .= [1., -2*α*cos(β), α^2]
            end
        elseif ξm == 1.
            Am = qxm
            Bm = qvm + ωm*qxm
            @. qh = (Am + Bm*t)*exp(-ωm*t)
            if method == :interp || method == :conv
                @. h = t*exp(-ωm*t)
            else
                α = exp(-ωm*Δt)
                num = [0., α*Δt, 0.]
                denom = [1., -2*α, α^2]
            end
        else
            βm = ωm*sqrt(ξm^2 - 1)
            Ξm = ξm*ωm
            Am = qxm
            Bm = (qvm + Ξm*qxm)/βm
            @. qh = (Am*cosh(βm*t) + Bm*sinh(βm*t))*exp(-Ξm*t)
            if method == :interp || method == :conv
                @. h = exp(-Ξm*t)*sinh(βm*t)/βm
            else
                α = exp(-ξ*ω0*Δt)
                γ = β*Δt
                num .= [0., α*sinh(γ)/m/β, 0.]
                denom .= [1., -2*α*cosh(γ), α^2]
            end
        end

        if method == :interp
            q[:, m] .= qh .+ duhamel_integral(Ln[m, :], h, t)
        elseif method == :conv
            q[:, m] .= qh .+ Δt*DSP.conv(Ln[m, :], h)[1:nt]
        else
            q[:, m] .= qh .+ Δt*DSP.filter(num, denom, Ln[m, :])
        end
    end

    u = Φm*q';
    du = gradient(u, t);
    ddu = gradient(du, t);

    return ModalTimeSolution(u, du, ddu)
end

"""
    impulse_response(K::Matrix{Float64}, M::Matrix{Float64}, ξn, t, n = size(K, 1); ismat = false)

Compute the impulse response of a multi-degrees of freedom (Mdof) system using the modal approach

# Inputs
* `K`: Stiffness matrix
* `M`: Mass matrix
* `ξn`: Damping ratios
* `t`: Time points at which to evaluate the response
* `n`: Number of modes to retain in the modal basis
* `ismat`: Flag to indicate if the output should be a matrix

# Output
* `sol`: ModalImpulseSolution
"""
function impulse_response(K::Matrix{T}, M::AbstractMatrix{T}, ξn, t::AbstractVector, n = size(K, 1); ismat = false) where {T <: Real}
    λ, Φ = eigen(K, M)
    fn = .√λ[1:n]/2π
    Φn = Φ[:, 1:n]
    ndofs = size(K, 1)
    nt = length(t)

    if !isa(ξn, Array)
        ξn = fill(ξn, n)
    elseif length(ξn) != n
        error("The number of damping ratios must be equal to n")
    end

    h = [similar(K) for _ in 1:nt]
    hsdof = similar(t)

    for (fm, ξm, Φm) in zip(fn, ξn, eachcol(Φn))
        # The modes are mass-normalized
        sdof = Sdof(1., fm, ξm)
        prob = SdofFreeTimeProblem(sdof, [0., 1.], t)
        hsdof .= solve(prob).u

        for (i, hi) in enumerate(hsdof)
            h[i] .+= Φm*hi*Φm'
        end
    end

    if ismat
        return ModalImpulseSolution(reshape(reduce(hcat, h), ndofs, ndofs, :))
    end

    return ModalImpulseSolution(h)
end