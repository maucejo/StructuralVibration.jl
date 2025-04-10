"""
    FreeModalTimeProblem(K, M, ξn, n = size(K, 1); ismodal = false)

Structure containing data for the modal time solver

**Fields**
* `K::VecOrMat{Real}`:
    * If `ismodal = false` then `Matrix{Real}`: Stiffness matrix
    * If `ismodal = true` then `Vector{Real}`: Natural angular frequencies
* `M::AbtractMatrix`:
    * If `ismodal = false`: Mass matrix
    * If `ismodal = true`: Mass-normalized mode shapes
* `ξn`: Damping ratios
* `u0::Tuple`: Initial conditions
    * `u0[1]`: Initial displacement (or modal displacement)
    * `u0[2]`: Initial velocity (or modal velocity)
* `t::AbstractRange`: Time points at which to evaluate the response
* `F::Matrix{Real}`: External force matrix (or modal force matrix)
* `n::Int`: Number of modes to retain in the modal basis
* `ismodal::Bool`: Flag to indicate if the problem contains modal data
"""
@show_data struct FreeModalTimeProblem{Tk <: Real, Tm <: AbstractMatrix, Tx <:Real, Tu <: AbstractVector, Tt <: AbstractRange}
    K::VecOrMat{Tk}
    M::Tm
    ξn::Vector{Tx}
    u0::Tuple{Tu, Tu}
    t::Tt
    n::Int
    ismodal::Bool

    function FreeModalTimeProblem(K::VecOrMat{Tk}, M::Tm, ξn::Union{Tx, Vector{Tx}}, u0::Tuple{Tu, Tu}, t::Tt, n = size(K, 1); ismodal = false) where {Tk, Tm, Tx, Tu, Tt}
        if !isa(ξn, Array)
            ξn = fill(ξn, n)
        elseif length(ξn) != n
            error("The number of damping ratios must be equal to n")
        end

        new{Tk, Tm, Tx, Tu, Tt}(K, M, ξn, u0, t, n, ismodal)
    end
end

"""
    HarmonicModalTimeProblem(K, M, ξn, u0, t, F, ω = 0., n = size(K, 1); ismodal = false)

Structure containing data for the modal time solver for computing the forced response due to an harmonic excitation

**Constructor**
* `K::VecOrMat{Real}`:
    * If `ismodal = false` then `Matrix{Real}`: Stiffness matrix
    * If `ismodal = true` then `Vector{Real}`: Natural angular frequencies
* `M::AbtractMatrix`:
    * If `ismodal = false`: Mass matrix
    * If `ismodal = true`: Mass-normalized mode shapes
* `ξn`: Damping ratios
* `u0::Tuple`: Initial conditions
    * `u0[1]`: Initial displacement (or modal displacement)
    * `u0[2]`: Initial velocity (or modal velocity)
* `t::AbstractRange`: Time points at which to evaluate the response
* `F::Real`: External force matrix (or modal participation factors)
* `freq::Real`: Excitation frequency
* `n::Int`: Number of modes to retain in the modal basis
* `ismodal::Bool`: Flag to indicate if the problem contains modal data

**Fields**
* `K`: Stiffness matrix (or modal stiffness matrix)
* `M`: Mass matrix (or mass-normalized mode shapes)
* `ξn`: Damping ratios
* `u0`: Initial conditions
    * `u0[1]`: Initial displacement (or modal displacement)
    * `u0[2]`: Initial velocity (or modal velocity)
* `t`: Time points at which to evaluate the response
* `F`: Amplitude vector (or modal participation vector)
* `ω`: Excitation angular frequency
* `n`: Number of modes to retain in the modal basis
* `ismodal`: Flag to indicate if the problem contains modal data
"""
@show_data struct HarmonicModalTimeProblem{Tk <: Real, Tm <: AbstractMatrix, Tx <: Real, TF <: AbstractVector, Tf <: Real, Tu <: AbstractVector, Tt <: AbstractRange}
    K::VecOrMat{Tk}
    M::Tm
    ξn::Vector{Tx}
    F::TF
    ω::Tf
    u0::Tuple{Tu, Tu}
    t::Tt
    n::Int
    ismodal::Bool

    function HarmonicModalTimeProblem(K::VecOrMat{Tk}, M::Tm, ξn::Union{Tx, Vector{Tx}}, F::TF, freq::Tf, u0::Tuple{Tu, Tu}, t::Tt, n = size(K, 1); ismodal = false) where {Tk, Tm, Tx, TF, Tf, Tu, Tt}
        if !isa(ξn, Array)
            ξn = fill(ξn, n)
        elseif length(ξn) != n
            error("The number of damping ratios must be equal to n")
        end

        new{Tk, Tm, Tx, TF, Tf, Tu, Tt}(K, M, ξn, F, 2π*freq, u0, t, n, ismodal)
    end
end

"""
    ForcedModalTimeProblem(K, M, ξn, u0, t, F, n = size(K, 1); ismodal = false)

Structure containing data for modal time solver for computing the forced response due to an arbitrary excitation

**Fields**
* `K::VecOrMat{Real}`:
    * If `ismodal = false` then `Matrix{Real}`: Stiffness matrix
    * If `ismodal = true` then `Vector{Real}`: Natural angular frequencies
* `M::AbtractMatrix`:
    * If `ismodal = false`: Mass matrix
    * If `ismodal = true`: Mass-normalized mode shapes
* `ξn`: Damping ratios
* `F::Matrix{Real}`: External force matrix (or modal participation factors)
* `u0::Tuple`: Initial conditions
    * `u0[1]`: Initial displacement (or modal displacement)
    * `u0[2]`: Initial velocity (or modal velocity)
* `t::AbstractRange`: Time points at which to evaluate the response
* `n::Int`: Number of modes to retain in the modal basis
* `ismodal::Bool`: Flag to indicate if the problem contains modal data
"""
@show_data struct ForcedModalTimeProblem{Tk <: Real, Tm <: AbstractMatrix, Tx <: Real, Tu <: AbstractVector, Tt <: AbstractRange}
    K::VecOrMat{Tk}
    M::Tm
    ξn::Vector{Tx}
    F::Matrix{Tk}
    u0::Tuple{Tu, Tu}
    t::Tt
    n::Int
    ismodal::Bool

    function ForcedModalTimeProblem(K::VecOrMat{Tk}, M::Tm, ξn::Union{Tx, Vector{Tx}}, F::Matrix{Tk}, u0::Tuple{Tu, Tu}, t::Tt, n = size(K, 1); ismodal = false) where {Tk, Tm, Tx, Tu, Tt}
        if !isa(ξn, Array)
            ξn = fill(ξn, n)
        elseif length(ξn) != n
            error("The number of damping ratios must be equal to n")
        end

        new{Tk, Tm, Tx, Tu, Tt}(K, M, ξn, F, u0, t, n, ismodal)
    end
end

"""
    ModalTimeSolution(u, du, ddu)

Structure containing problem solutions

**Fields**
* `u::Matrix{Real}`: Displacement matrix
* `du::Matrix{Real}`: Velocity matrix
* `ddu::Matrix{Real}`: Acceleration matrix
"""
@show_data struct ModalTimeSolution{T <: Real}
    u::Matrix{T}
    du::Matrix{T}
    ddu::Matrix{T}
end

"""
    ModalImpulseSolution(u)

Structure containing the impulse response of a multi-degrees of freedom (Mdof) system

**Fields**
* `u`: Impulse response matrix
"""
@show_data struct ModalImpulseSolution{T <: Real}
    u::Union{Array{T, 3}, Vector{Matrix{T}}}
end


"""
    solve(prob::FreeModalTimeProblem)

Compute the free response of a multi-degrees of freedom (Mdof) system using the modal approach.

**Inputs**
* `prob`: Structure containing the parameters of the Mdof problem

**Output**
* `sol`: ModalTimeSolution structure containing the response of the system at the given time points
"""
function solve(prob::FreeModalTimeProblem)
    (; K, M, ξn, u0, t, n, ismodal) = prob
    x0, v0 = u0

    # Modal analysis
    if !ismodal
        ωn, Φn = eigenmode(K, M, n)
        # Note: The mode shapes are mass-normalized, so Mn = I

        # Modal initial conditions
        qx = Φn'*M*x0
        qv = Φn'*M*v0
    else
        Φn = M[:, 1:n]
        if K isa Vector
            ωn = K[1:n]
        else
            ωn = diag(K)[1:n]
        end
        qx = x0[1:n]
        qv = v0[1:n]
    end

    # Modal coordinate calculation
    q = similar(qx, length(t), n)
    dq = similar(q)
    ddq = similar(q)

    cache = (
        xh = similar(t),
        vh = similar(t),
        ah = similar(t)
    )

    for (m, (ωm, ξm, qxm, qvm))  in enumerate(zip(ωn, ξn, qx, qv))
        free_response_sdof!(cache, ωm, ξm, qxm, qvm, t)
        q[:, m] .= cache.xh
        dq[:, m] .= cache.vh
        ddq[:, m] .= cache.ah
    end

    # Computation of the displacement
    u = Φn*q'
    du = Φn*dq'
    ddu = Φn*ddq'

    return ModalTimeSolution(u, du, ddu)
end

"""
    solve(prob::HarmonicModalTimeProblem)

Compute the forced response of a multi-degrees of freedom (Mdof) system due to an harmonic excitation using the modal approach.

**Inputs**
* `prob`: Structure containing the parameters of the Mdof problem

**Output**
* `sol`: Solution structure containing the response of the system at the given time points
    * `u`: Displacement
    * `du`: Velocity
    * `ddu`: Acceleration
"""
function solve(prob::HarmonicModalTimeProblem)
    (; K, M, ξn, F, ω, u0, t, n, ismodal) = prob
    x0, v0 = u0

    if size(F, 2) ≠ 1
        error("The external force amplitude must be a vector")
    end

    if !ismodal
        # Modal analysis
        ωn, Φn = eigenmode(K, M, n)
        # Note: The mode shapes are mass-normalized, so Mn = I

        # Modal initial conditions
        qx = Φn'*M*x0
        qv = Φn'*M*v0

        # Modal participation factor
        Ln = Φn'*F
    else
        Φn = M[:, 1:n]
        if K isa Vector
            ωn = K[1:n]
        else
            ωn = diag(K)[1:n]
        end
        qx = x0[1:n]
        qv = v0[1:n]

        # Modal participation factor
        Ln = F[1:n]
    end

    # Modal coordinate calculation
    q = similar(qx, length(t), n)
    dq = similar(q)
    ddq = similar(q)

    cache = (
        x = similar(t),
        v = similar(t),
        a = similar(t),
        xh = similar(t),
        vh = similar(t),
        ah = similar(t)
    )

    for (m, (ωm, ξm, Lm,  qxm, qvm))  in enumerate(zip(ωn, ξn, Ln, qx, qv))
        harmonic_response_sdof!(cache, Lm, ω, qxm, qvm, 1., ωm, ξm, t, :force)
        q[:, m] .= cache.x
        dq[:, m] .= cache.v
        ddq[:, m] .= cache.a
    end

    u = Φn*q';
    du = Φn*dq';
    ddu = Φn*ddq';

    return ModalTimeSolution(u, du, ddu)
end

"""
    solve(prob::ForcedModalTimeProblem)

Compute the forced response of a multi-degrees of freedom (Mdof) system due to an arbitrary excitation using the modal approach.

**Inputs**
* `prob`: Structure containing the parameters of the Mdof problem
* `method`: Method to compute the Duhamel's integral
    * `:filt`: Filtering using the Z-transform of the impulse response (default)
    * `:interp`: Interpolation + Gaussian quadrature
    * `:conv`: Convolution

**Output**
* `sol`: ModalTimeSolution structure containing the response of the system at the given time points
"""
function solve(prob::ForcedModalTimeProblem; method = :filt)

    (; K, M, ξn, F, u0, t, n, ismodal) = prob
    x0, v0 = u0

    if !ismodal
        # Modal analysis
        ωn, Φn = eigenmode(K, M, n)

        # Note: The mode shapes are mass-normalized, so Mn = I

        # Modal initial conditions
        qx = Φn'*M*x0
        qv = Φn'*M*v0

        # Modal participation factor
        Ln = Φn'*F
    else
        # Mode shapes and natural frequencies are already provided
        Φn = M[:, 1:n]
        if K isa Vector
            ωn = K[1:n]
        else
            ωn = diag(K)[1:n]
        end

        # Modal initial condition
        qx = x0[1:n]
        qv = v0[1:n]

         # Modal participation factor
        Ln = F[1:n, :]
    end

    # Modal coordinate calculation
    q = similar(x0, length(t), n)
    dq = similar(q)
    ddq = similar(q)

    # Cache for the free response
    cache_h = (
        xh = similar(t),
        vh = similar(t),
        ah = similar(t),
    )

    # Cache for the particular solution
    cache_p = (
        xp = similar(t),
        vp = similar(t),
        ap = similar(t),
        h = similar(t),
        num = similar(t, 3),
        denom = similar(t, 3)
    )

    for (m, (ωm, ξm, Lm, qxm, qvm)) in enumerate(zip(ωn, ξn, eachrow(Ln), qx, qv))

        free_response_sdof!(cache_h, ωm, ξm, qxm, qvm, t)
        forced_response_sdof!(cache_p, Lm, 1., ωm, ξm, t, :force, method)

        q[:, m] .= cache_h.xh .+ cache_p.xp
        dq[:, m] .= cache_h.vh .+ cache_p.vp
        ddq[:, m] .= cache_h.ah .+ cache_p.ap
    end

    u = Φn*q'
    du = Φn*dq'
    ddu = Φn*ddq'

    return ModalTimeSolution(u, du, ddu)
end

"""
    impulse_response(K::Matrix{Float64}, M::Matrix{Float64}, ξn, t, n = size(K, 1); ismat = false)

Compute the impulse response of a multi-degrees of freedom (Mdof) system using the modal approach

**Inputs**
* `K`:
    * If `ismodal = false`: Stiffness matrix
    * If `ismodal = true`: Squared natural angular frequencies
* `M`:
    * If `ismodal = false`: Mass matrix
    * If `ismodal = true`: Mass-normalized mode shapes
* `ξn`: Damping ratios
* `t`: Time points at which to evaluate the response
* `n`: Number of modes to retain in the modal basis
* `ismodal::Bool`: Flag to indicate if the problem contains modal data
* `ismat::Bool`: Flag to indicate if the output should be a matrix

**Output**
* `sol`: ModalImpulseSolution
    * `u`: Impulse response matrix
"""
function impulse_response(K, M, ξn, t, n = size(K, 1); ismodal = false, ismat = false)
    if !ismodal
        ωn, Φn = eigenmode(K, M, n)
        fn = ωn/2π
        ndofs = size(K, 1)
    else
        if K isa Vector
            fn = K[1:n]/2π
        else
            fn = diag(K)[1:n]/2π
        end
        Φn = M[:, 1:n]
        ndofs = size(Φn, 1)
    end

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
        hsdof .= impulse_response(sdof, t)

        for (i, hi) in enumerate(hsdof)
            h[i] .= Φm*hi*Φm'
        end
    end

    if ismat
        return ModalImpulseSolution(reshape(reduce(hcat, h), ndofs, ndofs, :))
    end

    return ModalImpulseSolution(h)
end