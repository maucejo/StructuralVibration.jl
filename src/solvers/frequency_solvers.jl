# Internal helpers to allow scalar frequencies (and vector forces for a single frequency)
_freq_vector(freq::Real) = [freq]
_freq_vector(freq::AbstractVector{<:Real}) = freq

_force_matrix(F::AbstractVector) = reshape(F, :, 1)
_force_matrix(F::AbstractMatrix) = F

"""
    ModalFRFProblem(ωn, ξn, freq, ϕo, ϕe)

Structure containing the data feeding the modal solver for calculating an FRF

**Fields**
* `ωn::Vector{Real}`: Natural angular frequencies
* `ξn::Vector{Real}`: Modal damping ratios
* `freq::Union{Real, AbstractVector}`: Frequencies of interest (a scalar is accepted for a single frequency)
* `ϕe::AbstractMatrix`: Mode shapes at excitation points
* `ϕo::AbstractMatrix`: Mode shapes at observation points

**Note**
The mode shapes must be mass-normalized
"""
@show_data struct ModalFRFProblem{T <: Real, Tf <: AbstractVector{<:Real}, Tpo <: AbstractMatrix, Tpe <: AbstractMatrix}
    ωn::Vector{T}
    ξn::Vector{T}
    freq::Tf
    ϕo::Tpo
    ϕe::Tpe

    function ModalFRFProblem(ωn::Vector{T}, ξn::Union{T, Vector{T}}, freq::Union{Real, AbstractVector{<:Real}}, ϕo::Tpo, ϕe::Tpe) where {T, Tpo, Tpe}
        freq = _freq_vector(freq)
        Tf = typeof(freq)
        if !isa(ξn, Array)
            ξn = fill(ξn, length(ωn))
        elseif length(ξn) != length(ωn)
            throw(DimensionMismatch("The number of damping ratios must be equal to the number of resonance frequencies"))
        end

        new{T, Tf, Tpo, Tpe}(ωn, ξn, freq, ϕo, ϕe)
    end
end

"""
    DirectFRFProblem(K, M, C, freq, So, Se)

Structure containing the data feeding the direct solver for calculating an FRF

**Fields**
* `K::AbstractMatrix`: Stiffness matrix
* `M::AbstractMatrix`: Mass matrix
* `C::AbstractMatrix`: Damping matrix
* `freq::Union{Real, AbstractVector}`: Frequencies of interest (a scalar is accepted for a single frequency)
* `So::AbstractMatrix`: Selection matrix for observation points
* `Se::AbstractMatrix`: Selection matrix for excitation points
"""
@show_data struct DirectFRFProblem{Tk <: AbstractMatrix, Tm <: AbstractMatrix, Tc <: AbstractMatrix, Tf <: AbstractVector{<:Real}, Tso <: AbstractMatrix, Tse <: AbstractMatrix}
    K::Tk
    M::Tm
    C::Tc
    freq::Tf
    So::Tso
    Se::Tse

    function DirectFRFProblem(K::Tk, M::Tm, C::Tc, freq::Union{Real, AbstractVector{<:Real}}, So::Tso = I(size(K, 1)), Se::Tse = I(size(K, 1))) where {Tk, Tm, Tc, Tso, Tse}
        freqv = _freq_vector(freq)
        return new{Tk, Tm, Tc, typeof(freqv), Tso, Tse}(K, M, C, freqv, So, Se)
    end
end

"""
    ModalFreqProblem(ωn, ξn, Fn, freq, ϕo)

Structure containing the data feeding the modal solver for calculating the frequency response by modal approach

**Fields**
* `ωn::Vector{Real}`: Natural angular frequencies
* ξn::Vector{Real}: Modal damping ratios
* `Ln::AbstractVecOrMat`: Modal participation factors (a vector is accepted for a single frequency)
* `freq::Union{Real, AbstractVector}`: Frequencies of interest (a scalar is accepted for a single frequency)
* `ϕo::Matrix{Real}`: Mode shapes at observation points
"""
@show_data struct ModalFreqProblem{T <: Real, Tln <: AbstractMatrix, Tf <: AbstractVector{<:Real}}
    ωn::Vector{T}
    ξn::Vector{T}
    Ln::Tln
    freq::Tf
    ϕo::Matrix{T}

    function ModalFreqProblem(ωn::Vector{T}, ξn::Union{T, Vector{T}}, Ln::AbstractVecOrMat, freq::Union{Real, AbstractVector{<:Real}}, ϕo::Matrix{T}) where T
        freq = _freq_vector(freq)
        Ln = _force_matrix(Ln)
        Tf = typeof(freq)
        Tln = typeof(Ln)
        if !isa(ξn, Array)
            ξn = fill(ξn, length(ωn))
        elseif length(ξn) != length(ωn)
            throw(DimensionMismatch("The number of damping ratios must be equal to the number of resonance frequencies"))
        end

        new{T, Tln, Tf}(ωn, ξn, Ln, freq, ϕo)
    end
end

"""
    DirectFreqProblem(K, M, C, F, freq, So)

Structure containing the data feeding the direct solver for calculating the modal frequencies

**Fields**
* `K::AbstractMatrix`: Stiffness matrix
* `M::AbstractMatrix`: Mass matrix
* `C::AbstractMatrix`: Damping matrix
* `F::AbstractVecOrMat`: Force matrix (a vector is accepted for a single frequency)
* `freq::Union{Real, AbstractVector}`: Frequencies of interest (a scalar is accepted for a single frequency)
* `So::AbstractMatrix`: Selection matrix for observation points
"""
@show_data struct DirectFreqProblem{Tk <: AbstractMatrix, Tm <: AbstractMatrix, Tc <: AbstractMatrix, TF <: AbstractMatrix, Tf <: AbstractVector{<:Real}, Ts <: AbstractMatrix}
    K::Tk
    M::Tm
    C::Tc
    F::TF
    freq::Tf
    So::Ts

    function DirectFreqProblem(K::Tk, M::Tm, C::Tc, F::AbstractVecOrMat, freq::Union{Real, AbstractVector{<:Real}}, So::Ts = I(size(K, 1))) where {Tk, Tm, Tc, Ts}
        freqv = _freq_vector(freq)
        Fm = _force_matrix(F)
        return new{Tk, Tm, Tc, typeof(Fm), typeof(freqv), Ts}(K, M, C, Fm, freqv, So)
    end
end

"""
    FRFSolution(u)

Structure containing the solution of the frequency response problem

**Fields**
* `u`: Transfer function matrix
"""
@show_data struct FRFSolution{T <: Complex}
    u::Array{T, 3}
end

"""
    FrequencySolution(u)

Structure containing the solution of the frequency response problem

**Fields**
* `u`: Frequency response matrix
"""
@show_data struct FrequencySolution{T <: Complex}
    u::Matrix{T}
end

"""
    solve(prob::DirectFRFProblem; type = :dis, progress = true)
    solve(prob::ModalFRFProblem; type = :dis, progress = true)

Computes the FRF matrix by direct or modal approach

**Inputs**
* `prob`: Structure containing the problem data
* `type`: Type of FRF to compute
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

**Output**
* `sol`: Solution of the problem
    * `u`: FRF matrix
"""
function solve(prob::ModalFRFProblem; type = :dis, progress = true)
    # Initialisation
    (; ωn, ξn, freq, ϕo, ϕe) = prob
    n = length(ωn)
    ne = size(ϕe, 1)
    no = size(ϕo, 1)
    nf = length(freq)

    FRF = similar(ωn, Complex{eltype(ωn)}, no, ne, nf)
    M = Diagonal(similar(ωn, Complex{eltype(ωn)}, n))
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(nf, desc = "FRF calculation - Modal approach...", showspeed = true)
    @inbounds for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        @. M[indm] = 1/(ωn^2 - ω^2 + 2im*ξn*ωn*ω)
        FRF[:, :, f] .= ϕo*M*ϕe'

        if type == :vel
            FRF[:, :, f] .*= 1im*ω
        elseif type == :acc
            FRF[:, :, f] .*= -ω^2
        end
    end

    return FRFSolution(FRF)
end

function solve(prob::DirectFRFProblem; type = :dis, progress = true)

    # Initialisation
    (; K, M, C, freq, So, Se) = prob
    no = size(So, 1)
    ne = size(Se, 1)
    Ndofs = size(K, 1)
    nf = length(freq)

    FRF = similar(K, Complex{eltype(K)}, no, ne, nf)
    D = similar(K, Complex{eltype(K)}, Ndofs, ne)

    ωf = 2π*freq
    p = Progress(nf, desc = "FRF calculation - Direct method...", showspeed = true)
    @inbounds for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        D .= (K + 1im*ω*C  - ω^2*M)\Se'
        FRF[:, :,f] .= So*D

        if type == :vel
            FRF[:, :, f] .*= 1im*ω
        elseif type == :acc
            FRF[:, :, f] .*= -ω^2
        end
    end

    return FRFSolution(FRF)
end

"""
    solve(prob::DirectFreqProblem; type = :dis, progress = true)
    solve(prob::ModalFreqProblem; type = :dis, progress = true)

Computes the frequency response by direct or modal approach

**Inputs**
* `prob`: Structure containing the problem data
* `type`: Type of response to compute
    * `:dis`: Displacement (default)
    * `:vel`: Velocity
    * `:acc`: Acceleration
* `progress`: Show progress bar (default = true)

**Output**
* `sol`: Solution of the problem
    * `u`: Response spectrum matrix
"""
function solve(prob::ModalFreqProblem; type = :dis, progress = true)
    # Initialisation
    (; ωn, ξn, Ln, freq, ϕo) = prob
    no, n = size(ϕo)
    nf = length(freq)

    ωf = 2π*freq

    y = similar(ωn, Complex{eltype(ωn)}, no, nf)
    M = Diagonal(similar(ωn, Complex{eltype(ωn)}, n))
    indm = diagind(M)

    p = Progress(nf, desc = "Frequency Response - Modal approach...", showspeed = true)
    @inbounds for (f, (ω, L)) in enumerate(zip(ωf, eachcol(Ln)))
        progress ? next!(p) : nothing

        @. M[indm] = 1/(ωn^2 - ω^2 + 2im*ξn*ωn*ω)
        y[:, f] .= ϕo*M*L

        if type == :vel
            y[:, f] .*= 1im*ω
        elseif type == :acc
            y[:, f] .*= -ω^2
        end
    end

    return FrequencySolution(y)
end

"""
    solve(prob::DirectFreqProblem, type = :dis, progress = true)
    solve(prob::ModalFreqProblem, type = :dis, progress = true)

Computes the frequency response by direct or modal approach

**Inputs**
* `prob`: Structure containing the problem data
* `type`: Type of response to compute
    * `:dis`: Displacement (default)
    * `:vel`: Velocity
    * `:acc`: Acceleration
* `progress`: Show progress bar (default = true)

**Output**
* `sol`: Solution of the problem
    * `u`: Response spectrum matrix
"""
function solve(prob::DirectFreqProblem; type = :dis, progress = true)
    # Initialisation
    (; K, M, C, F, freq, So) = prob
    Ndofs = size(K, 1)
    no = size(So, 1)
    nf = length(freq)

    ωf = 2π*freq

    y = similar(K, Complex{eltype(K)}, no, nf)
    D = similar(K, Complex{eltype(K)}, Ndofs)

    p = Progress(nf, desc = "Frequency Response - Direct method...", showspeed = true)
    @inbounds for (f, (ω, Fe)) in enumerate(zip(ωf, eachcol(F)))
        progress ? next!(p) : nothing

        D .= (K + 1im*ω*C  - ω^2*M)\Fe
        y[:, f] .= So*D

        if type == :vel
            y[:, f] .*= 1im*ω
        elseif type == :acc
            y[:, f] .*= -ω^2
        end
    end

    return FrequencySolution(y)
end