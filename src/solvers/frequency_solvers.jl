"""
    ModalFRFProblem(ωn, ξn, freq, ϕo, ϕe)

Structure containing the data feeding the modal solver for calculating an FRF

# Fields
* ωn : Resonance frequencies
* ξn : Modal damping ratios
* freq : Frequencies of interest
* ϕe : Mode shapes at excitation points
* ϕo : Mode shapes at observation points

# note
The mode shapes must be mass-normalized
"""
@show_data struct ModalFRFProblem{T <: Real, Tf <: AbstractVector, Tp <: AbstractMatrix}
    ωn::Vector{T}
    ξn::Vector{T}
    freq::Tf
    ϕo::Tp
    ϕe::Tp

    function ModalFRFProblem(ωn::Vector{T}, ξn::Union{T, Vector{T}}, freq::Tf, ϕo::Tp, ϕe::Tp) where {T, Tf, Tp}
        if !isa(ξn, Array)
            ξn = fill(ξn, length(ωn))
        elseif length(ξn) != length(ωn)
            throw(DimensionMismatch("The number of damping ratios must be equal to the number of resonance frequencies"))
        end

        new{T, Tf, Tp}(ωn, ξn, freq, ϕo, ϕe)
    end
end

"""
    DirectFRFProblem(K, M, C, freq, So, Se)

Structure containing the data feeding the direct solver for calculating an FRF

# Fields
* K: Stiffness matrix
* M: Mass matrix
* C: Damping matrix
* freq: Frequencies of interest
* So: Selection matrix for observation points
* Se: Selection matrix for excitation points
"""
@show_data struct DirectFRFProblem{Tk <: AbstractMatrix, Tm <: AbstractMatrix, Tc <: AbstractMatrix, Tf <: AbstractVector, Ts <: AbstractMatrix}
    K::Tk
    M::Tm
    C::Tc
    freq::Tf
    So::Ts
    Se::Ts

    DirectFRFProblem(K::Tk, M::Tm, C::Tc, freq::Tf, So::Ts = I(size(K, 1)), Se::Ts = I(size(K, 1))) where {Tk, Tm, Tc, Tf, Ts} = new{Tk, Tm, Tc, Tf, Ts}(K, M, C, freq, So, Se)
end

"""
    ModalFreqProblem(ωn, ξn, Fn, freq, ϕo)

Structure containing the data feeding the modal solver for calculating the frequency reponse by modal approach

# Fields
* ωn: Resonance frequencies
* ξn: Modal damping ratios
* Fn: Modal force matrix
* freq: Frequencies of interest
* ϕo: Mode shapes at observation points
"""
@show_data struct ModalFreqProblem{T <: Real, Tfn <: AbstractMatrix, Tf <: AbstractVector}
    ωn::Vector{T}
    ξn::Vector{T}
    Fn::Tfn
    freq::Tf
    ϕo::Matrix{T}

    function ModalFreqProblem(ωn::Vector{T}, ξn::Union{T, Vector{T}}, Fn::Tfn, freq::Tf, ϕo::Matrix{T}) where {T, Tfn, Tf}
        if !isa(ξn, Array)
            ξn = fill(ξn, length(ωn))
        elseif length(ξn) != length(ωn)
            throw(DimensionMismatch("The number of damping ratios must be equal to the number of resonance frequencies"))
        end

        new{T, Tfn, Tf}(ωn, ξn, Fn, freq, ϕo)
    end
end

"""
    DirectFreqProblem(K, M, C, F, freq, So)

Structure containing the data feeding the direct solver for calculating the modal frequencies

# Fields
* K: Stiffness matrix
* M: Mass matrix
* C: Damping matrix
* F: Force matrix
* freq: Frequencies of interest
* So: Selection matrix for observation points
"""
@show_data struct DirectFreqProblem{Tk <: AbstractMatrix, Tm <: AbstractMatrix, Tc <: AbstractMatrix, TF <: AbstractMatrix, Tf <: AbstractVector, Ts <: AbstractMatrix}
    K::Tk
    M::Tm
    C::Tc
    F::TF
    freq::Tf
    So::Ts

    DirectFreqProblem(K::Tk, M::Tm, C::Tc, F::TF, freq::Tf, So::Ts = I(size(K, 1))) where {Tk, Tm, Tc , TF, Tf, Ts} = new{Tk, Tm, Tc, TF, Tf, Ts}(K, M, C, F, freq, So)
end

"""
    FRFSolution(u)

Structure containing the solution of the frequency response problem

# Fields
* u: Transfer function matrix
"""
@show_data struct FRFSolution{T <: Complex}
    u::Union{Array{T, 3}, Vector{Matrix{T}}}
end

"""
    FrequencySolution(u)

Structure containing the solution of the frequency response problem

# Fields
* u: Frequency response matrix
"""
@show_data struct FrequencySolution{T <: Complex}
    u::Matrix{T}
end

"""
    solve(m::ModalFRFProblem, type = :dis; ismat = false)

Computes the FRF matrix by modal approach

# Parameter
* `m`: Structure containing the problem data
* `type`: Type of FRF to compute (:dis, :vel, :acc)
* `ismat`: Return the FRF matrix as a 3D array (default = false)
* `progress`: Show progress bar (default = true)

# Output
* `sol`: FRFSolution structure
"""
function solve(prob::ModalFRFProblem, type = :dis; ismat = false, progress = true)
    # Initialisation
    (; ωn, ξn, freq, ϕo, ϕe) = prob
    n = length(ωn)
    ne = size(ϕe, 1)
    no = size(ϕo, 1)
    nf = length(freq)

    FRF = [similar([], Complex{eltype(ωn)}, no, ne) for _ in 1:nf]
    M = Diagonal(similar([], Complex{eltype(ωn)}, n))
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(nf, color = :black, desc = "FRF calculation - Modal approach...", showspeed = true)
    @inbounds for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        @. M[indm] = 1/(ωn^2 - ω^2 + 2im*ξn*ωn*ω)
        FRF[f] .= ϕo*M*ϕe'

        if type == :vel
            FRF[f] .*= 1im*ω
        elseif type == :acc
            FRF[f] .*= -ω^2
        end
    end

    if ismat
        return FRFSolution(reshape(reduce(hcat, FRF), no, ne, :))
    end

    return FRFSolution(FRF)
end

"""
    solve(m::DirectFRF, type = :dis; ismat = false)

Computes the FRF matrix by direct method

# Parameter
* m: Structure containing the problem data
* type: Type of FRF to compute (:dis, :vel, :acc)
* ismat: Return the FRF matrix as a 3D array (default = false)

# Output
* sol: FRFSolution structure
"""
function solve(m::DirectFRFProblem, type = :dis; ismat = false, progress = true)
    # Initialisation
    (; K, M, C, freq, So, Se) = m
    no = size(So, 1)
    ne = size(Se, 1)
    Ndofs = size(K, 1)
    nf = length(freq)

    FRF = [similar([], Complex{eltype(ωn)}, no, ne) for _ in 1:nf]
    D = similar([], Complex{eltype(ωn)}, Ndofs, Ndofs)

    ωf = 2π*freq
    p = Progress(nf, color = :black, desc = "FRF calculation - Direct method...", showspeed = true)
    @inbounds for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        D .= (K + 1im*ω*C  - ω^2*M)\I
        FRF[f] .= So*D*Se'

        if type == :vel
            FRF[f] .*= 1im*ω
        elseif type == :acc
            FRF[f] .*= -ω^2
        end
    end

    if ismat
        return FRFSolution(reshape(reduce(hcat, FRF), no, ne, :))
    end

    return FRFSolution(FRF)
end

"""
    solve(m::ModalFreqProblem, type = :dis)

Computes the frequency response by modal approach

# Inputs
* m: Structure containing the problem data
* type: Type of response to compute (:dis, :vel, :acc)

# Output
* sol: FrequencySolution structure
"""
function solve(m::ModalFreqProblem, type = :dis; progress = true)
    # Initialisation
    (; ωn, ξn, Fn, freq, ϕo) = m
    no, n = size(ϕo)
    nf = length(freq)

    ωf = 2π*freq

    y = similar([], Complex{eltype(ωn)}, no, nf)
    D = Diagonal(similar([], Complex{eltype(ωn)}, n))
    indm = diagind(M)

    p = Progress(nf, color = :black, desc = "Frequency Response - Modal approach...", showspeed = true)
    @inbounds for (f, (ω, F)) in enumerate(zip(ωf, eachcol(Fn)))
        progress ? next!(p) : nothing

        @. M[indm] = 1/(ωn^2 - ω^2 + 2im*ξn*ωn*ω)
        y[:, f] .= ϕo*M*F

        if type == :vel
            y[:, f] .*= 1im*ω
        elseif type == :acc
            y[:, f] .*= -ω^2
        end
    end

    return FrequencySolution(y)
end

"""
    solve(m::DirectFreqProblem, type = :dis)

Computes the frequency response by direct method

# Inputs
* m: Structure containing the problem data
    * K: Stiffness matrix
    * M: Mass matrix
    * C: Damping matrix
    * F: Force matrix
    * freq: Frequencies of interest
    * So: Selection matrix for observation points
* type: Type of response to compute (:dis, :vel, :acc)

# Output
* sol: FrequencySolution structure
"""
function solve(m::DirectFreqProblem, type = :dis; progress = true)
    # Initialisation
    (; K, M, C, F, freq, So) = m
    Ndofs = size(K, 1)
    nf = length(freq)

    ωf = 2π*freq

    y = similar([], Complex{eltype(ωn)}, no, nf)
    D = similar([], Complex{eltype(ωn)}, Ndofs, Ndofs)

    p = Progress(nf, color = :black, desc = "Frequency Response - Direct method...", showspeed = true)
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