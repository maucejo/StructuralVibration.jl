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
@with_kw struct ModalFRFProblem
    ωn :: Vector{Float64}
    ξn :: Vector{Float64}
    freq
    ϕo
    ϕe

    function ModalFRFProblem(ωn, ξn, freq, ϕo, ϕe)
        if !isa(ξn, Array)
            ξn = fill(ξn, length(ωn))
        elseif length(ξn) != length(ωn)
            error("The number of damping ratios must be equal to the number of resonance frequencies")
        end

        new(ωn, ξn, freq, ϕo, ϕe)
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
@with_kw struct DirectFRFProblem
    K
    M
    C
    freq
    So
    Se

    DirectFRFProblem(K, M, C, freq, So = I(size(K, 1)), Se = I(size(K, 1))) = new(K, M, C, freq, So, Se)
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
@with_kw struct ModalFreqProblem
    ωn :: Vector{Float64}
    ξn :: Vector{Float64}
    Fn
    freq
    ϕo

    function ModalFreqProblem(ωn, ξn, Fn, freq, ϕo)
        if !isa(ξn, Array)
            ξn = fill(ξn, length(ωn))
        elseif length(ξn) != length(ωn)
            error("The number of damping ratios must be equal to the number of resonance frequencies")
        end

        new(ωn, ξn, Fn, freq, ϕo)
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
@with_kw struct DirectFreqProblem
    K
    M
    C
    F
    freq
    So

    DirectFreqProblem(K, M, C, F, freq, So = I(size(K, 1))) = new(K, M, C, F, freq, So)
end

"""
    FRFSolution(u)

Structure containing the solution of the frequency response problem

# Fields
* u: Transfer function matrix
"""
@with_kw struct FRFSolution
    u :: Union{Array{ComplexF64, 3}, Vector{Matrix{ComplexF64}}}
end

"""
    FrequencySolution(u)

Structure containing the solution of the frequency response problem

# Fields
* u: Frequency response matrix
"""
@with_kw struct FrequencySolution
    u :: Matrix{ComplexF64}
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

    FRF = [undefs(ComplexF64, no, ne) for _ in 1:nf]
    M = Diagonal(undefs(ComplexF64, n))
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

    FRF = [undefs(ComplexF64, no, ne) for _ in 1:nf]
    D = undefs(ComplexF64, Ndofs, Ndofs)

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

    y = undefs(ComplexF64, no, nf)
    M = Diagonal(undefs(ComplexF64, n))
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

    y = undefs(ComplexF64, no, nf)
    D = undefs(ComplexF64, Ndofs, Ndofs)

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