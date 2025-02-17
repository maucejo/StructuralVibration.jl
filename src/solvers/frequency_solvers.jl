"""
    ModalFRFProblem(ωₙ, ξₙ, freq, ϕₒ, ϕₑ)

Structure containing the data feeding the modal solver for calculating an FRF

# Fields
* ωₙ : Resonance frequencies
* ξₙ : Modal damping ratios
* freq : Frequencies of interest
* ϕₑ : Mode shapes at excitation points
* ϕₒ : Mode shapes at observation points

# Note
The mode shapes must be mass-normalized
"""
@with_kw struct ModalFRFProblem
    ωₙ :: Vector{Float64}
    ξₙ :: Vector{Float64}
    freq
    ϕₒ
    ϕₑ

    function ModalFRFProblem(ωₙ, ξₙ, freq, ϕₒ, ϕₑ)
        if !isa(ξₙ, Array)
            ξₙ = fill(ξₙ, length(ωₙ))
        elseif length(ξₙ) != length(ωₙ)
            error("The number of damping ratios must be equal to the number of resonance frequencies")
        end

        new(ωₙ, ξₙ, freq, ϕₒ, ϕₑ)
    end
end

"""
    DirectFRFProblem(K, M, C, freq, Sₒ, Sₑ)

Structure containing the data feeding the direct solver for calculating an FRF

# Fields
* K: Stiffness matrix
* M: Mass matrix
* C: Damping matrix
* freq: Frequencies of interest
* Sₒ: Selection matrix for observation points
* Sₑ: Selection matrix for excitation points
"""
@with_kw struct DirectFRFProblem
    K
    M
    C
    freq
    Sₒ
    Sₑ

    DirectFRFProblem(K, M, C, freq, Sₒ = I(size(K, 1)), Sₑ = I(size(K, 1))) = new(K, M, C, freq, Sₒ, Sₑ)
end

"""
    ModalFreqProblem(ωₙ, ξₙ, Fₙ, freq, ϕₒ)

Structure containing the data feeding the modal solver for calculating the frequency reponse by modal approach

# Fields
* ωₙ: Resonance frequencies
* ξₙ: Modal damping ratios
* Fₙ: Modal force matrix
* freq: Frequencies of interest
* ϕₒ: Mode shapes at observation points
"""
@with_kw struct ModalFreqProblem
    ωₙ :: Vector{Float64}
    ξₙ :: Vector{Float64}
    Fₙ
    freq
    ϕₒ

    function ModalFreqProblem(ωₙ, ξₙ, Fₙ, freq, ϕₒ)
        if !isa(ξₙ, Array)
            ξₙ = fill(ξₙ, length(ωₙ))
        elseif length(ξₙ) != length(ωₙ)
            error("The number of damping ratios must be equal to the number of resonance frequencies")
        end

        new(ωₙ, ξₙ, Fₙ, freq, ϕₒ)
    end
end

"""
    DirectFreqProblem(K, M, C, F, freq, Sₒ)

Structure containing the data feeding the direct solver for calculating the modal frequencies

# Fields
* K: Stiffness matrix
* M: Mass matrix
* C: Damping matrix
* F: Force matrix
* freq: Frequencies of interest
* Sₒ: Selection matrix for observation points
"""
@with_kw struct DirectFreqProblem
    K
    M
    C
    F
    freq
    Sₒ

    DirectFreqProblem(K, M, C, F, freq, Sₒ = I(size(K, 1))) = new(K, M, C, F, freq, Sₒ)
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
    (; ωₙ, ξₙ, freq, ϕₒ, ϕₑ) = prob
    Nₘ = length(ωₙ)
    Nₑ = size(ϕₑ, 1)
    Nₒ = size(ϕₒ, 1)
    Nf = length(freq)

    FRF = [undefs(ComplexF64, Nₒ, Nₑ) for _ in 1:Nf]
    M = Diagonal(undefs(ComplexF64, Nₘ))
    indm = diagind(M)

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "FRF calculation - Modal approach...", showspeed = true)
    @inbounds for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        @. M[indm] = 1/(ωₙ^2 - ω^2 + 2im*ξₙ*ωₙ*ω)
        FRF[f] .= ϕₒ*M*ϕₑ'

        if type == :vel
            FRF[f] .*= 1im*ω
        elseif type == :acc
            FRF[f] .*= -ω^2
        end
    end

    if ismat
        return FRFSolution(reshape(reduce(hcat, FRF), Nₒ, Nₑ, :))
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
    (; K, M, C, freq, Sₒ, Sₑ) = m
    Nₒ = size(Sₒ, 1)
    Nₑ = size(Sₑ, 1)
    Ndofs = size(K, 1)
    Nf = length(freq)

    FRF = [undefs(ComplexF64, Nₒ, Nₑ) for _ in 1:Nf]
    D = undefs(ComplexF64, Ndofs, Ndofs)

    ωf = 2π*freq
    p = Progress(Nf, color = :black, desc = "FRF calculation - Direct method...", showspeed = true)
    @inbounds for (f, ω) in enumerate(ωf)
        progress ? next!(p) : nothing

        D .= (K + 1im*ω*C  - ω^2*M)\I
        FRF[f] .= Sₒ*D*Sₑ'

        if type == :vel
            FRF[f] .*= 1im*ω
        elseif type == :acc
            FRF[f] .*= -ω^2
        end
    end

    if ismat
        return FRFSolution(reshape(reduce(hcat, FRF), Nₒ, Nₑ, :))
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
    (; ωₙ, ξₙ, Fₙ, freq, ϕₒ) = m
    Nₒ, Nₘ = size(ϕₒ)
    Nf = length(freq)

    ωf = 2π*freq

    y = undefs(ComplexF64, Nₒ, Nf)
    M = Diagonal(undefs(ComplexF64, Nₘ))
    indm = diagind(M)

    p = Progress(Nf, color = :black, desc = "Frequency Response - Modal approach...", showspeed = true)
    @inbounds for (f, (ω, F)) in enumerate(zip(ωf, eachcol(Fₙ)))
        progress ? next!(p) : nothing

        @. M[indm] = 1/(ωₙ^2 - ω^2 + 2im*ξₙ*ωₙ*ω)
        y[:, f] .= ϕₒ*M*F

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
    * Sₒ: Selection matrix for observation points
* type: Type of response to compute (:dis, :vel, :acc)

# Output
* sol: FrequencySolution structure
"""
function solve(m::DirectFreqProblem, type = :dis; progress = true)
    # Initialisation
    (; K, M, C, F, freq, Sₒ) = m
    Ndofs = size(K, 1)
    Nf = length(freq)

    ωf = 2π*freq

    y = undefs(ComplexF64, Nₒ, Nf)
    D = undefs(ComplexF64, Ndofs, Ndofs)

    p = Progress(Nf, color = :black, desc = "Frequency Response - Direct method...", showspeed = true)
    @inbounds for (f, (ω, Fₑ)) in enumerate(zip(ωf, eachcol(F)))
        progress ? next!(p) : nothing

        D .= (K + 1im*ω*C  - ω^2*M)\Fₑ
        y[:, f] .= Sₒ*D

        if type == :vel
            y[:, f] .*= 1im*ω
        elseif type == :acc
            y[:, f] .*= -ω^2
        end
    end

    return FrequencySolution(y)
end