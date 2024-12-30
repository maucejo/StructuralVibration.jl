"""
Structure containing the data feeding the modal solver for calculating an FRF

# Parameters
* ωₙ : Resonance frequencies
* ξₙ : Modal damping ratios
* ϕₑ : Mode shapes at excitation points
* ϕₒ : Mode shapes at observation points
* freq : Frequencies of interest

# Note
The mode shapes must be mass-normalized
"""
struct ModalFRF
    ωₙ :: Vector{Float64}
    ξₙ
    ϕₑ
    ϕₒ
    freq
end

"""
Structure containing the data feeding the direct solver for calculating an FRF

# Parameters
* K : Stiffness matrix
* M : Mass matrix
* C : Damping matrix
* exc_dofs : Degrees of freedom of excitation
* obs_dofs : Degrees of freedom of observation
* freq : Frequencies of interest
"""
struct DirectFRF
    K :: Matrix{Float64}
    M :: Matrix{Float64}
    C :: Matrix{Float64}
    exc_dofs
    obs_dofs
    freq
end

"""
    frf(m::ModalFRF)
    frf(m::DirectFRF)

Computes the FRF matrix by modal or direct approach

# Parameter
* m : Structure containing the problem data

# Output
* FRF : FRF matrix
"""
function frf(m::ModalFRF, type = :dis)
    # Initialisation
    (; ωₙ, ξₙ, ϕₑ, ϕₒ, freq) = m
    Nₑ = size(ϕₑ, 1)
    Nₒ = size(ϕₒ, 1)
    Nf = length(freq)

    FRF = [Matrix{ComplexF64}(undef, Nₒ, Nₑ) for _ in 1:Nf]

    ωf = 2π*freq
    p = Progress(Nf, color=:black, barlen=75, showspeed=true)
    @inbounds for (f, ω) in enumerate(ωf)
        next!(p)
        M = Diagonal(@. 1/(ωₙ^2 - ω^2 + 2im*ξₙ*ωₙ*ω))
        FRF[f] = ϕₑ*M*ϕₒ'

        if type == :vel
            FRF[f] *= 1im*ω
        elseif type == :acc
            FRF[f] *= -ω^2
        end
    end

    return FRF
end