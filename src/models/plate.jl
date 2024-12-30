"""
Structure containing the data of a homogeneous and isotropic bending plate

# Problem parameters
* L : Length [m]
* b : Width [m]
* E : Young's modulus [Pa]
* ρ : Density [kg/m³]
* ν : Poisson's ratio

# Model parameters
* L : Length [m]
* b : Width [m]
* m : Surface mass [kg/m²]
* D : Bending stiffness [N.m]
"""
@with_kw struct Plate
    L::Float64
    b::Float64
    m::Float64
    D::Float64

    function Plate(L::Float64, b::Float64, h::Float64, E::Float64, ρ::Float64, ν::Float64)

        m = ρ*h
        D = E*h^3/12/(1. - ν^2.)

        return new(L, b, m, D)
    end
end

"""
    eigval(p::Plate, fmax)

Computes the natural frequencies of a simply supported plate up to fmax

# Parameters
    * p : Structure containing the data related to the plate
    * fmax : Maximum frequency for calculating the modal shapes [Hz]

# Outputs
    * ωₘₙ : Natural frequencies calculated up to ωmax = 2π*fmax [Hz]
    * kₘₙ : Matrix of modal wave numbers
"""
function eigval(p::Plate, fmax)
   (; L, b, m, D) = p

    c = sqrt(D/m)
    ωₘₐₓ = 2π*fmax

    m = 1
    n = 1
    kₘ = m*π/L
    kₙ = n*π/b
    ωᵢ = c*(kₘ^2 + kₙ^2)

    ωₘₙ = Float64[]
    kₘₙ = Float64[]
    ind = Int64[]
    # Boucle sur m
    while ωᵢ ≤ ωₘₐₓ
        # Boucle sur n
        while ωᵢ ≤ ωₘₐₓ
            push!(ωₘₙ,  ωᵢ)
            append!(kₘₙ, [kₘ, kₙ])
            append!(ind, [m, n])
            n += 1
            kₘ = m*π/L
            kₙ = n*π/b
            ωᵢ = c*(kₘ^2 + kₙ^2)
        end

        m += 1
        n = 1
        kₘ = m*π/L
        kₙ = n*π/b
        ωᵢ = c*(kₘ^2 + kₙ^2)
    end

    kₘₙ = reshape(kₘₙ, (2, Int(length(kₘₙ)/2)))
    ind = reshape(ind, (2, Int(length(ind)/2)))
    pos = sortperm(ωₘₙ)

    return ωₘₙ[pos], kₘₙ[:, pos], ind[:, pos]
end

"""
    eigmode(p::Plate, kₘₙ, X, Y)

Computes the mass-normalized mode shapes of a simply supported plate

# Parameters
    * p : Structure containing the data related to the plate
    * kₘₙ : Matrix of modal wave numbers
    * (X, Y): Coordinates of the points where the mode shapes are calculated

# Output
    * ϕ: Mass-normalized mode shapes
"""
@views function eigmode(p::Plate, kₘₙ, X, Y)
    (; L, b, m) = p
    Mₙ = m*L*b/4

    return sin.(X*kₘₙ[1, :]').*sin.(Y*kₘₙ[2, :]')./sqrt(Mₙ)
end