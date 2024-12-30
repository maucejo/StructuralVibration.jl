abstract type BarRod end

"""
Structure containing the data of a homogeneous and isotropic longitudinal bar

# Problem parameters
* L: Length [m]
* S: Cross-section area [m²]
* E: Young's modulus [Pa]
* ρ: Mass density [kg/m³]

# Model parameters
* L : Length [m]
* m : Line mass [kg/m]
* D : Stiffness coefficient [Pa]
"""
@with_kw struct Bar <: BarRod
    L::Float64
    m::Float64
    D::Float64

    function Bar(L::Float64, S::Float64, E::Float64, ρ::Float64)
        m = ρ*S
        D = E*S

        return new(L, m, D)
    end
end

"""
Structure containing the data of a homogeneous and isotropic torsional bar


# Problem parameters
* L: Length [m]
* I: Second-moment of area [m⁴]
* J: Torsion constant [m⁴]
* G: Shear modulus [Pa]
* ρ: Mass density [kg/m³]

# Model parameters
* L : Length [m]
* m : Line mass [kg/m]
* D : Stiffness coefficient [Pa]
"""
@with_kw struct Rod <: BarRod
    L :: Float64
    m::Float64
    D::Float64

    function Rod(L::Float64, I::Float64, J::Float64, G::Float64, ρ::Float64)
        m = ρ*I
        D = G*J

        return new(L, m, D)
    end
end

"""
    eigval(b::Bar, fmax, bc)
    eigval(b::Rod, fmax, bc)

Computes the natural frequencies of a longitudinal or torsional bar up to fmax

# Parameters
* p: Structure containing the bar data
* fmax: Maximum frequency for calculating the mode shapes [Hz]
* bc: Boundary conditions
    * :CC : Clamped - Clamped
    * :CF : Clamped - Free
    * :FF : Free - Free

# Outputs
* ωₙ: Natural frequencies calculated up to ωmax = 2π*fmax [Hz]
* kₙ: Vector of modal wavenumbers
"""
function eigval(b::BarRod, fmax, bc = :CC)
    (; L, m, D) = b

    c = sqrt(D/m)
    ωmax = 2π*fmax

    ωₙ = Float64[]
    kₙ = Float64[]
    if bc == :CC
        n = 1
        kᵢ = n*π/L
        ωᵢ = c*kᵢ
        while ωᵢ ≤ ωmax
            push!(ωₙ, ωᵢ)
            push!(kₙ, kᵢ)
            n += 1
            kᵢ = n*π/L
            ωᵢ = c*kᵢ
        end
    elseif bc === :CF
        n = 1
        kᵢ = (2n - 1)π/L
        ωᵢ = c*kᵢ
        while ωᵢ ≤ ωmax
            push!(ωₙ, ωᵢ)
            push!(kₙ, kᵢ)
            n += 1
            kᵢ = (2n - 1)π/L
            ωᵢ = c*kᵢ
        end
    elseif bc == :FF
        n = 0
        kᵢ = 0.
        ωᵢ = 0.
        while ωᵢ ≤ ωmax
            push!(ωₙ, ωᵢ)
            push!(kₙ, kᵢ)
            n += 1
            kᵢ = n*π/L
            ωᵢ = c*kᵢ
        end
    else
        error("Boundary conditions not implemented")
    end

    return ωₙ, kₙ
end

"""
    eigmode(b::Bar, kₙ, x, bc)
    eigmode(b::Rod, kₙ, x, bc)

Computes the mass-normalized mode shapes of a longitudinal or torsional bar

# Parameters
* b: Structure containing the bar data
* kₙ: Array of modal wavenumbers
* x: Coordinates of calculation points of the mode shapes
* bc: Boundary conditions
    * :CC : Clamped - Clamped
    * :CF : Clamped - Free
    * :FF : Free - Free

# Output
* ϕ: Mass-normalized mode shapes
"""
function eigmode(b::BarRod, kₙ, x, bc = :CC)
    (; L, m) = b

    if !isa(x, Array)
        x = collect(x);
    end

    if bc == :CC
        Mₙ = m*L/2

        return sin.(x*kₙ')./sqrt(Mₙ)
    elseif bc == :CF
        Mₙ = m*L/2

        return sin.(x*kₙ')./sqrt(Mₙ)
    elseif bc == :FF
        Mₙ = m*L.*ones(length(n))./2
        Mₙ[1] *= 2.

        return cos.(x*kₙ')./sqrt.(Mₙ')
    else
        error("Boundary conditions not implemented")
    end
end