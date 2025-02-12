abstract type OneDStructure end
abstract type BarRodString <: OneDStructure end

"""
    Bar(L, S, E, ρ)

Structure containing the data of a homogeneous and isotropic longitudinal bar

# Constructor parameters
* `L`: Length [m]
* `S`: Cross-section area [m²]
* `E`: Young's modulus [Pa]
* `ρ`: Mass density [kg/m³]

# Fields
* `L`: Length [m]
* `m`: Line mass [kg/m]
* `D`: Stiffness coefficient [Pa]
"""
@with_kw struct Bar <: BarRodString
    L::Float64
    m::Float64
    D::Float64

    Bar(L::Float64, S::Float64, E::Float64, ρ::Float64) = new(L, ρ*S, E*S)
end

"""
    Rod(L, I, J, G, ρ)

Structure containing the data of a homogeneous and isotropic torsional bar

# Constructor parameters
* L: Length [m]
* I: Second-moment of area [m⁴]
* J: Torsion constant [m⁴]
* G: Shear modulus [Pa]
* ρ: Mass density [kg/m³]

# Fields
* `L`: Length [m]
* `m`: Line mass [kg/m]
* `D`: Stiffness coefficient [Pa]
"""
@with_kw struct Rod <: BarRodString
    L::Float64
    m::Float64
    D::Float64

    Rod(L::Float64, I::Float64, J::Float64, G::Float64, ρ::Float64) = new(L, ρ*I, G*J)
end

"""
    Strings(L, m, D)

Structure containing the data of a homogeneous and isotropic string

# Fields
* `L`: Length [m]
* `m`: Linear mass density [kg/m]
* `D`: Tension [N]
"""
@with_kw struct Strings <: BarRodString
    L :: Float64
    m :: Float64
    D :: Float64
end

"""
    Beam(L, S, I, E, ρ)

Structure containing the data of a homogeneous and isotropic bending beam

# Constructor parameters
* `L`: Length [m]
* `S`: Cross-section area [m²]
* `I`: Second moment of area [m⁴]
* `E`: Young's modulus [Pa]
* `ρ`: Density [kg/m³]

# Fields
* `L`: Length [m]
* `M`: Linear mass density [kg/m]
* `D`: Bending stiffness [N.m²]
"""
@with_kw struct Beam <: OneDStructure
    L :: Float64
    m :: Float64
    D :: Float64

    Beam(L::Float64, S::Float64, I::Float64, E::Float64, ρ::Float64) = new(L, ρ*S, E*I)
end

"""
    modefreq(b::Bar, fmax, bc)
    modefreq(r::Rod, fmax, bc)
    modefreq(s::Strings, fmax, bc)

Computes the natural frequencies of a longitudinal or torsional bar up to fmax

# Inputs
* `p`: Structure containing the bar data
* `fmax`: Maximum frequency for calculating the mode shapes [Hz]
* `bc`: Boundary conditions
    * :CC: Clamped - Clamped
    * :CF: Clamped - Free
    * :FF: Free - Free

# Outputs
* `ωₙ`: Natural frequencies calculated up to ωmax = 2π*fmax [Hz]
* `kₙ`: Vector of modal wavenumbers
"""
function modefreq(b::BarRodString, fmax, bc = :CC)
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
    modefreq(b::Beam, fmax, bc)

Computes the natural frequencies of a beam in bending up to fmax

# Inputs
* `p`: Structure containing the data related to the beam
* `fmax`: Maximum frequency for calculating the modal shapes [Hz]
* `bc`: Boundary conditions
    * :SS: Simply Supported - Simply Supported
    * :CC: Clamped - Clamped
    * :SC: Simply Supported - Clamped
    * :CF: Clamped - Free
    * :SF: Simply Supported - Free
    * :FF: Free - Free

# Outputs
* `ωₙ`: Natural frequencies calculated up to ωmax = 2π*fmax [Hz]
* `kₙ`: Vector of modal wave numbers
"""
function modefreq(b::Beam, fmax, bc = :SS)
    (; L, m, D) = b

    c = sqrt(D/m)
    ωmax = 2π*fmax

    ωₙ = Float64[]
    kₙ = Float64[]
    if bc == :SS
        n = 1
        kᵢ = n*π/L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωmax
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = n*π/L
            ωᵢ = c*kᵢ^2
        end
    elseif bc == :CC
        append!(kₙ, [4.73, 7.85, 11]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 4
        kᵢ = (2n + 1)π/2L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωmax
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (2n + 1)π/2L
            ωᵢ = c*kᵢ.^2
        end
    elseif bc == :CS
        append!(kₙ, [3.92, 7.07, 10.2]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 4
        kᵢ = (4n + 1)π/4L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωmax
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (4n + 1)π/4L
            ωᵢ = c*kᵢ^2
        end
    elseif bc == :CF
        append!(kₙ, [1.87,  4.73, 7.85, 11]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 5
        kᵢ = (2n + 1)π/2L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωmax
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (2n + 1)π/2L
            ωᵢ = c*kᵢ^2
        end
    elseif bc == :SF
        append!(kₙ, [3.92, 7.07, 10.2]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 4
        kᵢ = (4n + 1)π/4L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωmax
            push!(kₙ,kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (4n + 1)π/4L
            ωᵢ = c*kᵢ^2
        end
    elseif bc == :FF
        append!(kₙ, [0., 0., 4.73, 7.85, 11]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 4
        kᵢ = (2n + 1)π/2L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωmax
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (2n + 1)π/2L
            ωᵢ = c*kᵢ^2
        end
    else
        error("Boundary conditions not implemented")
    end

    return ωₙ, kₙ
end

"""
    modeshape(b::Bar, kₙ, x, bc)
    modeshape(b::Rod, kₙ, x, bc)

Computes the mass-normalized mode shapes of a longitudinal or torsional bar

# Inputs
* `b`: Structure containing the bar data
* `kₙ`: Array of modal wavenumbers
* `x`: Coordinates of calculation points of the mode shapes
* `bc`: Boundary conditions
    * :CC: Clamped - Clamped
    * :CF: Clamped - Free
    * :FF: Free - Free

# Output
* `ϕ`: Mass-normalized mode shapes
"""
function modeshape(b::BarRodString, kₙ, x, bc = :CC)
    (; L, m) = b

    if isa(x, Number)
        x = [x]
    else !isa(x, Array)
        x = collect(x)
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

"""
    modeshape(b::Beam, kₙ, x, bc)

Calculates the mass-normalized mode shapes of a beam in bending

# Inputs
* `b`: Structure containing the data related to the beam
* `kₙ`: Vector of modal wave numbers
* `x`: Coordinates of the points where the mode shapes are calculated
* `bc`: Boundary conditions
    * :SS: Simply Supported - Simply Supported
    * :CC: Clamped - Clamped
    * :CS: Clamped - Simply Supported
    * :CF: Clamped - Free
    * :SF: Simply Supported - Free
    * :FF: Free - Free

# Output
* `ϕ`: Mass-normalized mode shapes
"""
function modeshape(b::Beam, kₙ, x, bc = :SS)
    (; L, m) = b

    if isa(x, Number)
        x = [x]
    else !isa(x, Array)
        x = collect(x)
    end

    if bc == :SS
        Mₙ = m*L/2.
        ϕₙ = sin.(x*kₙ')
    elseif bc == :CC
        Mₙ = @. m*((-kₙ*L*cos(2kₙ*L) + kₙ*L*cosh(2kₙ*L) + cosh(kₙ*L)^2*sin(2kₙ*L) + 2cos(kₙ*L)*sinh(kₙ*L) - 2sin(kₙ*L)*(cosh(kₙ*L) + 2kₙ*L*sinh(kₙ*L)) - cos(kₙ*L)^2*sinh(2kₙ*L))/(2kₙ*(cos(kₙ*L) - cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L))^2))

        ϕₙ = @. ((cosh(x*kₙ') - cos(x*kₙ'))/(cosh(kₙ'*L) - cos(kₙ'*L))) - ((sinh(x*kₙ') - sin(x*kₙ'))/(sinh(kₙ'*L) - sin(kₙ'*L)))
    elseif bc == :CS
        Mₙ = @. m*((-kₙ*L*cos(2kₙ*L) + kₙ*L*cosh(2kₙ*L) + cosh(kₙ*L)^2*sin(2kₙ*L) + 2cos(kₙ*L)*sinh(kₙ*L) - 2sin(kₙ*L)*(cosh(kₙ*L) + 2kₙ*L*sinh(kₙ*L)) - cos(kₙ*L)^2*sinh(2kₙ*L))/(2kₙ*(cos(kₙ*L) - cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L))^2))

        ϕₙ = @. ((cosh(x*kₙ') - cos(x*kₙ'))/(cosh(kₙ'*L) - cos(kₙ'*L))) - ((sinh(x*kₙ') - sin(x*kₙ'))/(sinh(kₙ'*L) - sin(kₙ'*L)))

    elseif bc == :CF
        Mₙ = @. m*((-2kₙ*L*cos(2kₙ*L) - 7sin(2kₙ*L) + cosh(2kₙ*L)*(2kₙ*L + 3sin(2kₙ*L)) - 2cosh(kₙ*L)*(3sin(kₙL) + sin(3kₙ*L)) - sin(4kₙ*L) + 2(3(cos(kₙ*L) + cos(3kₙ*L)) - 4kₙ*L*sin(kₙ*L))*sinh(kₙ*L) + 6cos(kₙ*L)^2*sinh(2kₙ*L))/(4kₙ*(cos(kₙ*L) + cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L)).^2))

        ϕₙ = @. ((cosh(x*kₙ') - cos(x*kₙ'))/(cosh(kₙ'*L) + cos(kₙ'*L))) - ((sinh(x*kₙ') - sin(x*kₙ'))/(sinh(kₙ'*L) + sin(kₙ'*L)))
    elseif bc == :SF
        Mₙ = @. m*((-3/tan(kₙ*L) + 3/tanh(kₙ*L) + kₙ*L/(sin(kₙ*L)^2) - kₙ*L/(sinh(kₙ*L)^2))/2kₙ);

        ϕₙ = @. (sin(x*kₙ')/sin(kₙ'*L)) + (sinh(x*kₙ')/sinh(kₙ'*L))
    elseif bc == :FF
        Mₙ = @. m*((-kₙ*L*cos(2kₙ*L) + kₙ*L*cosh(2kₙ*L) + 6cosh(kₙ*L)*sin(kₙ*L) - 3cosh(kₙ*L)^2*sin(2kₙ*L) - 2(3cos(kₙ*L) + 2kₙ*L*sin(kₙ*L))*sinh(kₙ*L) + 3cos(kₙ*L)^2*sinh(2kₙ*L))/(2kₙ*(cos(kₙ*L) - cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L))^2))
        Mₙ[1:2] .= m*L

        ϕₙ = @. ((cosh(x*kₙ') + cos(x*kₙ'))/(cosh(kₙ'*L) - cos(kₙ'*L))) - ((sinh(x*kₙ') + sin(x*kₙ'))/(sinh(kₙ'*L) - sin(kₙ'*L)))

        ϕₙ[:, 1] .= 1.
        ϕₙ[:, 2] = x .- L/2
    else
        error("Boundary conditions not implemented")
    end

    return @. ϕₙ/sqrt(Mₙ')
end