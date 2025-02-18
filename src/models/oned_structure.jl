abstract type OneDStructure end
abstract type WaveEquation <: OneDStructure end

"""
    Bar(L, S, E, ρ)

Structure containing the data of a homogeneous and isotropic longitudinal bar

**Constructor parameters**
* `L`: Length [m]
* `S`: Cross-section area [m²]
* `E`: Young's modulus [Pa]
* `ρ`: Mass density [kg/m³]

**Fields**
* `L`: Length [m]
* `m`: Line mass [kg/m]
* `D`: Stiffness coefficient [Pa]
"""
@with_kw struct Bar <: WaveEquation
    L::Float64
    m::Float64
    D::Float64

    Bar(L::Float64, S::Float64, E::Float64, ρ::Float64) = new(L, ρ*S, E*S)
end

"""
    Rod(L, I, J, G, ρ)

Structure containing the data of a homogeneous and isotropic torsional bar

**Constructor parameters**
* L: Length [m]
* I: Second-moment of area [m⁴]
* J: Torsion constant [m⁴]
* G: Shear modulus [Pa]
* ρ: Mass density [kg/m³]

**Fields**
* `L`: Length [m]
* `m`: Line mass [kg/m]
* `D`: Stiffness coefficient [Pa]
"""
@with_kw struct Rod <: WaveEquation
    L::Float64
    m::Float64
    D::Float64

    Rod(L::Float64, I::Float64, J::Float64, G::Float64, ρ::Float64) = new(L, ρ*I, G*J)
end

"""
    Strings(L, S, D, ρ)

Structure containing the data of a homogeneous and isotropic string

**Constructor parameters**
* `L`: Length [m]
* `S`: Cross-section area [m²]
* `D`: Tension [N]
* `ρ`: Mass density [kg/m³]

**Fields**
* `L`: Length [m]
* `m`: Linear mass density [kg/m]
* `D`: Tension [N]
"""
@with_kw struct Strings <: WaveEquation
    L :: Float64
    m :: Float64
    D :: Float64

    Strings(L::Float64, S::Float64, D::Float64, ρ::Float64) = new(L, ρ*S, D)
end

"""
    Beam(L, S, I, E, ρ)

Structure containing the data of a homogeneous and isotropic bending beam

**Constructor parameters**
* `L`: Length [m]
* `S`: Cross-section area [m²]
* `I`: Second moment of area [m⁴]
* `E`: Young's modulus [Pa]
* `ρ`: Density [kg/m³]

**Fields**
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
    modefreq(model::Bar, fmax, bc = :CC)
    modefreq(model::Rod, fmax, bc = :CC)
    modefreq(model::Strings, fmax, bc = :CC)
    modefreq(model::Beam, fmax, bc = :SS)

Computes the natural frequencies of a longitudinal or torsional bar up to fmax

**Inputs**
* `model`: Structure containing the bar data
* `fmax`: Maximum frequency for calculating the mode shapes [Hz]
* `bc`: Boundary conditions
    * For all OneDStructure
        * `:CC`: Clamped - Clamped
        * `:CF`: Clamped - Free
        * `:FF`: Free - Free
    * For beams
        * `:SS`: Simply Supported - Simply Supported
        * `:SC`: Simply Supported - Clamped
        * `:SF`: Simply Supported - Free

**Outputs**
* `ωₙ`: Natural frequencies calculated up to ωmax = 2π*fmax [Hz]
* `kₙ`: Vector of modal wavenumbers
"""
function modefreq(model::WaveEquation, fmax, bc = :CC)
    (; L, m, D) = model

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

function modefreq(model::Beam, fmax, bc = :SS)
    (; L, m, D) = model

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
    modeshape(model::Bar, kₙ, x, bc = :CC)
    modeshape(model::Rod, kₙ, x, bc = :CC)
    modeshape(model::Strings, kₙ, x, bc = :CC)
    modeshape(model::Beam, kₙ, x, bc = :SS)

Computes the mass-normalized mode shapes of a longitudinal or torsional bar

**Inputs**
* `model`: Structure containing the bar data
* `kₙ`: Array of modal wavenumbers
* `x`: Coordinates of calculation points of the mode shapes
* `bc`: Boundary conditions
    * For all OneDStructure
        * `:CC`: Clamped - Clamped
        * `:CF`: Clamped - Free
        * `:FF`: Free - Free
    * For beams
        * `:SS`: Simply Supported - Simply Supported
        * `:SC`: Simply Supported - Clamped
        * `:SF`: Simply Supported - Free

**Output**
* `ϕ`: Mass-normalized mode shapes
"""
function modeshape(model::WaveEquation, kₙ, x, bc = :CC)

    (; L, m) = model

    if x isa Number
        x = [x]
    else !(x isa Array)
        x = collect(x)
    end

    if bc == :CC
        # Modal mass
        M = m*L/2

        return sin.(x*kₙ')./sqrt(M)
    elseif bc == :CF
        # Modal mass
        M = m*L/2

        return sin.(x*kₙ')./sqrt(M)
    elseif bc == :FF
        # Modal mass
        n = length(kₙ)
        Mₙ = m*L.*ones(length(n))./2
        Mₙ[1] *= 2.

        return cos.(x*kₙ')./sqrt.(Mₙ')
    else
        error("Boundary conditions not implemented")
    end
end

function modeshape(model::Beam, kₙ, x, bc = :SS)

    (; L, m) = model

    if x isa Number
        x = [x]
    else !(x isa Array)
        x = collect(x)
    end

    if bc == :SS
        # Modal mass
        M = m*L/2.

        return sin.(x*kₙ')./sqrt(M)
    elseif bc == :CC
        Mₙ = @. m*(
                -kₙ*L*cos(2kₙ*L)
                + kₙ*L*cosh(2kₙ*L)
                + cosh(kₙ*L)^2*sin(2kₙ*L)
                + 2cos(kₙ*L)*sinh(kₙ*L)
                - 2sin(kₙ*L)*(cosh(kₙ*L) + 2kₙ*L*sinh(kₙ*L))
                - cos(kₙ*L)^2*sinh(2kₙ*L)
                )/(2kₙ*(cos(kₙ*L) - cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L))^2)

        return @. ((cosh(x*kₙ') - cos(x*kₙ'))/(cosh(kₙ'*L) - cos(kₙ'*L)) - (sinh(x*kₙ') - sin(x*kₙ'))/(sinh(kₙ'*L) - sin(kₙ'*L)))/sqrt(Mₙ')
    elseif bc == :CS
        Mₙ = @. m*(
                -kₙ*L*cos(2kₙ*L)
                + kₙ*L*cosh(2kₙ*L)
                + cosh(kₙ*L)^2*sin(2kₙ*L)
                + 2cos(kₙ*L)*sinh(kₙ*L)
                - 2sin(kₙ*L)*(cosh(kₙ*L) + 2kₙ*L*sinh(kₙ*L))
                - cos(kₙ*L)^2*sinh(2kₙ*L)
                )/(2kₙ*(cos(kₙ*L) - cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L))^2)

        return @. ((cosh(x*kₙ') - cos(x*kₙ'))/(cosh(kₙ'*L) - cos(kₙ'*L)) - (sinh(x*kₙ') - sin(x*kₙ'))/(sinh(kₙ'*L) - sin(kₙ'*L)))/sqrt(Mₙ')

    elseif bc == :CF
        Mₙ = @. m*(
                -2kₙ*L*cos(2kₙ*L)
                - 7sin(2kₙ*L)
                + cosh(2kₙ*L)*(2kₙ*L + 3sin(2kₙ*L))
                - 2cosh(kₙ*L)*(3sin(kₙ*L) + sin(3kₙ*L))
                - sin(4kₙ*L)
                + (6cos(kₙ*L) + 6cos(3kₙ*L) - 8kₙ*L*sin(kₙ*L))*sinh(kₙ*L)
                + 6cos(kₙ*L)^2*sinh(2kₙ*L)
                )/(4kₙ*(cos(kₙ*L) + cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L)).^2)

        return @. ((cosh(x*kₙ') - cos(x*kₙ'))/(cosh(kₙ'*L) + cos(kₙ'*L)) - (sinh(x*kₙ') - sin(x*kₙ'))/(sinh(kₙ'*L) + sin(kₙ'*L)))/sqrt(Mₙ')

    elseif bc == :SF
        Mₙ = @. m*(
            -3/tan(kₙ*L)
            + 3/tanh(kₙ*L)
            + kₙ*L/sin(kₙ*L)^2
            - kₙ*L/sinh(kₙ*L)^2)/2kₙ

        return @. (sin(x*kₙ')/sin(kₙ'*L) + sinh(x*kₙ')/sinh(kₙ'*L))/sqrt(Mₙ')

    elseif bc == :FF
        Mₙ = @. m*(
                -kₙ*L*cos(2kₙ*L)
                + kₙ*L*cosh(2kₙ*L)
                + 6cosh(kₙ*L)*sin(kₙ*L)
                - 3cosh(kₙ*L)^2*sin(2kₙ*L)
                - (6cos(kₙ*L) + 4kₙ*L*sin(kₙ*L))*sinh(kₙ*L)
                + 3cos(kₙ*L)^2*sinh(2kₙ*L)
            )/(2kₙ*(cos(kₙ*L) - cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L))^2)

        Mₙ[1:2] .= m*L

        ϕₙ = @. ((cosh(x*kₙ') + cos(x*kₙ'))/(cosh(kₙ'*L) - cos(kₙ'*L))) - ((sinh(x*kₙ') + sin(x*kₙ'))/(sinh(kₙ'*L) - sin(kₙ'*L)))

        ϕₙ[:, 1] .= 1.
        ϕₙ[:, 2] = x .- L/2

        return @. ϕₙ/sqrt(Mₙ')
    else
        error("Boundary conditions not implemented")
    end
end