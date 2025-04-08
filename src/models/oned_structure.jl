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
@show_data struct Bar{T <: Real} <: WaveEquation
    L::T
    m::T
    D::T

    Bar(L::T, S::T, E::T, ρ::T) where T = new{T}(L, ρ*S, E*S)
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
@show_data struct Rod{T <:Real} <: WaveEquation
    L::T
    m::T
    D::T

    Rod(L::T, I::T, J::T, G::T, ρ::T) where T = new{T}(L, ρ*I, G*J)
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
@show_data struct Strings{T <: Real} <: WaveEquation
    L::T
    m::T
    D::T

    Strings(L::T, S::T, D::T, ρ::T) where T = new{T}(L, ρ*S, D)
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
@show_data struct Beam{T <: Real} <: OneDStructure
    L::T
    m::T
    D::T

    Beam(L::T, S::T, I::T, E::T, ρ::T) where T = new{T}(L, ρ*S, E*I)
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
* `ωn`: Natural frequencies calculated up to ωmax = 2π*fmax [Hz]
* `kn`: Vector of modal wavenumbers
"""
function modefreq(model::WaveEquation, fmax, bc = :CC)
    (; L, m, D) = model

    c = sqrt(D/m)
    ωmax = 2π*fmax

    ωn = Float64[]
    kn = Float64[]
    if bc == :CC
        n = 1
        ki = n*π/L
        ωi = c*ki
        while ωi ≤ ωmax
            push!(ωn, ωi)
            push!(kn, ki)
            n += 1
            ki = n*π/L
            ωi = c*ki
        end
    elseif bc === :CF
        n = 1
        ki = (2n - 1)π/L
        ωi = c*ki
        while ωi ≤ ωmax
            push!(ωn, ωi)
            push!(kn, ki)
            n += 1
            ki = (2n - 1)π/L
            ωi = c*ki
        end
    elseif bc == :FF
        n = 0
        ki = 0.
        ωi = 0.
        while ωi ≤ ωmax
            push!(ωn, ωi)
            push!(kn, ki)
            n += 1
            ki = n*π/L
            ωi = c*ki
        end
    else
        error("Boundary conditions not implemented")
    end

    return ωn, kn
end

function modefreq(model::Beam, fmax, bc = :SS)
    (; L, m, D) = model

    c = sqrt(D/m)
    ωmax = 2π*fmax

    ωn = Float64[]
    kn = Float64[]
    if bc == :SS
        n = 1
        ki = n*π/L
        ωi = c*ki^2
        while ωi ≤ ωmax
            push!(kn, ki)
            push!(ωn, ωi)
            n += 1
            ki = n*π/L
            ωi = c*ki^2
        end
    elseif bc == :CC
        append!(kn, [4.73, 7.85, 11]./L)
        append!(ωn, c.*kn.^2)

        n = 4
        ki = (2n + 1)π/2L
        ωi = c*ki^2
        while ωi ≤ ωmax
            push!(kn, ki)
            push!(ωn, ωi)
            n += 1
            ki = (2n + 1)π/2L
            ωi = c*ki.^2
        end
    elseif bc == :CS
        append!(kn, [3.92, 7.07, 10.2]./L)
        append!(ωn, c.*kn.^2)

        n = 4
        ki = (4n + 1)π/4L
        ωi = c*ki^2
        while ωi ≤ ωmax
            push!(kn, ki)
            push!(ωn, ωi)
            n += 1
            ki = (4n + 1)π/4L
            ωi = c*ki^2
        end
    elseif bc == :CF
        append!(kn, [1.87,  4.73, 7.85, 11]./L)
        append!(ωn, c.*kn.^2)

        n = 5
        ki = (2n + 1)π/2L
        ωi = c*ki^2
        while ωi ≤ ωmax
            push!(kn, ki)
            push!(ωn, ωi)
            n += 1
            ki = (2n + 1)π/2L
            ωi = c*ki^2
        end
    elseif bc == :SF
        append!(kn, [3.92, 7.07, 10.2]./L)
        append!(ωn, c.*kn.^2)

        n = 4
        ki = (4n + 1)π/4L
        ωi = c*ki^2
        while ωi ≤ ωmax
            push!(kn,ki)
            push!(ωn, ωi)
            n += 1
            ki = (4n + 1)π/4L
            ωi = c*ki^2
        end
    elseif bc == :FF
        append!(kn, [0., 0., 4.73, 7.85, 11]./L)
        append!(ωn, c.*kn.^2)

        n = 4
        ki = (2n + 1)π/2L
        ωi = c*ki^2
        while ωi ≤ ωmax
            push!(kn, ki)
            push!(ωn, ωi)
            n += 1
            ki = (2n + 1)π/2L
            ωi = c*ki^2
        end
    else
        error("Boundary conditions not implemented")
    end

    return ωn, kn
end

"""
    modeshape(model::Bar, kn, x, bc = :CC)
    modeshape(model::Rod, kn, x, bc = :CC)
    modeshape(model::Strings, kn, x, bc = :CC)
    modeshape(model::Beam, kn, x, bc = :SS)

Computes the mass-normalized mode shapes of a longitudinal or torsional bar

**Inputs**
* `model`: Structure containing the bar data
* `kn`: Array of modal wavenumbers
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
function modeshape(model::WaveEquation, kn, x, bc = :CC)

    (; L, m) = model

    ϕ = similar(kn, length(x), length(kn))

    if bc == :CC || bc == :CF
        # Modal mass
        M = m*L/2

        for (i, xi) in enumerate(x)
            @. ϕ[i, :] = sin(xi*kn)/sqrt(M)
        end

    elseif bc == :FF
        # Modal mass
        n = length(kn)
        Mn = m*L.*ones(length(n))./2
        Mn[1] *= 2.

        for (i, xi) in enumerate(x)
            @. ϕ[i, :] = cos(xi*kn)/sqrt(Mn)
        end
    else
        throw(ArgumentError("Boundary conditions not implemented"))
    end

    return ϕ
end

function modeshape(model::Beam, kn, x, bc = :SS)

    (; L, m) = model

    ϕ = similar(kn, length(x), length(kn))

    if bc == :SS
        # Modal mass
        M = m*L/2.

        for (i, xi) in enumerate(x)
           @. ϕ[i, :] = sin(xi*kn)/sqrt(M)
        end
    elseif bc == :CC
        Mn = @. m*(
                -kn*L*cos(2kn*L)
                + kn*L*cosh(2kn*L)
                + cosh(kn*L)^2*sin(2kn*L)
                + 2cos(kn*L)*sinh(kn*L)
                - 2sin(kn*L)*(cosh(kn*L) + 2kn*L*sinh(kn*L))
                - cos(kn*L)^2*sinh(2kn*L)
                )/(2kn*(cos(kn*L) - cosh(kn*L))^2*(sin(kn*L) - sinh(kn*L))^2)

        for (i, xi) in enumerate(x)
            @. ϕ[i, :] = ((cosh(xi*kn) - cos(xi*kn))/(cosh(kn*L) - cos(kn*L)) - (sinh(xi*kn) - sin(xi*kn))/(sinh(kn*L) - sin(kn*L)))/sqrt(Mn)
        end

    elseif bc == :CS
        Mn = @. m*(
                -kn*L*cos(2kn*L)
                + kn*L*cosh(2kn*L)
                + cosh(kn*L)^2*sin(2kn*L)
                + 2cos(kn*L)*sinh(kn*L)
                - 2sin(kn*L)*(cosh(kn*L) + 2kn*L*sinh(kn*L))
                - cos(kn*L)^2*sinh(2kn*L)
                )/(2kn*(cos(kn*L) - cosh(kn*L))^2*(sin(kn*L) - sinh(kn*L))^2)

        for (i, xi) in enumerate(x)
            @. ϕ[i, :] = ((cosh(xi*kn) - cos(xi*kn))/(cosh(kn*L) - cos(kn*L)) - (sinh(xi*kn) - sin(xi*kn))/(sinh(kn*L) - sin(kn*L)))/sqrt(Mn)
        end

    elseif bc == :CF
        Mn = @. m*(
                -2kn*L*cos(2kn*L)
                - 7sin(2kn*L)
                + cosh(2kn*L)*(2kn*L + 3sin(2kn*L))
                - 2cosh(kn*L)*(3sin(kn*L) + sin(3kn*L))
                - sin(4kn*L)
                + (6cos(kn*L) + 6cos(3kn*L) - 8kn*L*sin(kn*L))*sinh(kn*L)
                + 6cos(kn*L)^2*sinh(2kn*L)
                )/(4kn*(cos(kn*L) + cosh(kn*L))^2*(sin(kn*L) - sinh(kn*L)).^2)

        for (i, xi) in enumerate(x)
            @. ϕ[i, :] = ((cosh(xi*kn) - cos(xi*kn))/(cosh(kn*L) + cos(kn*L)) - (sinh(xi*kn) - sin(xi*kn))/(sinh(kn*L) + sin(kn*L)))/sqrt(Mn)
        end

    elseif bc == :SF
        Mn = @. m*(
            -3/tan(kn*L)
            + 3/tanh(kn*L)
            + kn*L/sin(kn*L)^2
            - kn*L/sinh(kn*L)^2)/2kn

        for (i, xi) in enumerate(x)
            @. ϕ[i, :] = (sin(xi*kn)/sin(kn*L) + sinh(xi*kn)/sinh(kn*L))/sqrt(Mn)
        end

    elseif bc == :FF
        Mn = @. m*(
                -kn*L*cos(2kn*L)
                + kn*L*cosh(2kn*L)
                + 6cosh(kn*L)*sin(kn*L)
                - 3cosh(kn*L)^2*sin(2kn*L)
                - (6cos(kn*L) + 4kn*L*sin(kn*L))*sinh(kn*L)
                + 3cos(kn*L)^2*sinh(2kn*L)
            )/(2kn*(cos(kn*L) - cosh(kn*L))^2*(sin(kn*L) - sinh(kn*L))^2)

        Mn[1:2] .= m*L

        for (i, xi) in enumerate(x)
            @. ϕ[i, :] = ((cosh(xi*kn) + cos(xi*kn))/(cosh(kn*L) - cos(kn*L))) - ((sinh(xi*kn) + sin(xi*kn))/(sinh(kn*L) - sin(kn*L)))/sqrt(Mn)
        end

        ϕ[:, 1] .= 1/sqrt(Mn[1])
        ϕ[:, 2] = (x .- L/2)/sqrt(Mn[2])

    else
        throw(ArgumentError("Boundary conditions not implemented"))
    end

    return ϕ
end