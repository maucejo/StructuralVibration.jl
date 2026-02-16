abstract type TwoDStructure end

"""
    Plate(L, b, h, E, ρ, ν)

Structure containing the data of a homogeneous and isotropic bending plate

**Constructor parameters**
* `L`: Length [m]
* `b`: Width [m]
* `E`: Young's modulus [Pa]
* `ρ`: Density [kg/m³]
* `ν`: Poisson's ratio

**Fields**
* `L`: Length [m]
* `b`: Width [m]
* `m`: Surface mass [kg/m²]
* `D`: Bending stiffness [N.m]
"""
@show_data struct Plate{T <: Real} <: TwoDStructure
    L::T
    b::T
    m::T
    D::T

    function Plate(L::T, b::T, h::T, E::T, ρ::T, ν::T) where T
        m = ρ*h
        D = E*h^3/12/(1. - ν^2.)
        return new{T}(L, b, m, D)
    end
end

"""
    Membrane(L, b, m, D)

Structure containing the data of a homogeneous and isotropic rectangular membrane

**Constructor parameters**
* `L`: Length [m]
* `b`: Width [m]
* `m`: Surface mass [kg/m²]
* `T`: Tension per unit length [N/m]

**Fields**
* `L`: Length [m]
* `b`: Width [m]
* `m`: Surface mass [kg/m²]
* `D`: Tension per unit length [N/m]
"""
@show_data struct Membrane <: TwoDStructure
    L::Float64
    b::Float64
    m::Float64
    D::Float64

    function Membrane(L::Tm, b::Tm, m::Tm, T::Tm) where Tm
        return new(L, b, m, T)
    end
end

"""
    wave_parameters(model::Plate, freq)
    wave_parameters(model::Membrane, freq)

Computes the wave parameters of given model at a given frequency

**Inputs**
* `model`: Structure containing the bar data
* `freq`: Frequency at which the wave parameters are calculated [Hz]

**Outputs**
* `ω`: Angular frequency [rad/s]
* `c`: Wave speed [m/s]
* `k`: Wavenumber [rad/m]
* `λ`: Wavelength [m]
"""
function wave_parameters(model::Plate, freq)
    # Model parameters
    (; m, D) = model

    # Angular frequency
    ω = 2π*freq

    # Wavenumber
    k = (ω^2*m/D)^(1/4)

    # Wave speed
    c = ω/k

    # Wavelength
    λ = 2π/k

    return ω, c, k, λ
end

function wave_parameters(model::Membrane, freq)
    # Model parameters
    (; m, D) = model

    # Angular frequency
    ω = 2π*freq

    # Wave speed
    c = sqrt(D/m)

    # Wavenumber
    k = ω/c

    # Wavelength
    λ = 2π/k

    return ω, c, k, λ
end

"""
    modefreq(model::Plate, fmax)
    modefreq(model::Membrane, fmax)

Computes the natural frequencies of a simply supported rectangular plate or a clamped rectangular membrane up to fmax

**Inputs**
* `model`: Structure containing the data related to the plate
* `fmax`: Maximum frequency for calculating the modal shapes [Hz]

**Outputs**
* `ωmn`: Natural frequencies calculated up to ωmax = 2π*fmax [Hz]
* `kmn`: Matrix of modal wave numbers
* `ind`: Matrix of modal indices (m, n)
"""
function modefreq(model::TwoDStructure, fmax)
   (; L, b, m, D) = model

    c = sqrt(D/m)
    ωmax = 2π*fmax

    p = 1
    q = 1
    kp = p*π/L
    kq = q*π/b
    ωi = c*(kp^2 + kq^2)

    ωpq = typeof(D)[]
    kpq = typeof(D)[]
    ind = typeof(p)[]
    # Boucle sur m
    while ωi ≤ ωmax
        # Boucle sur n
        while ωi ≤ ωmax
            push!(ωpq,  ωi)
            append!(kpq, [kp, kq])
            append!(ind, [p, q])
            q += 1
            kp = p*π/L
            kq = q*π/b
            ωi = c*(kp^2 + kq^2)
        end

        p += 1
        q = 1
        kp = p*π/L
        kq = q*π/b
        ωi = c*(kp^2 + kq^2)
    end

    kpq = reshape(kpq, (2, Int(length(kpq)/2)))
    ind = reshape(ind, (2, Int(length(ind)/2)))
    pos = sortperm(ωpq)

    return ωpq[pos], kpq[:, pos], ind[:, pos]
end

"""
    modeshape(model::Plate, kpq, x, y)
    modeshape(model::Membrane, kpq, x, y)

Computes the mass-normalized mode shapes of a simply supported rectangular plate or a clamped rectangular membrane

**Inputs**
* `model`: Structure containing the data related to the structure
* `kmn`: Matrix of modal wave numbers
* `(x, y)`: Coordinates of the points where the mode shapes are calculated

**Output**
* `ϕ`: Mass-normalized mode shapes
"""
function modeshape(p::TwoDStructure, kpq, x, y)
    (; L, b, m) = p

    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end

    Mn = m*L*b/4

    ϕ = similar(kpq, length(x), size(kpq, 2))

    for (i, (xi, yi)) in enumerate(zip(x, y))
        @. ϕ[i, :] = sin(kpq[1, :]*xi).*sin(kpq[2, :]*yi)/sqrt(Mn)
    end

    return ϕ
end