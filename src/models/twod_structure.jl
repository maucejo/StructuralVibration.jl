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
"""
function modefreq(model::TwoDStructure, fmax)
   (; L, b, m, D) = model

    c = sqrt(D/m)
    ωmax = 2π*fmax

    m = 1
    n = 1
    km = m*π/L
    kn = n*π/b
    ωi = c*(km^2 + kn^2)

    ωmn = Float64[]
    kmn = Float64[]
    ind = Int64[]
    # Boucle sur m
    while ωi ≤ ωmax
        # Boucle sur n
        while ωi ≤ ωmax
            push!(ωmn,  ωi)
            append!(kmn, [km, kn])
            append!(ind, [m, n])
            n += 1
            km = m*π/L
            kn = n*π/b
            ωi = c*(km^2 + kn^2)
        end

        m += 1
        n = 1
        km = m*π/L
        kn = n*π/b
        ωi = c*(km^2 + kn^2)
    end

    kmn = reshape(kmn, (2, Int(length(kmn)/2)))
    ind = reshape(ind, (2, Int(length(ind)/2)))
    pos = sortperm(ωmn)

    return ωmn[pos], kmn[:, pos], ind[:, pos]
end

"""
    modeshape(model::Plate, kmn, x, y)
    modeshape(model::Membrane, kmn, x, y)

Computes the mass-normalized mode shapes of a simply supported rectangular plate or a clamped rectangular membrane

**Inputs**
* `model`: Structure containing the data related to the structure
* `kmn`: Matrix of modal wave numbers
* `(x, y)`: Coordinates of the points where the mode shapes are calculated

**Output**
* `ϕ`: Mass-normalized mode shapes
"""
function modeshape(p::TwoDStructure, kmn, x, y)
    (; L, b, m) = p

    if isa(x, Number)
        x = [x]
    else !isa(x, Array)
        x = collect(x)
    end

    if isa(y, Number)
        y = [y]
    else !isa(y, Array)
        y = collect(y)
    end

    Mn = m*L*b/4

    return sin.(x*kmn[1, :]').*sin.(y*kmn[2, :]')./sqrt(Mn)
end