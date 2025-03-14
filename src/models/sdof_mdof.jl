"""
    Sdof(m, ω0, ξ)

Structure containing the data of a sdof system

**Constructor**
* `m`: Mass [kg]
* `f0`: Natural frequency [Hz]
* `ξ`: Damping ratio

**Fields**
* `m`: Mass [kg]
* `ω0`: Natural frequency [rad/s]
* `ξ`: Damping ratio
"""
@with_kw struct Sdof
    m :: Float64
    ω0 ::Float64
    ξ :: Float64

    Sdof(m, f0, ξ) = new(m, 2π*f0, ξ)
end

"""
    Mdof(k, m, c = Float64[])

Structure containing the data for building a mdof system

**Fields**
* `k`: Stiffness coefficients of the spring elements
* `m`: Masses of the mdof system
* `c`: Damping coefficients of the viscous dampers
"""
@with_kw struct Mdof
    k::Vector{Float64}
    m::Vector{Float64}
    c::Vector{Float64}

    Mdof(k, m, c = Float64[]) = new(k, m, c)
end

"""
    MdofMesh(Elt, constrained_dofs, free_dofs)

Structure containing the data for building a mdof mesh

**Constructor**
* `model`: Mdof model
* `bc`: Boundary conditions
    * `:CC`: Clamped - Clamped
    * `:CF`: Clamped - Free
    * `:FF`: Free - Free

**Fields**
* `Elt`: Element connectivity matrix
* `constrained_dofs`: Constrained degrees of freedom
* `free_dofs`: Free degrees of freedom
"""
@with_kw struct MdofMesh
    Elt::Matrix{Int}
    constrained_dofs::Vector{Int}
    free_dofs::Vector{Int}

    function MdofMesh(model::Mdof, bc = :CC)
        Nelt = length(model.k)
        Nnodes = length(model.m)
        Elt = undefs(Nelt, 3)

        for i = 1:Nelt
            Elt[i, 1] = i
            Elt[i, 2] = i
            Elt[i, 3] = i + 1
        end

        dofs = 1:Nnodes
        if bc == :CC
            constrained_dofs = [1, Nnodes]
        elseif bc == :CF
            constrained_dofs = [1]
        elseif bc == :FF
            constrained_dofs = []
        end

        free_dofs = setdiff(dofs, constrained_dofs)

        return new(Elt, constrained_dofs, free_dofs)
    end
end

"""
    assembly(model::Mdof)

Assembly of the mass, stiffness and damping matrices of a mdof system

**Input**
* `model`: Mdof model

**Outputs**
* `K`: Stiffness matrix
* `M`: Mass matrix
* `C`: Damping matrix (if viscous dampers are defined in `model`)
"""
function assembly(model::Mdof)
    (; k, m, c) = model

    Nm = length(m)
    Nk = length(k)
    Nc = length(c)

    Nk != Nm - 1 ? error("The number of masses must be equal to the number of springs plus 1") : nothing

    M = Diagonal(m)
    K = zeros(Nm, Nm)

    if Nc == 0
        flag = false
    else
        flag = true

        Nc != Nm - 1 ? error("The number of masses must be equal to the number of dampers plus 1") : nothing

        Nk != Nc ? error("The number of springs must be equal to the number of dampers") : nothing

        C = zeros(Nm, Nm)
    end

    for (i, ki) in enumerate(k)
        K[i:i+1, i:i+1] .+= element_matrix(ki)

        flag ? C[i:i+1, i:i+1] .+= element_matrix(c[i]) : nothing
    end

    return flag ? (K, M, C) : (K, M)
end

"""
    element_matrix(coeff::Float64)

Elemental stiffness or damping matrix

**Input**
* `coeff`: Stiffness or damping coefficient

**Output**
* `Kₑ`: Elemental stiffness or damping matrix
"""
element_matrix(coeff::Float64) = coeff.*[1. -1.; -1. 1.]