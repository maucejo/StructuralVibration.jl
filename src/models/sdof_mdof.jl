"""
    Sdof(m, ω₀, ξ)

Structure containing the data of a sdof system

# Constructor
* `m`: Mass [kg]
* `f₀`: Natural frequency [rad/s]
* `ξ`: Damping ratio

# Fields
* `m`: Mass [kg]
* `ω₀`: Natural frequency [rad/s]
* `ξ`: Damping ratio
"""
@with_kw struct Sdof
    m :: Float64
    ω₀ ::Float64
    ξ :: Float64

    Sdof(m, f₀, ξ) = new(m, 2π*f₀, ξ)
end

"""

    Mdof(k, m, c)

Structure containing the data for building a mdof system

# Fields
* `k`: Stiffness coefficients of the spring elements
* `m`: Masses of the mdof system
* `c`: Damping coefficients of the dampers
"""
@with_kw struct Mdof
    k::Vector{Float64}
    m::Vector{Float64}
    c::Vector{Float64}

    Mdof(k, m, c = Float64[]) = new(k, m, c)
end

struct MdofMesh
    Elt::Matrix{Int}
    constrained_dofs::Vector{Int}
    free_dofs::Vector{Int}

    function MdofMesh(model::Mdof, Nelt, bc = :CC)
        Nnodes = Nelt + 1
        Elt = undefs(Nelt, 3)
        Nnodes = length(Mdof.m)

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

        free_dofs = setdiff(collect(1:2Nnodes), constrained_dofs)

        new(Elt, constrained_dofs, free_dofs)
    end
end

"""
    assembly(model::Mdof)

Assembly of the mass, stiffness and damping matrices of a mdof system

# Input
* `model`: Mdof model
* `mesh`: OneDMesh

# Output
* `K`: Stiffness matrix
* `M`: Mass matrix
* `C`: Damping matrix (if dampers are present)
"""
function assembly(model::Mdof)
    (; k, m, c) = model

    Nm = length(m)
    Nk = length(k)
    Nc = length(c)

    if Nk != Nm - 2
        error("The number of masses must be equal to the number of springs plus 2")
    end

    M = Diagonal(m)
    K = undefs(Nm, Nm)

    if Nc == 0
        flag = false
    else
        flag = true

        if Nc != Nm - 2
            error("The number of masses must be equal to the number of springs plus 2")
        end

        if Nk != Nc
            error("The number of springs must be equal to the number of dampers")
        end

        C = undefs(Nm, Nm)
    end

    for (i, ki) in enumerate(k)
        K[i:i+1, i:i+1] .= element_matrix(ki)
        if flag
            C[i:i+1, i:i+1] .= element_matrix(c[i])
        end
    end

    if flag
        return K, M, C
    else
        return K, M
    end
end

"""
    element_matrix(coeff::Float64)

Elemental stiffness or damping matrix

# Input
* `coeff`: Stiffness or damping coefficient

# Output
* `Kₑ`: Elemental stiffness or damping matrix
"""
element_matrix(coeff::Float64) = coeff.*[1. -1.; -1. 1.]