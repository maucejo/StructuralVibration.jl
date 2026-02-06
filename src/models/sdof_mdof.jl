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
@show_data struct Sdof{T <: Real}
    m::T
    ω0::T
    ξ::T

    Sdof(m::T, f0::T, ξ::T) where T = new{T}(m, 2π*f0, ξ)
end

"""
    Mdof(k, m, c = Float64[])

Structure containing the data for building a mdof system

**Fields**
* `k`: Stiffness coefficients of the spring elements
* `m`: Masses of the mdof system
* `c`: Damping coefficients of the viscous dampers
"""
@show_data struct Mdof{T <: Real}
    k::Vector{T}
    m::Vector{T}
    c::Vector{T}

    function Mdof(k::Vector{T}, m::Vector{T}, c::Vector{T} = T[]) where T
        nk = length(k)
        nm = length(m)
        nc = length(c)

        nk != nm - 1 ? throw(DimensionMismatch("The number of masses must be equal to the number of springs plus 1")) : nothing

        if nc != 0
            nc != nm - 1 ? throw(DimensionMismatch("The number of masses must be equal to the number of dampers plus 1")) : nothing

            nk != nc ? throw(DimensionMismatch("The number of springs must be equal to the number of dampers")) : nothing
        end

        return new{T}(k, m, c)
    end
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
@show_data struct MdofMesh
    Elt::Matrix{Int}
    constrained_dofs::Vector{Int}
    free_dofs::Vector{Int}

    function MdofMesh(model::Mdof, bc = :CC)
        nelt = length(model.k)
        nnodes = length(model.m)

        dofs = 1:nnodes
        Elt = similar(dofs, nelt, 3)

        for i = 1:nelt
            Elt[i, 1] = i
            Elt[i, 2] = i
            Elt[i, 3] = i + 1
        end

        if bc == :CC
            constrained_dofs = [1, nnodes]
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

    nm = length(m)
    nk = length(k)
    nc = length(c)

    nk != nm - 1 ? error("The number of masses must be equal to the number of springs plus 1") : nothing

    M = Diagonal(m)
    K = zeros(eltype(k), nm, nm)

    if nc == 0
        flag = false
    else
        flag = true
        C = zeros(eltype(k), nm, nm)
    end

    for (i, ki) in enumerate(k)
        K[i:i+1, i:i+1] .+= element_matrix(ki)

        flag ? C[i:i+1, i:i+1] .+= element_matrix(c[i]) : nothing
    end

    # return flag ? (K, M, C) : (K, M)
    return K, M, C
end

"""
    element_matrix(coeff)

Elemental stiffness or damping matrix

**Input**
* `coeff`: Stiffness or damping coefficient

**Output**
* `Ke`: Elemental stiffness or damping matrix
"""
element_matrix(coeff) = coeff.*[1. -1.; -1. 1.]