"""
    OneDMesh(model, xmin, Nelt, bc)

Construct a mesh for a beam with Nelt elements, length L and starting at xmin.

**Constructor parameters**
* `model`: Structure containing the data related to the 1D system
* `xmin`: starting position of the beam
* `Nelt`: number of elements
* `bc`: Boundary conditions type
    * `:CC`: Clamped - Clamped
    * `:FF`: Free - Free
    * `:CF`: Clamped - Free
    * `:SS`: Simply Supported - Simply Supported (specific to beam)
    * `:CS`: Clamped - Simply Supported (specific to beam)
    * `:SF`: Simply Supported - Free (specific to beam)

**Fields**
* `xmin`: Starting position of the beam
* `L`: Length of the beam
* `Nodes`: Nodes of the mesh
* `Elt`: Elements of the mesh
* `Ndof_per_node``: Number of degrees of freedom per node
* `elem_size`: Size of the elements
* `constrained_dofs`: Constrained degrees of freedom
* `free_dofs`: Free degrees of freedom
"""
@show_data struct OneDMesh{T <: Real}
    xmin::T
    L::T
    Nodes::Matrix{T}
    Elt::Matrix{Int}
    Ndof_per_node::Int
    elem_size::T
    constrained_dofs::Vector{Int}
    free_dofs::Vector{Int}

    function OneDMesh(model::OneDStructure, xmin::T, Nelt::Int, bc::Symbol = :CC) where T
        Nnodes = Nelt + 1
        Nodes = undefs(Nnodes, 2)
        Elt = undefs(Nelt, 3)
        elem_size = model.L/Nelt

        for i = 1:Nnodes
            Nodes[i, 1] = i
            Nodes[i, 2] = xmin + (i - 1)*elem_size
        end

        for i = 1:Nelt
            Elt[i, 1] = i
            Elt[i, 2] = i
            Elt[i, 3] = i + 1
        end

        if model isa Beam
            Ndof_per_node = 2
            dofs = 1:2Nnodes
            if bc == :SS
                constrained_dofs = [1, 2Nnodes - 1]
            elseif bc == :CC
                constrained_dofs = [1, 2, 2Nnodes - 1, 2Nnodes]
            elseif bc == :CS
                constrained_dofs = [1, 2, 2Nnodes - 1]
            elseif bc == :CF
                constrained_dofs = [1, 2]
            elseif bc == :SF
                constrained_dofs = [1]
            elseif bc == :FF
                constrained_dofs = Int[]
            else
                error("Boundary conditions not implemented")
            end
        elseif model isa WaveEquation
            Ndof_per_node = 1
            dofs = 1:Nnodes
            if bc == :CC
                constrained_dofs = [1, Nnodes]
            elseif bc == :CF
                constrained_dofs = [1]
            elseif bc == :FF
                constrained_dofs = Int[]
            else
                error("Boundary conditions not implemented")
            end
        end
        free_dofs = setdiff(dofs, constrained_dofs)

        new{T}(xmin, model.L, Nodes, Elt, Ndof_per_node, elem_size, constrained_dofs, free_dofs)
    end
end

"""
    assembly(model::OneDstructure, mesh::OneDMesh)

Compute the global stiffness and mass matrices for a 1D structure with a given mesh.

**Inputs**
* `model`: OneDStructure
    * `Beam`: Beam model
    * `WaveEquation`: Bar, Rod or String model
* `mesh`: OneDMesh

**Outputs**
* `K`: global stiffness matrix
* `M`: global mass matrix
"""
function assembly(model::OneDStructure, mesh::OneDMesh)
    # Compute elemental matrices
    ke, me = element_matrix(model, mesh.elem_size)

    (; Elt, Ndof_per_node) = mesh
    Nelt = size(Elt, 1)

    # Assemble global matrices
    Nddl = size(mesh.Nodes, 1)*Ndof_per_node

    K = zeros(Nddl, Nddl)
    M = zeros(Nddl, Nddl)
    ind = zeros(Int, 2Ndof_per_node)
    @inbounds @views for i = 1:Nelt
        ind .= (Ndof_per_node*Elt[i, 2:end]' .+ repeat((0:Ndof_per_node-1), 1, Ndof_per_node) .- Ndof_per_node .+ 1)[:];

        K[ind, ind] += ke
        M[ind, ind] += me
    end

    return K, M
end

"""
    element_matrix(model::Beam, h)
    element_matrix(model::WaveEquation, h)

Compute the elemental stiffness and mass matrices for a beam with a element size `h`.

**Inputs**
* `model`: Structure containing the data related to the 1D system
* `h`: element size

**Outputs**
* `ke`: elemental stiffness matrix
* `me`: elemental mass matrix
"""
function element_matrix(beam::Beam, h)
    # Constants
    kc = beam.D/h^3
    mc = beam.m*h/420.

    # Elemental stiffness matrix
    ke = kc.*[12. 6h -12. 6h;
              6h 4h^2 -6h 2h^2;
              -12. -6h 12. -6h;
              6h 2h^2 -6h 4h^2]

    # Elemental mass matrix
    me = mc.*[156. 22h 54. -13h;
              22h 4h^2 13h -3h^2;
              54. 13h 156. -22h;
              -13h -3h^2 -22h 4h^2]

    return ke, me
end

function element_matrix(we::WaveEquation, h)
    # Constants
    kc = we.D/h
    mc = we.m*h/6.

    # Elemental stiffness matrix
    ke = kc.*[1. -1.;
              -1. 1.]

    # Elemental mass matrix
    me = mc.*[2. 1.;
              1. 2.]

    return ke, me
end

"""
    selection_matrix(mesh, selected_dofs)

Compute the selection matrix for the selected dofs.

**Inputs**
* `mesh`: OneDmesh
* `selected_dofs`: Selected dofs


**Output**
- `S`: Selection matrix
"""
function selection_matrix(mesh::OneDMesh, selected_dofs)

    N = length(selected_dofs)
    S = zeros(N, length(mesh.free_dofs))
    for (i, dof) in enumerate(selected_dofs)
        pos = findall(mesh.free_dofs .== dof)
        if length(pos) != 0
            S[i, pos[1]] = 1.
        end
    end

    return S
end

"""
    apply_bc(A, mesh)

Apply boundary conditions to a given matrix

**Inputs**
* `A`: Matrix to apply the boundary conditions
* `mesh`: Mesh of the system

**Output**
* `A_bc`: Matrix with boundary conditions applied
"""
function apply_bc(A, mesh)

    (; free_dofs) = mesh

    return A[free_dofs, free_dofs]
end

"""
    eigenmode(K, M, n = size(K, 1))

Computes the eigenmodes of a system defined by its mass and stiffness matrices.

**Inputs**
* `K`: Stiffness matrix
* `M`: Mass matrix
* `n`: Number of modes to be keep in the modal basis

**Outputs**
* `ω`: Vector of natural frequencies
* `Φ`: Mass-normalized mode shapes
"""
function eigenmode(K, M, n = size(K, 1))

    λ, Φ = eigen(K, M)
    ω = @. √abs(λ[1:n])

    return ω, Φ[:, 1:n]
end

"""
    rayleigh_damping_matrix(K, M, α, β)
    rayleigh_damping_matrix(K, M, ω1, ω2, ξ1, ξ2)

Compute the Rayleigh damping matrix for a given stiffness and mass matrices

**Inputs**
* `K`: Stiffness matrix
* `M`: Mass matrix
* Construction parameters
    * Method 1
        * `α`: Mass proportional damping coefficient
        * `β`: Stiffness proportional damping coefficient
    * Method 2
        * `ω1`: First natural frequency
        * `ω2`: Second natural frequency
        * `ξ1`: Damping ratio for the first natural frequency
        * `ξ2`: Damping ratio for the second natural frequency

**Output**
* `C`: Rayleigh damping matrix
"""
function rayleigh_damping_matrix(K, M, α, β)
    return α*M + β*K
end

function rayleigh_damping_matrix(K, M, ω1, ω2, ξ1, ξ2)
    β = 2(ξ2*ω2 - ξ1*ω1)/(ω2^2 - ω1^2)
    α = 2ξ1*ω1 - ω1^2*β
    return α*M + β*K
end

"""
    modal_damping_matrix(M, ωn, ξn, Φn)

Compute the damping matrix C from modal parameters

**Inputs**
* `M`: Mass matrix
* `ωn`: Natural angular frequencies
* `ξn`: Damping ratios
* `Φn`: Mass-normalized mode shapes

**Output**
* `C`: Damping matrix
"""
function modal_damping_matrix(M, ωn, ξn, Φn)
    Cn = Diagonal(2ξn.*ωn)

    return M*Φn*Cn*Φn'*M
end
