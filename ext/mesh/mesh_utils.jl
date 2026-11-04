"""
    construct_mesh(nodes, elts)

Generate mesh points and connectivity from node coordinates and element connectivity.

**Inputs**
- `nodes::Matrix{Float64}`: The node coordinates.
- `elts::Vector{Vector{Int}}`: The element connectivity.

**Outputs**
- `points::Vector{Point}`: A list of mesh points.
- `tris::Vector{Triangle}`: A list of triangular elements.
- `quads::Vector{Quadrangle}`: A list of quadrilateral elements.

**Notes**
- This function is restricted to 2D meshes where elements are either triangles or quadrilaterals.
"""
function build_mesh(nodes, elts)
    points, tris, quads = construct_mesh(nodes, elts)
    return SimpleMesh(points, [tris; quads])
end

"""
    deformed_grid(nodes, displacements; scale_factor = 20)

Generate a deformed grid of points by applying displacements to the original node coordinates.

**Inputs**
- `nodes::Matrix{Float64}`: The original node coordinates.
- `values::Matrix{Float64}`: The values to apply to the nodes.
- `scale_factor::Float64`: A factor to scale the values (default is 20).

**Outputs**
- `deformed_nodes::Vector{Point}`: A list of deformed points.
"""
function deformed_grid(nodes, values, scale_factor = 100)
    nn = size(nodes, 1)
    val = permutedims(reshape(values, (3, nn)), [2, 1])
    p = @. nodes + scale_factor*val

    return [Point(point...) for point in eachrow(p)]
end

"""
    detrend_mesh(nodes, elts)

Remove global trends from a mesh defined by node coordinates and element connectivity.

**Inputs**
- `nodes::Matrix{Float64}`: The original node coordinates.
- `elts::Vector{Vector{Int}}`: The element connectivity.

**Outputs**
- `det_nodes::Matrix{Float64}`: The detrended node coordinates.
- `det_mesh::SimpleMesh`: The detrended mesh.
"""
function detrend_mesh(nodes, elts)

    # Remove global trend
    det_nodes = nodes .- mean(nodes, dims = 1)  # Centering the nodes sample

    # Mesh alignment
    det_nodes .= detrend_plane(det_nodes, :xy)  # Detrending in the XY plane
    det_nodes .= detrend_plane(det_nodes, :yz)  # Detrending in the YZ plane
    det_nodes .= detrend_plane(det_nodes, :xz)  # Detrending in the ZX plane

    det_points, tris, quads = construct_mesh(det_nodes, elts)

    det_mesh = SimpleMesh(det_points, [tris; quads])

    return det_nodes, det_mesh
end

"""
    renumber_element_connectivity(elts, nodes_ID)

Renumber elements connectivity based on a given list of node IDs.

**Inputs**
- `elts::Vector{Vector{Int}}`: The original element connectivity.
- `nodes_ID::Vector{Int}`: The list of node IDs.

**Outputs**
- `elts_new::Vector{Vector{Int}}`: The renumbered element connectivity.

**Notes**
- This function assumes that each node ID in `nodes_ID` corresponds to a unique node in the original mesh and that the element connectivity in `elts` references these node IDs.
- The renumbering is necessary if the elements ID are not contiguous or if they do not start from 1, which can cause issues when using Meshes.jl for visualization or further processing.
"""
function renumber_element_connectivity(elts, nodes_ID)
    # Renumber elements based on nodes_ID
    elts_new = [zeros(Int, length(elt)) for elt in elts]
    for (i, elt) in enumerate(elts)
        for (j, nID) in enumerate(elt)
            new_nID = only(findall(nodes_ID .== nID))
            elts_new[i][j] = new_nID
        end
    end

    return elts_new
end