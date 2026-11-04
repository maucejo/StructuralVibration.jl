"""
    construct_mesh(nodes, elts)

Construct a mesh from node coordinates and element connectivity.

**Inputs**
- `nodes::Matrix{Real}`: A matrix of node coordinates (n_nodes x 3).
- `elts::Vector{Vector{Int}}`: A vector of element connectivity, where each element is a vector of node indices.

**Outputs**
- `points`: A vector of `Meshes.Point` objects representing the node coordinates.
- `tris`: A vector of `Meshes.Triangle` objects representing the triangular elements.
- `quads`: A vector of `Meshes.Quadrangle` objects representing the quadrilateral elements.
"""
function construct_mesh(nodes, elts)

    points = [Meshes.Point(node...) for node in eachrow(nodes)]
    connectivity = Tuple.(elts)

    nnode_per_elt = length.(elts)

    faces_tri = connectivity[nnode_per_elt .== 3]
    faces_quad = connectivity[nnode_per_elt .== 4]

    tris = connect.(faces_tri, Meshes.Triangle)
    quads = connect.(faces_quad, Meshes.Quadrangle)

    return points, tris, quads
end