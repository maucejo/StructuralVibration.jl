function build_mesh end
function deformed_grid end
function detrend_mesh end
function renumber_element_connectivity end

"""
    detrend_plane(nodes, plane = :xy)

Remove the planar trend from a set of 3D points based on the specified plane.

**Inputs**
- `nodes::Matrix{Float64}`: The original node coordinates.
- `plane::Symbol`: The plane to detrend against (:xy, :xz, or :yz).

**Outputs**
- `detrended_nodes::Matrix{Float64}`: The detrended node coordinates.
"""
function detrend_plane(nodes::Matrix{T}, plane = :xy) where T <: Real
    plane in (:xy, :xz, :yz) || throw(ArgumentError("Invalid plane specified. Choose from :xy, :xz, or :yz."))

    perm = if plane == :xy
        [1, 2, 3]
    elseif plane == :xz
        [1, 3, 2]
    elseif plane == :yz
        [2, 3, 1]
    end
    x, y, z = eachcol(nodes[:, perm])

    A = [x y one.(x)]
    coeffs = A \ z

    trend = A * coeffs
    detrended_z = z .- trend

    detrended_nodes = copy(nodes)
    detrended_nodes[:, perm[end]] .= detrended_z

    return detrended_nodes
end