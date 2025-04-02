"""
    dofs_selection(mesh, X, dof_type)

Select the dofs corresponding to the closest nodes to the positions X.

# Inputs
- mesh: Structure mesh
- `X`: Selected positions
- `dof_type`: Type of the selected degree of freedom (only for beams)
    * :trans: Transverse displacement
    * :rot: Rotation

# Outputs
- `S`: Selection matrix
- `dofs`: dofs corresponding to the closest nodes

# Example
```julia-repl
julia> mesh = BeamMesh(0., 1., 10)
julia> dofs, S = dofs_selection(mesh, [0.1, 0.2])
```
"""
function dofs_selection(mesh :: Mesh, X, dof_type = :trans)
    N = length(X)
    dofs = zeros(Int, N)
    S = zeros(N, length(mesh.free_dofs))
    @inbounds for (i, Xi) in enumerate(X)
        d = @. abs(mesh.Nodes[:, 2] - Xi)
        if mesh.Ndof_per_node == 2
            if dof_type == :trans
                dofs[i] = 2argmin(d) - 1
            elseif dof_type == :rot
                dofs[i] = 2argmin(d)
            end
        elseif mesh.Ndof_per_node == 1
            dofs[i] = argmin(d)
        end

        pos = findall(mesh.free_dofs .== dofs[i])
        if length(pos) != 0
            S[i, pos[1]] = 1.
        end
    end

    return S, dofs
end

"""
    gradient(f::Vector{Float64}; h = 1.)

Compute the gradient of a vector `f` given a step size `h`.

# Inputs
- `f::Vector{Float64}`: Function values
- `h`: Step size

# Output
- `df`: Gradient of the vector `f`
"""
function gradient(f::Vector{Float64}; h = 1.)
    n = length(f)
    df = zeros(n)

    # First point
    df[1] = (f[2] - f[1])/h

    # Last point
    df[end] = (f[end] - f[end-1])/h

    # Central difference
    df[2:end-1] = (f[3:end] - f[1:end-2])/2h

    return df
end
