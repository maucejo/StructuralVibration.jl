using StructuralVibration
@usingany LazyGrids, Meshes, CairoMakie

## Helper function to create a mesh
function create_quad_mesh(xp::Vector, yp::Vector, nx::Int, ny::Int)
    # Number of quadrangles
    n_quads = (nx - 1) * (ny - 1)
    quads = Vector{Vector{Int}}(undef, n_quads)

    quad_idx = 1
    for j in 1:(ny-1)
        for i in 1:(nx-1)
            # Indices of the nodes of the quadrangle in counter-clockwise order
            n1 = (j-1) * nx + i
            n2 = (j-1) * nx + i + 1
            n3 = j * nx + i + 1
            n4 = j * nx + i

            quads[quad_idx] = [n1, n2, n3, n4]
            quad_idx += 1
        end
    end

    return quads
end

## Create a plate
# Dimensions
Lp = 0.6
bp = 0.4
hp = 1e-3

# Material parameters
E = 2.1e11
ρ = 7800.
ν = 0.33

plate = Plate(Lp, bp, hp, E, ρ, ν)

## Create a grid
nx = 10
ny = 10
x, y = ndgrid(range(0., Lp, nx), range(0., bp, ny))
xp = x[:]
yp = y[:]
zp = zeros(length(xp))

nodes = [xp yp zp]
elts = create_quad_mesh(xp, yp, nx, ny)
mesh = build_mesh(nodes, elts)

fig_mesh = viz_mesh(mesh, zlim = (-0.15, 0.15), title = "Mesh of the plate")

## Compute the second mode shape
ms = sin.(2π*xp/Lp) .* sin.(π*yp/bp)

## Compute deformed grid
ms_deformed = zeros(3size(nodes, 1))
ms_deformed[3:3:end] .= ms

dpoints = deformed_grid(nodes, ms_deformed, 0.15)

fig_deformed = viz_mesh(SimpleMesh(dpoints, mesh.topology.connec), color = ms, colormap = :jet1, alpha = 1., zlim = (-0.15, 0.15), title = "Deformed mesh of the plate")

## Animate the mesh
dpoints_poster = deformed_grid(nodes, ms_deformed, 0.1)
poster_fig = viz_mesh(SimpleMesh(dpoints_poster, mesh.topology.connec), color = ms, colormap = :jet1, alpha = 1., zlim = (-0.2, 0.2), title = "Second mode shape of a simply supported rectangular plate")

animate_mesh(nodes, elts, ms_deformed, scale_factor = 0.1, zlim = (-0.2, 0.2), title = "Second mode shape of a simply supported rectangular plate", framerate = 24, filename = "animated_plate.mp4", color = ms)