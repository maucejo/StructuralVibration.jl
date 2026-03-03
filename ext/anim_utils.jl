"""
    create_figure(; title, xlabel, ylabel, zlabel)

Create a 3D figure with specified axis labels.

**Inputs**
- `xlabel::String`: Label for the x-axis (default is "x").
- `ylabel::String`: Label for the y-axis (default is "y").
- `zlabel::String`: Label for the z-axis (default is "z").

**Outputs**
- `fig::Figure`: The created figure.
- `ax::Axis3`: The 3D axis of the figure.
"""
function create_figure(; xlabel = "x", ylabel = "y", zlabel = "z")
    fig = Figure(figure_padding = 35)
    ax = Axis3(
                fig[1, 1],
                aspect = :data,
                xlabel = xlabel,
                ylabel = ylabel,
                zlabel = zlabel,
                xticksvisible = false,
                yticksvisible = false,
                zticksvisible = false,
            )

    hidespines!(ax)

    return fig, ax
end

"""
    viz_mesh(mesh; xlabel, ylabel, zlabel, title, color, colormap, alpha, showsegments, segmentsize, segmentcolor, xlim, ylim, zlim)

Visualize a 3D mesh with customizable options.

**Inputs**
- `mesh`: The mesh to visualize.
- `xlabel::String`: Label for the x-axis (default is "x").
- `ylabel::String`: Label for the y-axis (default is "y").
- `zlabel::String`: Label for the z-axis (default is "z").
- `title::String`: Title of the plot (default is an empty string).
- `color`: Color mapping for the mesh (default is :white).
- `colormap::Symbol`: Colormap to use for coloring the mesh (default is :jet1).
- `alpha::Real`: Transparency level of the mesh (default is 0.5).
- `showsegments::Bool`: Whether to show element segments (default is true).
- `segmentsize::Real`: Size of the segments if shown (default is 0.5).
- `segmentcolor`: Color of the segments if shown (default is :black).
- `xlim::Tuple{Real, Real}`: Limits for the x-axis (default is no limits).
- `ylim::Tuple{Real, Real}`: Limits for the y-axis (default is no limits).
- `zlim::Tuple{Real, Real}`: Limits for the z-axis (default is no limits).

**Output**
- `fig::Figure`: The figure containing the visualized mesh.
"""
function viz_mesh(mesh; xlabel = "x", ylabel = "y", zlabel = "z", title = "",color = :white, colormap = :jet1, alpha = 0.5, showsegments = true, segmentsize = 0.5, segmentcolor = :black, xlim = nothing, ylim = nothing, zlim = nothing)

    fig, ax = create_figure(xlabel = xlabel, ylabel = ylabel, zlabel = zlabel)
    ax.title = title
    viz!(ax, mesh, color = color, colormap = colormap, alpha = alpha, showsegments = showsegments, segmentsize = segmentsize, segmentcolor = segmentcolor)

    if xlim !== nothing
        xlims!(ax, xlim[1], xlim[2])
    end
    if ylim !== nothing
        ylims!(ax, ylim[1], ylim[2])
    end
    if zlim !== nothing
        zlims!(ax, zlim[1], zlim[2])
    end

    return fig
end

"""
    animate_mesh(nodes, elts, values; title, xlabel, ylabel, zlabel, color, colormap, alpha, showsegments, scale_factor, frequency, framerate, filename, xlim, ylim, zlim)

    animate_mesh(nodes, elts, values, x; xlabel, ylabel, zlabel, color, colormap, alpha, showsegments, unit_x, quantity_y, scale_factor, framerate, title_update_framerate, filename, xlim, ylim, zlim)

Create an animation of a deformed mesh.

**Inputs**
- `nodes::Matrix{Real}`: Original node coordinates.
- `elts::Vector{Vector{Int}}`: Element connectivity.
- `values::VecOrMat{Real}`: Values to apply to the nodes.
- `x::AbstractVector{Real}`: Parameter values corresponding to each set of values (e.g., time or frequency).
- `title::String`: Title of the animation (default is "Animated mesh").
- `xlabel::String`: Label for the x-axis (default is "x").
- `ylabel::String`: Label for the y-axis (default is "y").
- `zlabel::String`: Label for the z-axis (default is "z").
- `color`: Color mapping for the mesh (default is the values).
- `colormap::Symbol`: Colormap to use for coloring the mesh (default is :jet1).
- alpha::Real: Transparency level of the mesh (default is 1.0).
- `showsegments::Bool`: Whether to show element segments (default is true).
- `segmentsize::Real`: Size of the segments if shown (default is 0.5).
- `segmentcolor`: Color of the segments if shown (default is :black).
- `scale_factor::Real`: Factor to scale the values (default is 100).
- `framerate::Int`: Frame rate of the animation in frames per second (default is 35).
- `frequency::Real`: Frequency of the animated motion in Hz (default is framerate).
- `title_update_framerate::Int`: How often (in frames) to update the title of the animation (default is 10).
- `filename::String`: Name of the output video file (default is "animate_mesh.mp4").
- `xlim::Tuple{Real, Real}`: Limits for the x-axis (default is no limits).
- `ylim::Tuple{Real, Real}`: Limits for the y-axis (default is no limits).
- `zlim::Tuple{Real, Real}`: Limits for the z-axis (default is no limits).

**Outputs**
- Animation saved as a video file showing the deformed mesh over time.

**Notes**
- When `values` is a vector of values and frequency is specified, the animation simulates a harmonic motion of the mesh based on the provided values and frequency. This is useful for visualizing mode shapes or operational deflection shapes in structural vibration analysis.

- When `values` is a matrix and `x` is provided, the animation shows the evolution of the mesh deformation as a function of the parameter `x`, which can represent time, frequency, or any other relevant parameter. The title of the animation updates dynamically every 10 frames to reflect the current value of `x`.
"""
function animate_mesh(nodes::Matrix{Td}, elts::Vector{Vector{Tn}}, values::Vector{Td}; title = "Animated mesh", xlabel = "x", ylabel = "y", zlabel = "z", color = nothing, colormap = :jet1, alpha = 1., showsegments = true, segmentsize = 0.5, segmentcolor = :black, scale_factor = 100., framerate::Int = 35, frequency = float(framerate), filename = "animate_mesh.mp4", xlim = nothing, ylim = nothing, zlim = nothing) where {Tn <: Int, Td <: Real}

    # Static mesh
    points, tris, quads = construct_mesh(nodes, elts)
    smesh = SimpleMesh(points, [tris; quads])

    # Animation info
    ω = 2π*frequency
    t = range(0., 1., length = 500)
    y = cos.(ω*t)

    fig, ax = create_figure(xlabel = xlabel, ylabel = ylabel, zlabel = zlabel)
    ax.title = title

    # Set axis limits if provided
    if !isnothing(xlim)
        xlims!(ax, xlim[1], xlim[2])
    end
    if !isnothing(ylim)
        ylims!(ax, ylim[1], ylim[2])
    end
    if !isnothing(zlim)
        zlims!(ax, zlim[1], zlim[2])
    end

    # Determine color data (n values)
    n = size(nodes, 1)
    if isnothing(color)
        # Extract magnitude from values (3n vector)
        color_data = [norm(values[3*i-2:3*i]) for i in 1:n]
    else
        if color isa Vector
            n == length(color) || throw(ArgumentError("Color vector length must match number of nodes"))
        end
        color_data = color
    end

    # Create observables
    dmesh = Observable(smesh)
    col = Observable(y[1].*color_data)

    # Initial frame
    viz!(ax, dmesh, color = col, colormap = colormap, alpha = alpha, showsegments = showsegments, segmentsize = segmentsize, segmentcolor = segmentcolor)

    record(fig, filename, eachindex(t), framerate = framerate) do i
        # Deformed mesh at time t[i]
        dpoints = deformed_grid(nodes, y[i].*values, scale_factor)
        col[] = y[i].*color_data
        dmesh[] = SimpleMesh(dpoints, [tris; quads])
    end
end

function animate_mesh(nodes::Matrix{Tv}, elts::Vector{Vector{Tn}}, values::Matrix{Tv}, x::AbstractVector{Tv}; xlabel = "x", ylabel = "y", zlabel = "z", color = nothing, colormap = :jet1, alpha = 1., showsegments = true, segmentsize = 0.5, segmentcolor = :black, unit_x = "s", quantity_y = "Values", scale_factor = 100, framerate = 35, title_update_framerate = 10, filename = "animate_mesh.mp4", xlim = nothing, ylim = nothing, zlim = nothing) where {Tn <: Int, Tv <: Real}

    # Static mesh
    points, tris, quads = construct_mesh(nodes, elts)
    smesh = SimpleMesh(points, [tris; quads])

    fig, ax = create_figure(xlabel = xlabel, ylabel = ylabel, zlabel = zlabel)

    # Set axis limits if provided
    if !isnothing(xlim)
        xlims!(ax, xlim[1], xlim[2])
    end
    if !isnothing(ylim)
        ylims!(ax, ylim[1], ylim[2])
    end
    if !isnothing(zlim)
        zlims!(ax, zlim[1], zlim[2])
    end

    # Determine color data (n values)
    n = size(nodes, 1)
    nx = length(x)
    if isnothing(color)
        # Extract magnitude from values (3n vector)
        color_data = [norm(values[3*i-2:3*i, 1]) for i in 1:n]
    else
        color_data = if color isa Matrix
            nr, nc = size(color)
            nr == n || throw(ArgumentError("Color matrix row count must match number of nodes"))
            nc == length(x) || throw(ArgumentError("Color matrix column count must match length of x"))
            color[:, 1]
        else
            color
        end
    end

    # Create observables
    dmesh = Observable(smesh)
    col = Observable(color_data)
    function update_title!(ax, val)
         ax.title = quantity_y*" @ "*string(round(val, digits = 4))*" "*unit_x
    end
    title_obs = Observable(update_title!(ax, x[1]))

    # Initial frame
    viz!(ax, dmesh, color = col, colormap = colormap, alpha = alpha, showsegments = showsegments, segmentsize = segmentsize, segmentcolor = segmentcolor)

    record(fig, filename, eachindex(x), framerate = framerate) do i
        # Deformed mesh at time t[i]
        dpoints = deformed_grid(nodes, values[:, i], scale_factor)
        col[] = values[:, i]
        if rem(i, title_update_framerate) == 0
            title_obs[] = update_title!(ax, x[i])
        end
        dmesh[] = SimpleMesh(dpoints, [tris; quads])
    end
end