"""
    fast_cmap()

Return the Paraview Fast Colormap as a color gradient.

**Output**
- `fast_cmap::ColorGradient`: The Paraview Fast Colormap as a color gradient that can be used for visualizations.

**Reference**
[1] F. Samsel, W.A. Scott, K. Moreland, A New Default Colormap for ParaView, in IEEE Computer Graphics and Applications, vol. 44, no. 04, pp. 150-160, 2024.

**Link to colormap json file**

https://github.com/kennethmoreland-com/kennethmoreland-com.github.io/blob/f0f50dcfda94c8e4c9c4927f78a98ff419fea450/content/color-advice/fast/fast.json
"""
function fast_cmap()
    # Paraview Fast Colormap
    fast_colors = [
        RGBf(0.056399999999999992,0.056399999999999992, 0.46999999999999997),
        RGBf(0.24300000000000013, 0.46035000000000043, 0.81000000000000005),
        RGBf(0.35681438265435211, 0.74502464853631423, 0.95436770289372197),
        RGBf(0.68820000000000003, 0.93000000000000005, 0.91790999999999989),
        RGBf(0.89949595512059022, 0.944646394975174, 0.7686567142818399),
        RGBf(0.92752075996107142, 0.62143890917391775, 0.31535705838676426),
        RGBf(0.80000000000000004, 0.35200000000000009, 0.15999999999999998),
        RGBf(0.58999999999999997, 0.076700000000000129, 0.11947499999999994)
    ]

    return fast_colors
end

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
- `colormap::Symbol`: Colormap to use for coloring the mesh (default is `fast_cmap()`).
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
function viz_mesh(mesh; xlabel = "x", ylabel = "y", zlabel = "z", title = "",color = :white, colormap = fast_cmap(), alpha = 1., showsegments = true, segmentsize = 0.5, segmentcolor = :black, xlim = nothing, ylim = nothing, zlim = nothing)

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
- `colormap::Symbol`: Colormap to use for coloring the mesh (default is `fast_cmap()`).
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
function animate_mesh(nodes::Matrix{Tn}, elts::Vector{Vector{Te}}, values::Vector{Tv}; title = "Animated mesh", xlabel = "x", ylabel = "y", zlabel = "z", color = nothing, colormap = fast_cmap(), alpha = 1., showsegments = true, segmentsize = 0.5, segmentcolor = :black, scale_factor = 100., framerate::Int = 35, frequency = float(framerate), filename = "animate_mesh.mp4", xlim = nothing, ylim = nothing, zlim = nothing) where {Tn <: Real, Te <: Int, Tv <: Real}

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
        color_data = [mean(values[3*i-2:3*i]) for i in 1:n]
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

function animate_mesh(nodes::Matrix{Tn}, elts::Vector{Vector{Te}}, values::Matrix{Tv}, x::AbstractArray; xlabel = "x", ylabel = "y", zlabel = "z", color = nothing, colormap = fast_cmap(), alpha = 1., showsegments = true, segmentsize = 0.5, segmentcolor = :black, unit_x = "s", quantity_y = "Values", scale_factor = 100, framerate = 35, title_update_framerate = 10, filename = "animate_mesh.mp4", xlim = nothing, ylim = nothing, zlim = nothing) where {Tn <: Real, Te <: Int, Tv <: Real}

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
        color_data = permutedims(reduce(hcat, [vec(mean(values[3*i-2:3*i, :], dims = 1)) for i in 1:n]), [2, 1])
    else
        color_data = if color isa Matrix
            nr, nc = size(color)
            nr == n || throw(ArgumentError("Color matrix row count must match number of nodes"))
            nc == nx || throw(ArgumentError("Color matrix column count must match length of x"))
            color
        else
            color
        end
    end

    # Create observables
    dmesh = Observable(smesh)
    if color_data isa Matrix
        col = Observable(color_data[:, 1])
    else
        col = Observable(color_data)
    end

    function update_title!(ax, val)
         ax.title = quantity_y*" @ "*string(round(val, digits = 4))*" "*unit_x
    end
    title_obs = Observable(update_title!(ax, x[1]))

    # Initial frame
    viz!(ax, dmesh, color = col, colormap = colormap, alpha = alpha, showsegments = showsegments, segmentsize = segmentsize, segmentcolor = segmentcolor)

    record(fig, filename, eachindex(x), framerate = framerate) do i
        # Deformed mesh at time t[i]
        dpoints = deformed_grid(nodes, values[:, i], scale_factor)

        if color_data isa Matrix
            col[] = color_data[:, i]
        else
            col[] = color_data
        end

        if rem(i, title_update_framerate) == 0
            title_obs[] = update_title!(ax, x[i])
        end
        dmesh[] = SimpleMesh(dpoints, [tris; quads])
    end
end