"""
    theme_choice(name::Symbol; fonts = Makie.theme(:fonts),
                 titlesize = 20., labelsize = 18., ticklabelsize = 14.)

Choose the theme for the plots.

**Inputs**
* `name`: Name of the theme
    * `:makie`
    * `:sv`
* `fonts`: Fonts of the figure (default: Makie.theme(:fonts))
* `titlesize`: Title size (default: 20.)
* `labelsize`: Label size (default: 18.)
* `ticklabelsize`: Tick label size (default: 14.)

**Output**
* `theme`: Theme
"""
function theme_choice(name::Symbol; fonts = Makie.theme(:fonts), titlesize = 20., labelsize = 18., ticklabelsize = 14.)
    linestyle = [:solid, :dash, :dashdot, :dashdotdot, :dot]
    if name == :makie
        colorline = Makie.wong_colors()
    else name == :sv
        colorline = [:blue, :red, :green, :magenta, :orange, :cyan, :black]
    end

    theme = Theme(
        fonts = fonts,
        palette = (color = colorline, linestyle = linestyle),
        Lines = (cycle = Cycle([:color, :linestyle], covary = true), ),
        Axis = (
            titlesize = titlesize,
            xlabelsize = labelsize,
            xlabelfont = :bold,
            ylabelsize = labelsize,
            ylabelfont = :bold,
            xticklabelsize = ticklabelsize,
            yticklabelsize = ticklabelsize,
        ),
        Axis3 = (
            xlabelsize = labelsize,
            xlabelfont = :bold,
            ylabelsize = labelsize,
            ylabelfont = :bold,
            zlabelsize = labelsize,
            zlabelfont = :bold,
            xticklabelsize = ticklabelsize,
            yticklabelsize = ticklabelsize,
            zticklabelsize = ticklabelsize,
        )
    )

    return theme
end

"""
    sv_plot(x, y...; lw = 1., xscale = identity, yscale = identity,
            axis_tight = true, title = " ", xlabel = "x", ylabel = "y",
            legend = (active = false, position = :right, orientation = :vertical, entry = " "))

Plot a 2D plot.

**Inputs**
* `x`: x-axis values
* `y`: y-axis values
* `lw`: linewidth
* `xscale`: x-axis scale (default: identity)
* `yscale`: y-axis scale (default: identity)
* `axis_tight`: Tight axis (default: true)
* `title`: Title of the plot (default: " ")
* `xlabel`: x-axis label
* `ylabel`: y-axis label
* `legend`: Legend parameters
    * active : Bool
    * position : Symbol
    * entry : String

**Output**
* `fig`: Figure
"""
function sv_plot(x, y...; lw = 1., xscale = identity, yscale = identity, axis_tight = true, title = " ",  xlabel = "x", ylabel = "y", legend = (active = false, position = :rt, orientation = :vertical, entry = " "))

    # Some checks
    t = typeof(legend)

    ny = length(y)

    if legend.active
        if !hasfield(t, :position)
            legend_position = :rt
        else
            legend_position = legend.position
        end

        if !hasfield(t, :orientation)
            legend_orientation = :vertical
        else
            legend_orientation = legend.orientation
        end

        if !hasfield(t, :entry)
            if ny == 1
                legend_entry = ["Data 1"]
            else
                legend_entry = ["Data $i" for i in 1:ny]
            end
        else
            if length(legend.entry) != ny
                error("The number of entries in the legend must be equal to the number of rows in y")
            end

            legend_entry = legend.entry
        end

        leg = (active = legend.active, position = legend_position, entry = legend_entry, orientation = legend_orientation)
    else
        if ny == 1
            legend_entry = [" "]
        else
            legend_entry = [" " for _ in 1:ny]
        end

        leg = (active = false, entry = legend_entry)
    end


    fig = Figure()
    ax = Axis(fig[1,1], xlabel = xlabel, ylabel = ylabel, xscale = xscale, yscale = yscale, title = title)

    if ny == 1
        lines!(ax, x, y[1], linewidth = lw, label = leg.entry)
    else
        for (yi, labeli) in zip(y, leg.entry)
            lines!(ax, x, yi, linewidth = lw, label = labeli)
        end
    end

    if leg.active
        axislegend(ax, position = leg.position, backgroundcolor = (:white, 0.5), orientation = leg.orientation)
    end

    if axis_tight
        xlims!(ax, minimum(x), maximum(x))
    end

    return fig
end

"""
    stabilization_plot(stab::StabilizationAnalysis, indicator; display_poles)

Plot stabilization diagram for EMA-MDOF pole stability analysis.

**Inputs**
* `stab`: EMA-MDOF stabilization data
* `indicator`: Indicator to plot
    * `:psif` : Power spectrum indicator function (default)
    * `:cmif` : Complex mode indicator function
* `display_poles`: Vector of Bool to choose which poles to display
    * `display_poles[1]` : Stable in frequency and damping (default: true)
    * `display_poles[2]` : Stable in frequency but not stable in damping (default: true)
    * `display_poles[3]` : Not stable in frequency (default: true)

**Output**
* `fig`: Figure
"""
function stabilization_plot(stab::StabilizationAnalysis, indicator = :psif; display_poles = [true, true, true])
    # Extract data for the selected indicator
    (; prob, poles, modefn, mode_stabfn, mode_stabdr) = stab

    # FRF post-processing - Frequency range reduction
    if prob isa EMAProblem
        (; frf, freq, type_frf) = prob
    elseif prob isa OMAProblem
        frf = prob.fullspec
        freq = prob.freq
        indicator = :psif
    end

    # Indicator calculation
    if indicator == :psif
        indicator_data = psif(frf)
        indicator_name = "PSIF"
    elseif indicator == :cmif
        indicator_data = cmif(frf, type = type_frf)
        indicator_name = "CMIF"
    else
        throw(ArgumentError("Indicator not available. Available indicators are :psif and :cmif"))
    end

    # Stabilization diagram data
    Nmodes = length(poles[end])
    max_orders = 1:Nmodes
    model_orders = one.(max_orders)*max_orders'
    fn = modefn[:]
    orders = model_orders[:]

    # Poles
    fig = Figure()
    ax_poles = Axis(
        fig[2,1],
        yticklabelcolor = :blue,
        ylabelcolor = :blue,
        ytickcolor = :blue,
        xgridvisible = false,
        ygridvisible = false,
        xlabel = "Frequency (Hz)",
        ylabel = "Model order"
    )
    # scatter!(ax_poles, poles_scatter)
    if display_poles[1]
        # Stable in frequency and damping
        isstab = (mode_stabfn[:] .&& mode_stabdr[:])

        scatter!(ax_poles, Point2f.(fn[isstab], orders[isstab]), color = :green, marker = :star4, label = "Stable freq. & damp.")
    end

    if display_poles[2]
        # Stable in frequency but not stable in damping
        isstabf = @. (mode_stabfn[:] && !mode_stabdr[:])

        scatter!(ax_poles, Point2f.(fn[isstabf], orders[isstabf]), color = :blue, label = "Stable freq.")
    end

    if display_poles[3]
        # Not stable in frequency
        isnotstab = .!mode_stabfn[:]

        scatter!(ax_poles, Point2f.(fn[isnotstab], orders[isnotstab]), color = :red, marker = :xcross, label = "Not stable")
    end

    xlims!(ax_poles, minimum(freq), maximum(freq))

    # Indicator
    ax_indicator = Axis(
        fig[2,1],
        leftspinecolor = :blue,
        rightspinecolor = :red,
        yticklabelcolor = :red,
        ylabelcolor = :red,
        ytickcolor = :red,
        ylabel = indicator_name,
        yaxisposition = :right,
        xticklabelsvisible = false,
        xticksvisible = false,
        yscale = log10
    )

    linkxaxes!(ax_poles, ax_indicator)

    if indicator == :psif
        lines!(ax_indicator, freq, indicator_data, color = (:gray, 0.5))
    else
        for cmifk in eachrow(indicator_data)
            lines!(ax_indicator, freq, cmifk, color = (:gray, 0.25))
        end
    end
    xlims!(ax_indicator, minimum(freq), maximum(freq))

    Legend(fig[1, 1], ax_poles, orientation = :horizontal)

    return fig
end

"""
    peaksplot(x, y; width = 5, min_prom = 0., max_prom = Inf,
              xlabel = "x", ylabel = "y")

Plot data with detected peaks.

**Inputs**
* `x`: x-axis values
* `y`: y-axis values
* `width`: Half-width of the peaks (default: 1)
* `min_prom`: Minimum peak prominence (default: 0.)
* `max_prom`: Maximum peak prominence (default: Inf)
* `xlabel`: x-axis label
* `ylabel`: y-axis label

**Output**
* `fig`: Figure
"""
function peaksplot(x, y; width = 1, min_prom = 0., max_prom = Inf, xlabel = "x", ylabel = "y")
    pks = findmaxima(y, width)
    pks = peakproms!(pks, min = min_prom, max = max_prom) |> peakwidths!

    fig = Figure()
    ax = Axis(fig[1,1], xlabel = xlabel, ylabel = ylabel)
    lines!(ax, x, y)
    scatter!(ax, x[pks.indices], y[pks.indices], color = :red, marker = :star4)
    xlims!(ax, minimum(x), maximum(x))

    return fig
end

"""
    bode_plot(freq, y...; lw = 1., xlabel = "Frequency (Hz)", xscale = identity,
              axis_tight = true, isdeg = false, layout = :vertical,
              ref_dB = 1., legend = (active = false, position = :rt, entry = " "))

Plot Bode diagram of a frequency response or a FRF.

**Inputs**
* `freq`: Frequency range of interest
* `y`: Frequency response or FRF
* `lw`: Line width
* `xlabel`: x-axis label
* `xscale`: x-axis scale (default: :log)
* `axis_tight`: Tight axis (default: false)
* `isdeg`: Phase in degrees (default: false)
* `layout`: Layout of the plot (default: :vertical)
* `ref_dB`: Reference value for magnitude (default: 1.)
* `legend`: Legend parameters (default: (active = false, position = :rt, entry = " "))

**Output**
* `fig`: Figure
"""
function bode_plot(freq, y...; lw = 1., xlabel = "Frequency (Hz)", xscale = identity, axis_tight = true, isdeg = false, layout = :vert, ref_dB = 1., legend = (active = false, position = :rt, entry = " "))

    # Some checks
    ny = length(y)

    t = typeof(legend)
    if !hasfield(t, :active)
        error("legend must be a NamedTuple with at least the fields active")
    end

    if legend.active
        if !hasfield(t, :position)
            legend_position = :rt
        else
            legend_position = legend.position
        end

        if !hasfield(t, :entry)
            if ny == 1
                legend_entry = ["Data 1"]
            else
                legend_entry = ["Data $i" for i in 1:ny]
            end
        else
            entry = legend.entry
        end

        leg = (active = legend.active, position = legend_position, entry = legend_entry)
    else
        if ny == 1
            legend_entry = [" "]
        else
            legend_entry = [" " for _ in 1:ny]
        end

        leg = (active = false, entry = legend_entry)
    end

    ylab1 = "Magnitude (dB)"
    if isdeg
        ylab2 = "Phase (deg)"
    else
        ylab2 = "Phase (rad)"
    end

    fig = Figure()
    if layout == :vert
        ax1 = Axis(fig[1,1], xscale = xscale)
        ax2 = Axis(fig[2,1], xscale = xscale)
        ax1.ylabel = ylab1
        ax2.xlabel = xlabel
        ax2.ylabel = ylab2
    else
        ax1 = Axis(fig[1,1], xscale = xscale)
        ax2 = Axis(fig[1,2], xscale = xscale)
        ax1.xlabel = xlabel
        ax1.ylabel = ylab1
        ax2.xlabel = xlabel
        ax2.ylabel = ylab2
    end

    ymag = similar(freq)
    ϕ = similar(freq)
    if ny == 1
        @. ymag = abs(y[1])
        if isdeg
            @. ϕ = rad2deg(angle(y[1]))
        else
            @. ϕ = angle(y[1])
        end
        lines!(ax1, freq, 20log10.(ymag/ref_dB), linewidth = lw, label = leg.entry[1])
        lines!(ax2, freq, unwrap(ϕ), linewidth = lw)
    else
        for (yi, labeli) in zip(y, leg.entry)
            @. ymag = abs(yi)
            if isdeg
                @. ϕ = rad2deg(angle(yi))
            else
                @. ϕ = angle(yi)
            end
            lines!(ax1, freq, 20log10.(ymag/ref_dB), linewidth = lw, label = labeli)
            lines!(ax2, freq, unwrap(ϕ), linewidth = lw, label = labeli)
        end
    end

    if layout == :vert
        if leg.active
            axislegend(ax1, position = leg.position, backgroundcolor = (:white, 0.5))
        end
    else
        if leg.active
            axislegend(ax1, position = leg.position, backgroundcolor = (:white, 0.5))
            axislegend(ax2, position = leg.position, backgroundcolor = (:white, 0.5))
        end
    end

    if axis_tight
        xlims!(ax1, minimum(freq), maximum(freq))
        xlims!(ax2, minimum(freq), maximum(freq))
    end

    return fig
end

"""
    nyquist_plot(y)

Plot Nyquist diagram

**Inputs**
* `y`: Complex data vector

**Output**
* `fig`: Figure
"""
function nyquist_plot(y::Vector{T}) where {T <: Complex}
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Real part", ylabel = "Imaginary part", aspect = DataAspect())
    lines!(ax, real.(y), imag.(y))

    return fig
end

"""
    nyquist_plot(freq, y, xlab = "Frequency (Hz)";
                 projection = false)

Plot Nyquist diagram in 3D

**Inputs**
* `freq`: Frequency range
* `y`: Complex vector
* `ylabel`: y-axis label
* `projection`: Projection of the curve on the xy, yz, and xz planes (default: false)
    * on the xy plane: (freq, real(y))
    * on the yz plane: (imag(y), freq)
    * on the xz plane: (real(y), imag(y))

**Output**
* `fig`: Figure
"""
function nyquist_plot(freq, y::Vector{T}, ylabel = "Frequency (Hz)"; projection = false) where {T <: Complex}

    fig = Figure()
    ax = Axis3(fig[1,1], xlabel = "Real part", ylabel = ylabel, zlabel = "Imaginary part", aspect = (1, 2, 1))

    yr = real.(y)
    yi = imag.(y)
    lines!(ax, yr, freq, yi)

    # Some checks
    minyr, maxyr = extrema(yr)
    if minyr*maxyr > 0.
        if maxyr > 0.
            minyr -= 0.1maxyr
        else
            maxyr -= 0.1minyr
        end
    end

    minyi, maxyi = extrema(yi)
    if minyi*maxyi > 0.
        if maxyi > 0.
            minyi -= 0.1maxyi
        else
            maxyi -= 0.1minyi
        end
    end

    minf, maxf = extrema(freq)
    if minf == 0.
        minf = -0.1
        α = 0.9
    else
        minf -= 0.1
        α = 1.1
    end

    ax.yreversed = true

    if projection
        lines!(ax, yr, yi, color = :black, linewidth = 0.5, transformation = (:xz, α*minf))
        lines!(ax, freq, yi, color = :black, linewidth = 0.5, transformation = (:yz, 1.09maxyr))
        lines!(ax, yr, freq, color = :black, linewidth = 0.5, transformation = (:yx, 1.09minyi))
    end

    xlims!(ax, 1.1minyr, 1.1maxyr)
    ylims!(ax, maxf, minf)
    zlims!(ax, 1.1minyi, 1.1maxyi)

    ax.zticklabelpad = 5.

    return fig
end

"""
    waterfall_plot(x, y, z; zmin = minimum(z), lw = 1.,
                   colorline = :auto, colmap = :viridis, colorband = (:white, 1.),
                   xlabel = "x", ylabel = "y", zlabel = "z", edge = true,
                   axis_tight = false, xlim = [minimum(x), maximum(x)],
                   ylim = [minimum(y), maximum(y)], zlim = [zmin, maximum(z)])

Plot a waterfall plot.

**Inputs**
* `x`: x-axis values
* `y`: y-axis values
* `z`: z-axis values
* `zmin`: minimum value of z-axis
* `lw::Real`: linewidth
* `colorline`: color of the lines
* `colmap`: Name of the colormap
* `colorband`: Tuple defining the color of the band
    * color : Color
    * alpha : Alpha value for transparency
* `xlabel`: x-axis label
* `ylabel`: y-axis label
* `zlabel`: z-axis label
* `edge`: Display edges (default: true)
* `axis_tight`: Tight axis (default: false)
* `xlim`: x-axis limits
* `ylim`: y-axis limits
* `zlim`: z-axis limits

**Output**
* `fig`: Figure
"""
function waterfall_plot(x, y, z; zmin = minimum(z), lw = 1., colorline = :auto, colmap = :viridis, colorband = (:white, 1.), xlabel = "x", ylabel = "y", zlabel = "z", edge = true, axis_tight = false, xlim = [minimum(x), maximum(x)], ylim = [minimum(y), maximum(y)], zlim = [zmin, maximum(z)])

    # Initialization
    ny = length(y)
    I₂ = ones(2)

    fig = Figure()
    ax = Axis3(fig[1,1], xlabel = xlabel, ylabel = ylabel, zlabel = zlabel)
    for (j, yv) in enumerate(reverse(y))
        idz = ny - j + 1
        zj = z[idz, :]
        lower = Point3f.(x, yv, zmin)
        upper = Point3f.(x, yv, zj)
        band!(ax, lower, upper, color = colorband)

        if edge
            edge_start = [Point3f(x[1], yv, zmin), Point3f(x[1], yv, zj[1])]
            edge_end = [Point3f(x[end], yv, zmin), Point3f(x[end], yv, zj[end])]
        end

        if colorline == :auto
            lines!(ax, upper, color = zj, colormap = colmap, linewidth = lw)

            if edge
                lines!(ax, edge_start, color = zj[1]*I₂, colormap = colmap, linewidth = lw)
                lines!(ax, edge_end, color = zj[end]*I₂, colormap = colmap, linewidth = lw)
            end
        else
            lines!(ax, upper, color = colorline, linewidth = lw)

            if edge
                lines!(ax, edge_start, color = colorline, linewidth = lw)
                lines!(ax, edge_end, color = colorline, linewidth = lw)
            end
        end
    end

    if axis_tight
        xlims!(ax, minimum(x), maximum(x))
        ylims!(ax, minimum(y), 1.01*maximum(y))
        zlims!(ax, zmin, maximum(z))
    else
        xlims!(ax, xlim[1], xlim[2])
        ylims!(ax, ylim[1], ylim[2])
        zlims!(ax, zlim[1], zlim[2])
    end

    ax.zticklabelpad = 5.

    return fig
end