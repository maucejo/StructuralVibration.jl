function bode_plot(freq, y, lw = 1.; colorline = :blue, xlab = "Frequency (Hz)", ylab = "Magnitude (dB)", title = "Bode Plot", xscale = :log, axis_tight = false, isdeg = false, layout = :vertical; ref_dB = 1.)

    fig = Figure()
    if layout == :vertical
        ax1 = Axis(fig[1,1])
        ax2 = Axis(fig[2,1])
        ax2.xlabel = xlab
        ax2.ylabel = ylab
    else
        fig = Figure()
        ax1 = Axis(fig[1,1])
        ax2 = Axis(fig[1,2])
        ax1.xlabel = xlab
        ax1.ylabel = ylab
        ax2.xlabel = xlab
        ax2.ylabel = ylab
    end
    lines!(ax, freq, y, color = colorline, linewidth = lw)
    ax.title = title
    return fig
end