"""
    bode_plot(freq, y; lw = 1., colorline = :blue, xlab = "Frequency (Hz)", ylab = "Magnitude (dB)", title = "Bode Plot", xscale = :log, axis_tight = false, isdeg = false, layout = :vertical, ref_dB = 1.)

Plot Bode diagram of a frequency response or a FRF.

# Inputs
- `freq`: Frequency range of interest
- `y`: Frequency response or FRF
- `lw`: Line width
- `colorline`: Line color
- `xlab`: x-axis label
- `xscale`: x-axis scale (default: :log)
- `axis_tight`: Tight axis (default: false)
"""
function bode_plot(freq, y; lw = 1., colorline = :blue, xlab = "Frequency (Hz)", xscale = identity, axis_tight = true, isdeg = false, layout = :vertical, ref_dB = 1.)
    ymag = abs.(y)
    ylab1 = "Magnitude (dB)"
    if isdeg
        ϕ = rad2deg.(angle.(y))
        ylab2 = "Phase (deg)"
    else
        ϕ = angle.(y)
        ylab2 = "Phase (rad)"
    end

    fig = Figure()
    if layout == :vertical
        ax1 = Axis(fig[1,1], xscale = xscale)
        ax2 = Axis(fig[2,1], xscale = xscale)
        ax1.ylabel = ylab1
        ax2.xlabel = xlab
        ax2.ylabel = ylab2
    else
        fig = Figure()
        ax1 = Axis(fig[1,1], xscale = xscale)
        ax2 = Axis(fig[1,2], xscale = xscale)
        ax1.xlabel = xlab
        ax1.ylabel = ylab1
        ax2.xlabel = xlab
        ax2.ylabel = ylab2
    end
    lines!(ax1, freq, 20log10.(ymag/ref_dB), color = colorline, linewidth = lw)
    lines!(ax2, freq, unwrap(ϕ), color = colorline, linewidth = lw)

    if axis_tight
        xlims!(ax1, minimum(freq), maximum(freq))
        xlims!(ax2, minimum(freq), maximum(freq))
    end

    return fig
end