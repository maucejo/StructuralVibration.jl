#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
#| output: false
using StructuralVibration, ShareAdd
using DSP: rect, hanning, hamming, tukey, cosine, lanczzos, triang, bartlett, gaussian, bartlett_hann, blackman, kaiser, dpss
@usingany CairoMakie
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| output: false
rectwin = rect(1024)
#
#
#
#| echo: false
fig_rect = Figure()
ax_rect = Axis(fig_rect[1, 1], title = "Rectangular window", xlabel = "Samples", ylabel = "Window function")
lines!(ax_rect, 1:1024, rectwin)

display(fig_rect);
#
#
#
#
#
