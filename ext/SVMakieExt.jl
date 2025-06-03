module SVMakieExt

using StructuralVibration, Makie, DSP

import StructuralVibration: sv_plot, bode_plot, nyquist_plot, waterfall_plot, theme_choice

include("./common_plots.jl")
end