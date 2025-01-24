module CairoMakieExt

using StructuralVibration, CairoMakie, DSP

import StructuralVibration: plot, bode_plot, nyquist_plot, waterfall_plot

include("./common_plots.jl")
end