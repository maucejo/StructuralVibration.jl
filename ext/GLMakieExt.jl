module GLMakieExt

using StructuralVibration, GLMakie, DSP

import StructuralVibration: sv_plot, bode_plot, nyquist_plot, waterfall_plot

include("./common_plots.jl")
end