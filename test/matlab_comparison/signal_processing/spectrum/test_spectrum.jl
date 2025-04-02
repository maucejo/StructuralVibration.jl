using StructuralVibration, DSP
@usingany GLMakie, MAT

# Load Matlab data for comparison
vars = matread("test/matlab_comparison/signal_processing/spectrum/test_spectrum.mat")
pxx_mat = vec(vars["pxx_mat"])
fmat = vec(vars["fmat"])
signal = vec(vars["signal"])

fs = 1024
bs = 512

spec, f = spectrum(signal, bs, hanning(bs), fs = fs, overlap = 0)

N = length(f)
y = [10log10.(pxx_mat[1:N])'; 20log10.(abs.(spec))']

lg = (active = true, entry = ["Matlab", "Julia"])
sv_plot(f, y, xlabel = "Frequency (Hz)", ylabel = "Amplitude (dB)", legend = lg)