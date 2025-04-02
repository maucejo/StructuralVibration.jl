using StructuralVibration, DSP
@usingany GLMakie, MAT

# Load Matlab data for comparison
vars = matread("test/matlab_comparison/signal_processing/tf/test_tf.mat")
Hmat = vec(vars["H"])
fmat = vec(vars["f"])
signal1 = vec(vars["signal1"])
signal2 = vec(vars["signal2"])

fs = 1024
bs = 512
H, f = tfestimation(signal1, signal2, bs, hamming(bs), fs = fs, overlap = 0.)[1:2]

N = length(f)
y = [20log10.(abs.(Hmat[1:N]))'; 20log10.(abs.(H))']

lg = (active = true, entry = ["Matlab", "Julia"])
sv_plot(f, y, xlabel = "Frequency (Hz)", ylabel = "Magnitude (dB)", legend = lg)