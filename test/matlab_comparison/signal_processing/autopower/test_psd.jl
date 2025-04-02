using StructuralVibration, DSP
@usingany GLMakie, MAT

# Load Matlab data for comparison
vars = matread("test/matlab_comparison/signal_processing/autopower/test_psd.mat")
pxx_mat = vec(vars["pxx_mat"])
fmat = vec(vars["fmat"])
signal = vec(vars["signal"])

fs = 1024
bs = 512
pxx, f = welch(signal, bs, hamming(512), fs = fs, overlap = 0.5)

N = length(f)
y = [10log10.(pxx_mat[1:N])'; 10log10.(pxx)']

lg = (active = true, entry = ["Matlab", "Julia"])
sv_plot(f, y, xlabel = "Frequency (Hz)", ylabel = "PSD (dB)", legend = lg)