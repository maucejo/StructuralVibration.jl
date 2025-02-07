using StructuralVibration, DSP
@usingany GLMakie, MAT

vars = matread("test/matlab_comparison/signal_processing/anti_aliasing/test_anti_aliasing.mat")

hf = vec(vars["hf"])

fs = 4096
fn = fs/2

# Filter design
order = 200     # Filter order
fb = 0.05*fn    # Filter bandwith = 2*fb
freq_filt = [(0., fn - fb) => 1]
filt_coeff = remez(order, freq_filt, Hz = fs, maxiter = 50)

# Filter visualization
z =  PolynomialRatio(filt_coeff, [1.0])
h, w = freqresp(z)
freq = w*fn/π

y = [transpose(hf); transpose(h)]

bode_plot(freq, y)

# Filtering
x = randn(fs)

y = filtfilt(filt_coeff, x)