using StructuralVibration, FFTW, DSP
import DSP
@usingany GLMakie

# Structural parameters
m = 1.
f0 = 25.
ξ = 0.1

sdof = Sdof(m, f0, ξ)

# Acquisition parameters
sample_rate = 256
block_size = 1024
fft_params = FFTParameters(sample_rate, block_size)

# Reference FRF
freq = fft_params.freq
prob_frf = SdofFRFProblem(sdof, freq)
H = solve(prob_frf).u

# Signal generation
F0 = 10.
nblocks = 5
tb = fft_params.t
dt = fft_params.dt
t = tb[1]:dt:(nblocks*(tb[end] + dt) - dt)
chirp = SweptSine(F0, tb[1], 0.8tb[end], freq[1], freq[end], zero_end = true)

x = repeat(excitation(chirp, tb), outer = nblocks)
# x = excitation(chirp, t)
prob = SdofForcedTimeProblem(sdof, x, [0., 0.], t)
y = solve(prob).u

# FRF estimation
win = tukey(block_size, 0.25)
H1 = tfestimate(x, y, block_size, win, fs = sample_rate, overlap = 0)[1]

fig_h = Figure()
ax_h = Axis(fig_h[1, 1], xlabel = "Frequency (Hz)", ylabel = "Amplitude (dB)")
lines!(ax_h, freq, 20log10.(abs.(H)), label = "Reference")
lines!(ax_h, freq, 20log10.(abs.(H1)), linestyle = :dash, label = "H1 estimator")
axislegend(ax_h, position = :rt, backgroundcolor = (:white, 0.5),)
xlims!(ax_h, 0., 50.)
ylims!(ax_h, -100., -70.)

## Spectrum
xs = excitation(chirp, tb)
Fx = rfft(xs)[1:length(freq)]/block_size
Fx[2:end] .*= 2.
prob_y = SdofFrequencyProblem(sdof, Fx, freq)
u = solve(prob_y).u

# Estimation
prob_spec = SdofForcedTimeProblem(sdof, xs, [0., 0.], tb)
ys = solve(prob_spec).u
uest = spectrum(ys, block_size, rect(block_size), fs = sample_rate)[1]

fig_u = Figure()
ax_u = Axis(fig_u[1, 1], xlabel = "Frequency (Hz)", ylabel = "Amplitude (dB)")
lines!(ax_u, freq, 20log10.(abs.(u)), label = "Reference")
lines!(ax_u, freq, 20log10.(abs.(uest)), linestyle = :dash, label = "Estimated")
axislegend(ax_u, position = :rt, backgroundcolor = (:white, 0.5),)
xlims!(ax_u, 0., 50.)
ylims!(ax_u, -110., -80.)

## Power spectral density
fs = 1000                          # Sampling frequency
bs = 100
# freq = rfftfreq(bs, fs)
                   # 1000 time samples

fft_params = FFTParameters(fs, bs, pow2 = false)
freq = fft_params.freq
t = (1:fs)/fs
# t = fft_params.t

f = [100 150]                      # 100Hz & 150Hz frequencies
A = [1; 2]                         # Amplitudes
x = sin.(2π*f.*t)*A + randn(1000); # Noise corrupted x signal
psd = DSP.welch_pgram(x, bs; fs = fs, window = hamming)
pxx_ref = DSP.power(psd)

pxx = welch(x, bs, hamming(bs), fs = fs)[1]

fig_pxx = Figure()
ax_pxx = Axis(fig_pxx[1, 1], xlabel = "Frequency (Hz)", ylabel = "Amplitude (dB)")
lines!(ax_pxx, freq, 10log10.(pxx_ref)[1:length(pxx)], label = "DSP.jl")
lines!(ax_pxx, freq, 10log10.(pxx), linestyle = :dash, label = "StructuralVibration.jl")
axislegend(ax_pxx, position = :rt, backgroundcolor = (:white, 0.5),)
xlims!(ax_pxx, 0., freq[end])

## Anti-aliasing
# Time parameters
Δt = 1e-4
t = 0.:Δt:10.
fs = 1/Δt

# Signal
nt = length(t)
y = randn(nt)

# filtered signal
yf = anti_alias(y, 2500., fs = fs)

freq = rfftfreq(nt, fs)

fig_a = Figure()
ax_a1 = Axis(fig_a[1, 1], xlabel = "Time (s)", ylabel = "Signal", title = "Original data")
ax_a2 = Axis(fig_a[2, 1], xlabel = "Frequency (Hz)", ylabel = "Amplitude (dB)")
lines!(ax_a1, t, y, label = "Original signal")
xlims!(ax_a1, 0., t[end])

lines!(ax_a2, freq, 10log10.(abs.(rfft(y)[1:length(freq)])), label = "FFT")
xlims!(ax_a2, 0., freq[end])

fig_b = Figure()
ax_b1 = Axis(fig_b[1, 1], xlabel = "Time (s)", ylabel = "Signal", title = "Filtered data")
ax_b2 = Axis(fig_b[2, 1], xlabel = "Frequency (Hz)", ylabel = "Amplitude (dB)")
lines!(ax_b1, t, yf, label = "Filtered signal")
xlims!(ax_b1, 0., t[end])

lines!(ax_b2, freq, 10log10.(abs.(rfft(yf)[1:length(freq)])), label = "FFT")
xlims!(ax_b2, 0., freq[end])