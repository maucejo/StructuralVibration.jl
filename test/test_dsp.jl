using StructuralVibration, FFTW
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

## Signal reference FRF
freq = fft_params.freq
prob_frf = SdofFRFProblem(sdof, freq)
H = solve(prob_frf).u

# Signal generation
F0 = 10.
t = fft_params.t
chirp = SweptSine(F0, t[1], 0.8t[end], freq[1], freq[end], zero_end = true)
x = excitation(chirp, t)
prob = SdofForcedTimeProblem(sdof, x, [0., 0.], t)
y = solve(prob).u

# FRF estimation
win = planck(block_size)
H1 = tfestimate(x, y, block_size, win, fs = sample_rate, overlap = 0.5)[1]

fig_h = Figure()
ax_h = Axis(fig_h[1, 1], xlabel = "Frequency (Hz)", ylabel = "Amplitude (dB)")
lines!(ax_h, freq, 20log10.(abs.(H)), label = "Reference")
lines!(ax_h, freq, 20log10.(abs.(H1)), linestyle = :dash, label = "Estimated")
axislegend(ax_h, position = :rt, backgroundcolor = (:white, 0.5),)
xlims!(ax_h, 0., 50.)
ylims!(ax_h, -100., -70.)

## Spectrum
Fx = rfft(x)[1:length(freq)]/block_size
Fx[2:end] .*= 2.
# Fx = spectrum(x, block_size, rect(block_size), fs = sample_rate, overlap = 0.)[1]
prob_y = SdofFrequencyProblem(sdof, Fx, freq)
u = solve(prob_y).u

# Estimation
uest = spectrum(y, block_size, rect(block_size), fs = sample_rate, overlap = 0.95)[1]

fig_u = Figure()
ax_u = Axis(fig_u[1, 1], xlabel = "Frequency (Hz)", ylabel = "Amplitude (dB)")
lines!(ax_u, freq, 20log10.(abs.(u)), label = "Reference")
lines!(ax_u, freq, 20log10.(abs.(uest)), linestyle = :dash, label = "Estimated")
axislegend(ax_u, position = :rt, backgroundcolor = (:white, 0.5),)
xlims!(ax_u, 0., 50.)
ylims!(ax_u, -110., -80.)

## Power spectral density
Fs = 1000                         # Sampling frequency
bs = 100
freq = rfftfreq(bs, Fs)
t = (1:Fs)/Fs                     # 1000 time samples
f = [100 150]                     # 100Hz & 150Hz frequencies
A = [1; 2]                        # Amplitudes
x = sin.(2π*f.*t)*A + randn(1000); # Noise corrupted x signal
psd = DSP.welch_pgram(x, bs; fs=Fs, window = hamming)
pxx_ref = power(psd)

pxx = welch(x, bs, hamming(bs), fs = Fs)[1]

fig_pxx = Figure()
ax_pxx = Axis(fig_pxx[1, 1], xlabel = "Frequency (Hz)", ylabel = "Amplitude (dB)")
lines!(ax_pxx, freq, 10log10.(pxx_ref), label = "DSP.jl")
lines!(ax_pxx, freq[1:length(pxx)], 10log10.(pxx), linestyle = :dash, label = "StructuralVibration.jl")
axislegend(ax_pxx, position = :rt, backgroundcolor = (:white, 0.5),)
xlims!(ax_pxx, 0., freq[length(pxx)])