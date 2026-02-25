using StructuralVibration, FFTW
import DSP

## Windows from DSP.jl
rectwin = rect(1024)
hannwin = hanning(1024)
hammwin = hamming(1024)
tukeywin = tukey(1024, 0.5)
coswin = cosine(1024)
lanczoswin = lanczos(1024)
triangwin = triang(1024)
bartlettwin = bartlett(1024)
barthannwin = bartlett_hann(1024)
gaussianwin = gaussian(1024, 0.2)
blackmanwin = blackman(1024)
kaiserwin = kaiser(1024, 0.5)
dpsswin = dpss(1024, 4, 1)

## Windows from StructuralVibration.jl
exponentialwin = exponential(1024, 0.1)
force_expwin = force(1024, 0.2, 0.01)
flattopwin = flattop(1024)
nuttallwin = nuttall(1024)
blackman_nuttallwin = blackman_nuttall(1024)
blackman_harriswin = blackman_harris(1024)
parzenwin = parzen(1024)
planckwin = planck(1024, 0.1)
flattriwin = flattri(1024, 0.25)



## Anti-aliasing filter
# Time parameters
Δt = 1e-4
t = 0.:Δt:10.
fs = 1/Δt

# Signal
nt = length(t)
y = randn(nt)

# filtered signal
yf = anti_alias(y, 2500., fs = fs)


## Response spectrum estimation
# Structural parameters
m = 1.
f0 = 25.
ξ = 0.1

sdof = Sdof(m, f0, ξ)

# Acquisition parameters
sample_rate = 256
block_size = 1024
fft_params = FFTParameters(sample_rate, block_size)
freq = fft_params.freq
t = fft_params.t

# Excitation signal generation
F0 = 10.
chirp = SweptSine(F0, t[1], 0.8t[end], freq[1], freq[end], zero_end = true)
x = excitation(chirp, t)

# Reference signal
prob = SdofForcedTimeProblem(sdof, [0., 0.], t, x)
y = solve(prob).u

# Frequency domain
Fx = Fx = rfft(x)[1:length(freq)]/block_size
Fx[2:end] .*= 2.

# Response spectrum - Reference
prob_y = SdofFrequencyProblem(sdof, freq, Fx)
u = solve(prob_y).u

# Spectrum stimation
uest = spectrum(y, block_size, rect, fs = sample_rate)[1]



## Autopower & Cross-power spectral density
# Acquisition parameters
fs = 1000
bs = 100

fft_params = FFTParameters(fs, bs)
freq = fft_params.freq
t = (1:fs)/fs

# Signal generation
f = [100 150]                      # 100Hz & 150Hz frequencies
A = [1; 2]                         # Amplitudes
x = sin.(2π*f.*t)*A + randn(1000)

# PSD estimation from DSP.jl - For comparison
psd = DSP.welch_pgram(x, bs; fs = fs, window = hamming)
pxx_ref = DSP.power(psd)

# PSD estimation from StructuralVibration.jl
pxx = welch(x, bs, hamming, fs = fs)[1]
pxx_csd = csd(x, x, bs, hamming, fs = fs)[1]


## Frequency response function estimation
# Structural parameters
m = 1.
f0 = 25.
ξ = 0.1

sdof = Sdof(m, f0, ξ)

# Acquisition parameters for one block
sample_rate = 256
block_size = 1024
fft_params = FFTParameters(sample_rate, block_size)

# Reference FRF
freq = fft_params.freq
prob_frf = SdofFRFProblem(sdof, freq)
H = solve(prob_frf).u

# Signal generation - Input signal
nblocks = 5
tb = fft_params.t
dt = fft_params.dt
t = tb[1]:dt:(nblocks*(tb[end] + dt) - dt)

F0 = 10.
chirp = SweptSine(F0, tb[1], 0.8tb[end], freq[1], freq[end], zero_end = true)
x = repeat(excitation(chirp, tb), outer = nblocks)

# Signal generation - Output signal
prob = SdofForcedTimeProblem(sdof, [0., 0.], t, x)
y = solve(prob).u

# FRF estimation
win(x) = tukey(x, 0.25)
H1 = tfestimate(x, y, block_size, win, fs = sample_rate)[1]