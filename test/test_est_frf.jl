using DSP, FFTW, Parameters, LinearAlgebra, Statistics
@usingany GLMakie

## Include files
includet("../src/utils/utils.jl")
includet("../src/models/sdof_mdof.jl")
includet("../src/solvers/sdof_solvers.jl")
includet("../src/models/excitation.jl")
includet("../src/estimation/signal_processing.jl")

## Sdof structure
m = 1.
f₀ = 10.
ξ = 0.01
sdof = Sdof(m, f₀, ξ)

## FFT Parameters
sample_rate = 128
block_size = 256
fft_params = FFTParameters(sample_rate, block_size)

## Signal reference FRF
t = fft_params.t
freq = fft_params.freq
prob = SdofFRFProblem(sdof, freq)
H = solve(prob).u

## Signal acquisition
F₀ = 1.
chirp = SweptSine(F₀, t[1], t[end], freq[1], freq[end], zero_end = true)
F = excitation(chirp, t)