using DSP, FFTW, Parameters, Interpolations, LinearAlgebra, Statistics
@usingany GLMakie

## Include files
includet("../src/utils/undefs.jl")
includet("../src/utils/windows.jl")
includet("../src/estimation/gradient.jl")
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
block_size = 512
fft_params = FFTParameters(sample_rate, block_size)

## Signal reference FRF
freq = fft_params.freq
prob_frf = SdofFRFProblem(sdof, freq)
H = solve(prob_frf).u

## Signal generation
F₀ = 1.
t = fft_params.t
chirp = SweptSine(F₀, t[1], 0.8t[end], freq[1], freq[end], zero_end = true)
x = excitation(chirp, t)
prob = SdofForcedTimeProblem(sdof, [0., 0.], t, x)
y = solve(prob).u

## FRF estimation
overlap_ratio = 0.
win = planck(block_size)
H1 = tfestimate(x, y, fft_params, win, type_frf = :h1, overlap_ratio = overlap_ratio)[1]