using Parameters, DSP, FFTW, LinearAlgebra, Interpolations
@usingany GLMakie

includet("../src/models/sdof.jl")
includet("../src/solvers/sdof_solvers.jl")
includet("../src/utils/excitation.jl")


## SDOF system
m = 1.
ω₀ = 2π*10.
ξ = 0.01
sdof = Sdof(m, ω₀, ξ)

## Excitation

#Time vector
Δt = 1e-3
t = 0.:Δt:10.

# Excitation
F₀ = 10.
rect = Rectangle(F₀, t[1], t[end])
F = excitation(rect, t)

## Check Duhamel's integral

# Exact solution
Ω₀ = ω₀*√(1 - ξ^2)
xexact = @. F₀*(Ω₀ - (Ω₀*cos(Ω₀*t) + ξ*ω₀*sin(Ω₀*t))*exp(-ξ*ω₀*t))/m/Ω₀/(Ω₀^2 + ξ^2*ω₀^2)

# Duhamel's integral
prob = SdofTimeProblem(sdof, F = F)
x = solve(prob, [0., 0.], t, rect).x

lines(t, xexact, color = :blue)
lines!(t, x, color = :red, linestyle = :dash)

## Frequency response

# Response calculation
freq = 1.:0.01:30.
prob_resp = SdofFrequencyProblem(sdof, freq, 1e-2ones(length(freq)))
y = solve(prob_resp).y
lines(freq, 20log10.(abs.(y)), color = :blue)

## FRF
prob_frf = SdofFrequencyProblem(sdof)
H = solve(prob_frf, freq).y
lines(freq, 20log10.(abs.(H)), color = :blue)
lines(freq, angle.(H), color = :blue)