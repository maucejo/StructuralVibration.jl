using Parameters, DSP, LinearAlgebra, Interpolations
@usingany GLMakie

includet("../src/models/sdof.jl")
includet("../src/utils/excitation.jl")

## SDOF system
m = 1.
ω₀ = 2π*10.
ξ = 0.01
sdof = SDOF(m, ω₀, ξ)

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
x = forced_response(sdof, F, t)

lines(t, xexact, color = :blue)
lines!(t, x, color = :red, linestyle = :dash)

## Frequency response

# Response calculation
freq = 1.:0.01:30.
y = freq_response(sdof, 2ones(length(freq)), freq)
lines(freq, 20log10.(abs.(y)), color = :blue)

## FRF
H = frf(sdof, freq)
lines(freq, 20log10.(abs.(H)), color = :blue)
lines(freq, angle.(H), color = :blue)