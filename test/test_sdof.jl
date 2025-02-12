using StructuralVibration
@usingany GLMakie

## SDOF system
m = 1.
f₀ = 10.
ξ = 0.01
sdof = Sdof(m, f₀, ξ)
(; ω₀) = sdof

## Excitation

#Time vector
Δt = 1e-2
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
prob = SdofForcedTimeProblem(sdof, [0., 0.], t, F)
x = solve(prob).u

sv_plot(t, [xexact x]', lw = 2., legend = (active = true, entry = ["Exact", "Duhamel"]))

## Frequency response

# Response calculation
freq = 1.:0.01:30.
prob_resp = SdofFrequencyProblem(sdof, freq, 1e-2ones(length(freq)))
y = solve(prob_resp).u
lines(freq, 20log10.(abs.(y)), color = :blue)

## FRF
prob_frf = SdofFRFProblem(sdof, freq)
H = solve(prob_frf).u
bode_plot(freq, H)