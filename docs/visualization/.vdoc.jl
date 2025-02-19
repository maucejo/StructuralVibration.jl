#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
@doc bode_plot
```
#
#
#
# Initialize a Sdof type
m = 1.
f₀ = 10.
ξ = 0.01
sdof = Sdof(m, f₀, ξ)

# Computation parameters
freq = 1.:0.01:30.

# Compute the FRF
prob_frf = SdofFRFProblem(sdof, freq)
H = solve(prob_frf).u
bode_plot(freq, H)
#
#
#
