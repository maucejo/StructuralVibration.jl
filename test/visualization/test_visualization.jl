using StructuralVibration
@usingany CairoMakie

# Initialize a Sdof type
m = 1.
f0 = 10.
ξ = 0.01
sdof = Sdof(m, f0, ξ)

# Computation parameters
freq = 1.:0.01:30.

# Compute the FRF
prob_frf = SdofFRFProblem(sdof, freq)
H = solve(prob_frf).u

# Bode plot
bode_plot(freq, H)

# Nyquist plot
nyquist_plot(H)

# Nyquist plot - 3D
nyquist_plot(freq, H, projection = true)



## Waterfall plot
x = range(0., 2π, 100)
y = range(0., 1., 5)

nx = length(x)
ny = length(y)
z = zeros(ny, nx)

for i in eachindex(y)
    z[i, :] = sin.(i*x/2.)
end

waterfall_plot(x, y, z, xlim = [-0.1, 2π + 0.1], ylim = [-0.1, 1.1])



## SV plot
x = range(0., 2π, 100)
z_sv = ntuple(i -> sin.(i*x/2), 5)

# SV plot
sv_plot(x, z_sv..., lw = 2., legend = (active = true,))


## Theming

# Makie theme
with_theme(theme_choice(:makie)) do
    sv_plot(x, z_sv..., lw = 2., legend = (active = true,), title = ":makie theme")
end

# SV theme
with_theme(theme_choice(:makie)) do
    sv_plot(x, z_sv..., lw = 2., legend = (active = true,), title = ":sv theme")
end