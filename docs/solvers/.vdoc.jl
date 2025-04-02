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
#| output: false
using StructuralVibration, ShareAdd
@usingany CairoMakie
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
@doc SdofFreeTimeProblem
```
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
@doc solve(prob::SdofFreeTimeProblem)
```
#
#
#
#
#
#| ouput: false

# Structural parameters
m = 1.
f0 = 1.

# Time vector
t = 0.:0.01:10.

# Initial conditions
u0 = [1, -2]

# Undamped system
sdof_nd = Sdof(m, f0, 0.)
prob_nd = SdofFreeTimeProblem(sdof_nd, u0, t)
x_nd = solve(prob_nd).u

# Underdamped system
sdof_ud = Sdof(m, f0, 0.1)
prob_ud = SdofFreeTimeProblem(sdof_ud, u0, t)
x_ud = solve(prob_ud).u

# Critically damped system
sdof_cd = Sdof(m, f0, 1.)
prob_cd = SdofFreeTimeProblem(sdof_cd, u0, t)
x_cd = solve(prob_cd).u

# Overdamped system
sdof_od = Sdof(m, f0, 1.2)
prob_od = SdofFreeTimeProblem(sdof_od, u0, t)
x_od = solve(prob_od).u;
#
#
#
#| echo: false

xmin, xmax = extrema([x_nd; x_ud; x_cd; x_od])

fig = Figure()
ax = Axis(fig[1, 1], ylabel="Displacement (m)", title = "Undamped system")
lines!(ax, t, x_nd)
xlims!(ax, 0., 10.)
ylims!(ax, xmin, xmax)

ax1 = Axis(fig[1, 2], title = "Underdamped system")
lines!(ax1, t, x_ud)
xlims!(ax1, 0., 10.)
ylims!(ax1, xmin, xmax)

ax2 = Axis(fig[2, 1], xlabel="Time (s)", ylabel="Displacement (m)", title = "Critically damped system")
lines!(ax2, t, x_cd)
xlims!(ax2, 0., 10.)
ylims!(ax2, xmin, xmax)

ax3 = Axis(fig[2, 2], xlabel="Time (s)", title = "Overdamped system")
lines!(ax3, t, x_od)
xlims!(ax3, 0., 10.)
ylims!(ax3, xmin, xmax)

display(fig);
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
@doc SdofHarmonicTimeProblem
```
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
@doc solve(prob::SdofHarmonicTimeProblem)
```
#
#
#
#
#
#| output: false

# Structural parameters
m = 1.
f0 = 1.
ξ = 0.1

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

# Harmonic excitation - Force amplitude + excitation frequency
F0 = 10.
f = 2.

# Instantiation of the problem
t = 0.:0.01:10.
u0 = [1, -2]
prob = SdofHarmonicTimeProblem(sdof, F0, f, u0, t, :force)

# Solve the problem
x_harmo = solve(prob).u
#
#
#
#| echo: false

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Displacement (m)", title = "Harmonic excitation")
lines!(ax, t, x_harmo)
xlims!(ax, 0., 10.)

display(fig);
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
@doc SdofForcedTimeProblem
```
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
@doc solve(prob::SdofForcedTimeProblem)
```
#
#
#
#
#
# Structural parameters
m = 1.
f0 = 1.
ξ = 0.1

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

# Haversine excitation signal
F0 = 10.
tstart = 0.5
duration = 2.
haversine = HaverSine(F0, tstart, duration)

t = 0.:0.01:10.
F = excitation(haversine, t)

# Instantiation of the problem
u0 = [1, -2]
prob = SdofForcedTimeProblem(sdof, F, u0, t, :force)

# Solve the problem
x_arb_filt = solve(prob, method = :filt).u
x_arb_interp = solve(prob, method = :interp).u
x_arb_conv = solve(prob, method = :conv).u;
#
#
#
#| echo: false

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Displacement (m)", title = "Arbitrary excitation - Haversine")
lines!(ax, t, x_arb_filt, label = "Filtering")
lines!(ax, t, x_arb_interp, linestyle = :dash, label = "Interp. + Quad.")
lines!(ax, t, x_arb_conv, linestyle = :dashdot , label = "Convolution")
axislegend(ax, position = :rt,
           backgroundcolor = (:white, 0.5))
xlims!(ax, 0., 10.)

display(fig);
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
@doc impulse_response(sdof, t)
```
#
#
#
#| output: false

# Structural parameters
m = 1.
f0 = 1.
ξ = 0.1

t = 0.:0.01:10.

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

h = impulse_response(sdof, t)
#
#
#
#| echo: false

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Impulse response", title = "Impulse response")
lines!(ax, t, h)
xlims!(ax, 0., 10.)

display(fig);
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
@doc srs
```
#
#
#
#| output: false

# Definition of the base acceleration
t = 0.:1e-5:0.5
base_acc = HalfSine(10., 0., 1e-3)

# Compute the SRS
f = 50.:25:1e3
srs_basic = srs(base_acc, f, t)
srs_rec_int = srs(base_acc, f, t, alg = :RecursiveInt)
srs_rec_filt = srs(base_acc, f, t, alg = :RecursiveFilt)
srs_smallwood = srs(base_acc, f, t, alg = :Smallwood)
#
#
#
#| echo: false

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Frequency (Hz)", ylabel="SRS", title = "Shock response spectrum")
lines!(ax, f, srs_basic, label = "Basic")
lines!(ax, f, srs_rec_filt, label = "Recursive filtering")
lines!(ax, f, srs_rec_int, label = "Recursive integration")
lines!(ax, f, srs_smallwood, label = "Smallwood")
axislegend(ax, position = :rb, backgroundcolor = (:white, 0.5),)
xlims!(ax, 50, 1000)

display(fig);
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
@doc SdofFRFProblem
```
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
@doc solve(prob::SdofFRFProblem)
```
#
#
#
#
#
#| output: false

# Structural parameters
m = 1.
f0 = 10.
ξ = 0.01

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

# Definition of the frequency range
f = 1.:0.1:30.

# Instantiation of the problem - Force excitation
prob = SdofFRFProblem(sdof, f)

# Solve the problem
H = solve(prob).u
#
#
#
#| echo: false

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Frequency (Hz)", ylabel = "FRF (m/N)", title = "Frequency Response Function")
lines!(ax, f, abs.(H))
xlims!(ax, 1., 30.)

display(fig);
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
@doc SdofFrequencyProblem
```
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
@doc solve(prob::SdofFRFProblem)
```
#
#
#
#
#
#| output: false

# Structural parameters
m = 1.
f0 = 10.
ξ = 0.01

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

# Definition of the frequency range
f = 1.:0.1:30.

# Instantiation of the problem - Force excitation
F = fill(10., length(f))
prob = SdofFrequencyProblem(sdof, f, F)

# Solve the problem
y = solve(prob).u
#
#
#
#| echo: false

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Frequency (Hz)", ylabel = "Displacement (m)", title = "Response spectrum")
lines!(ax, f, abs.(y))
xlims!(ax, 1., 30.)

display(fig);
#
#
#
