using StructuralVibration
@usingany GLMakie

## Free Response
m = 1.
f0 = 1.
ξ = 0.1

# Time vector
t = 0.:0.01:10.

# Initial conditions
u0 = [1., -2.]

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
x_od = solve(prob_od).u

xmin, xmax = extrema([x_nd; x_ud; x_cd; x_od])

fig = Figure()
ax = Axis(fig[1, 1], ylabel="Displacement (m)", title = "Undamped system")
lines!(ax, t, x_nd)
xlims!(ax, minimum(t), maximum(t))
ylims!(ax, xmin, xmax)

ax1 = Axis(fig[1, 2], title = "Underdamped system")
lines!(ax1, t, x_ud)
xlims!(ax1, minimum(t), maximum(t))
ylims!(ax1, xmin, xmax)

ax2 = Axis(fig[2, 1], xlabel="Time (s)", ylabel="Displacement (m)", title = "Critically damped system")
lines!(ax2, t, x_cd)
xlims!(ax2, minimum(t), maximum(t))
ylims!(ax2, xmin, xmax)

ax3 = Axis(fig[2, 2], xlabel="Time (s)", title = "Overdamped system")
lines!(ax3, t, x_od)
xlims!(ax3, minimum(t), maximum(t))
ylims!(ax3, xmin, xmax)

fig

## Harmonic Response
sdof = Sdof(m, f0, ξ)
F0 = 10.
f = 2.

# Instantiation of the problem
t = 0.:0.01:10.
u0 = [1., -2.]
prob = SdofHarmonicTimeProblem(sdof, F0, f, u0, t, :force)

# Solve the problem
x_harmo = solve(prob).u

fig_harmo = Figure()
ax_harmo = Axis(fig_harmo[1, 1], xlabel="Time (s)", ylabel="Displacement (m)", title = "Harmonic excitation")
lines!(ax_harmo, t, x_harmo)
xlims!(ax_harmo, minimum(t), maximum(t))

fig_harmo

## Arbibtrary excitation
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
x_arb_conv = solve(prob, method = :conv).u

fig_arb = Figure()
ax_arb = Axis(fig_arb[1, 1], xlabel="Time (s)", ylabel="Displacement (m)", title = "Arbitrary excitation - Haversine")
lines!(ax_arb, t, x_arb_filt, label = "Filtering")
lines!(ax_arb, t, x_arb_interp, linestyle = :dash, label = "Interp. + Quad.")
lines!(ax_arb, t, x_arb_conv, linestyle = :dashdot , label = "Convolution")
axislegend(ax_arb, position = :rt,
           backgroundcolor = (:white, 0.5))
xlims!(ax_arb, minimum(t), maximum(t))

fig_arb

## Impulse response
m = 1.
f0 = 1.
ξ = 0.1

t = 0.:0.01:10.

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)
h = impulse_response(sdof, t)

fig_h = Figure()
ax_h = Axis(fig_h[1, 1], xlabel = "Time (s)", ylabel = "Impulse response", title = "Impulse response")
lines!(ax_h, t, h)
xlims!(ax_h, minimum(t), maximum(t))

fig_h

## SRS
t = 0.:1e-5:0.5
base_acc = HalfSine(10., 0., 1e-3)

# Compute the SRS
f = 50.:25:1e3
srs_basic = srs(base_acc, f, t)
srs_rec_int = srs(base_acc, f, t, alg = :RecursiveInt)
srs_rec_filt = srs(base_acc, f, t, alg = :RecursiveFilt)
srs_smallwood = srs(base_acc, f, t, alg = :Smallwood)

fig_srs = Figure()
ax_srs = Axis(fig_srs[1, 1], xlabel="Frequency (Hz)", ylabel="SRS", title = "Shock response spectrum")
lines!(ax_srs, f, srs_basic, label = "Basic")
lines!(ax_srs, f, srs_rec_filt, label = "Recursive filtering")
lines!(ax_srs, f, srs_rec_int, label = "Recursive integration")
lines!(ax_srs, f, srs_smallwood, label = "Smallwood")
axislegend(ax_srs, position = :rb, backgroundcolor = (:white, 0.5),)
xlims!(ax_srs, minimum(f), maximum(f))

fig_srs

## FRF
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

fig_frf = Figure()
ax_frf = Axis(fig_frf[1, 1], xlabel = "Frequency (Hz)", ylabel = "FRF (m/N)", title = "Frequency Response Function")
lines!(ax_frf, f, abs.(H))
xlims!(ax_frf, minimum(f), maximum(f))

fig_frf

## Frequency spectrum
m = 1.
f0 = 10.
ξ = 0.01

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

# Definition of the frequency range
f = 1.:0.1:30.

# Instantiation of the problem - Force excitation
F = fill(10., length(f))
prob = SdofFrequencyProblem(sdof, F, f)

# Solve the problem
y = solve(prob).u

fig_y = Figure()
ax_y = Axis(fig_y[1, 1], xlabel = "Frequency (Hz)", ylabel = "Displacement (m)", title = "Response spectrum")
lines!(ax_y, f, abs.(y))
xlims!(ax_y, minimum(f), maximum(f))