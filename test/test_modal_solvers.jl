using StructuralVibration, LinearAlgebra
@usingany GLMakie

## Matrices
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05
t = 0.:1e-2:30.

## Free response
u0 = ([0.2, 0.1], zeros(2))
prob = FreeModalTimeProblem(K, M, ξ, u0, t)
x_free = solve(prob).u

ω, Φ = eigenmode(K, M)
u0_modal = (Φ'*M*[0.2, 0.1], zeros(2))
prob_modal = FreeModalTimeProblem(ω, Φ, ξ, u0_modal, t, ismodal = true)
x_free_modal = solve(prob_modal).u

fig = Figure()
ax_1 = Axis(fig[1, 1], ylabel = "Displacement (m)", title = "Free response")
ax_2 = Axis(fig[2, 1], xlabel = "Time (s)", ylabel = "Displacement (m)")
lines!(ax_1, t, x_free[1, :], label = "x₁ - prob")
lines!(ax_1, t, x_free_modal[1, :], label = "x₁ - prob_modal", linestyle = :dash)
axislegend(ax_1, position = :rt,
           backgroundcolor = (:white, 0.5))
xlims!(ax_1, minimum(t), maximum(t))

lines!(ax_2, t, x_free[2, :], label = "x₂ - prob")
lines!(ax_2, t, x_free_modal[2, :], label = "x₂ - prob_modal", linestyle = :dash)
axislegend(ax_2, position = :rb,
           backgroundcolor = (:white, 0.5))
xlims!(ax_2, minimum(t), maximum(t))

## Harmonic response
F = [1., 2.]
freq = 0.5
u0 = ([0., 1e-4], zeros(2))
prob_harmo = HarmonicModalTimeProblem(K, M, ξ, F, 2π*freq, u0, t)
x_harmo = solve(prob_harmo).u

ω, Φ = eigenmode(K, M)
u0_modal = (Φ'*M*[0., 1e-4], zeros(2))
Fmodal = Φ'*F
prob_harmo_modal = HarmonicModalTimeProblem(ω, Φ, ξ, Fmodal, 2π*freq, u0_modal, t, ismodal = true)
x_harmo_modal = solve(prob_harmo_modal).u

fig_harmo = Figure()
ax_harmo_1 = Axis(fig_harmo[1, 1], ylabel = "Displacement (m)", title = "Free response")
ax_harmo_2 = Axis(fig_harmo[2, 1], xlabel = "Time (s)", ylabel = "Displacement (m)")
lines!(ax_harmo_1, t, x_harmo[1, :], label = "x₁ - prob")
lines!(ax_harmo_1, t, x_harmo_modal[1, :], label = "x₁ - prob_modal", linestyle = :dash)
axislegend(ax_harmo_1, position = :rt,
           backgroundcolor = (:white, 0.5))
xlims!(ax_harmo_1, minimum(t), maximum(t))

lines!(ax_harmo_2, t, x_harmo[2, :], label = "x₂ - prob")
lines!(ax_harmo_2, t, x_harmo_modal[2, :], label = "x₂ - prob_modal", linestyle = :dash)
axislegend(ax_harmo_2, position = :rb,
           backgroundcolor = (:white, 0.5))
xlims!(ax_harmo_2, minimum(t), maximum(t))

## Forced response
u0 = (zeros(2), zeros(2))

F0 = 10.
tstart = 2.
duration = 5.
haversine = HaverSine(F0, tstart, duration)
F0 = excitation(haversine, t)
F = zeros(2, length(t))
F[1, :] .= F0
prob_forced = ForcedModalTimeProblem(K, M, ξ, F, u0, t)

ω, Φ = eigenmode(K, M)
u0_modal = (zeros(2), zeros(2))
Fmodal = Φ'*F
prob_forced_modal = ForcedModalTimeProblem(ω, Φ, ξ, Fmodal, u0_modal, t, ismodal = true)

x_forced = solve(prob_forced).u
x_forced_modal = solve(prob_forced_modal).u

fig_forced = Figure()
ax_forced_1 = Axis(fig_forced[1, 1], ylabel = "Displacement (m)", title = "Forced response")
ax_forced_2 = Axis(fig_forced[2, 1], xlabel = "Time (s)", ylabel = "Displacement (m)")
lines!(ax_forced_1, t, x_forced[1, :], label = "x₁ - prob")
lines!(ax_forced_1, t, x_forced_modal[1, :], label = "x₁ - prob_modal", linestyle = :dash)
axislegend(ax_forced_1, position = :rt,
           backgroundcolor = (:white, 0.5))
xlims!(ax_forced_1, minimum(t), maximum(t))

lines!(ax_forced_2, t, x_forced[2, :], label = "x₂ - prob")
lines!(ax_forced_2, t, x_forced_modal[2, :], label = "x₂ - prob_modal", linestyle = :dash)
axislegend(ax_forced_2, position = :rt,
           backgroundcolor = (:white, 0.5))
xlims!(ax_forced_2, minimum(t), maximum(t))

fig_forced

## Impulse response
h = impulse_response(K, M, ξ, t, ismat = true)

fig_h = Figure()
ax_h11 = Axis(fig_h[1, 1], ylabel = "Impulse response")
ax_h12 = Axis(fig_h[1, 2])
ax_h21 = Axis(fig_h[2, 1], xlabel = "Time (s)", ylabel = "Impulse response")
ax_h22 = Axis(fig_h[2, 2], xlabel = "Time (s)")
lines!(ax_h11, t, h[1, 1, :], label = "h₁₁")
xlims!(ax_h11, minimum(t), maximum(t))
axislegend(ax_h11, position = :rt,
           backgroundcolor = (:white, 0.5))

lines!(ax_h12, t, h[1, 2, :], label = "h₁₂")
xlims!(ax_h12, minimum(t), maximum(t))
axislegend(ax_h12, position = :rt,
           backgroundcolor = (:white, 0.5))

lines!(ax_h21, t, h[2, 1, :], label = "h₂₁")
xlims!(ax_h21, minimum(t), maximum(t))
axislegend(ax_h21, position = :rt,
           backgroundcolor = (:white, 0.5))

lines!(ax_h22, t, h[2, 2, :], label = "h₂₂")
xlims!(ax_h22, minimum(t), maximum(t))
axislegend(ax_h22, position = :rt,
           backgroundcolor = (:white, 0.5))

fig_h