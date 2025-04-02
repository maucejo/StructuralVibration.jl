using StructuralVibration, LinearAlgebra
@usingany GLMakie

## Matrices
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

t = 0.:1e-2:30.

## Free response
F_free = zeros(2, length(t))
u0 = ([0.2, 0.1], zeros(2))
prob_free = DirectTimeProblem(K, M, C, F_free, u0, t)
x_free_gα = solve(prob_free).u
x_free_cd = solve(prob_free, CentralDiff()).u
x_free_rk = solve(prob_free, RK4()).u

prob_free_modal =  FreeModalTimeProblem(K, M, ξ, u0, t)
x_free_modal = solve(prob_free_modal).u

fig = Figure(size = (600, 600))
ax_1 = Axis(fig[1, 1], ylabel = "Displacement (m)", title = "Free response")
ax_2 = Axis(fig[2, 1], xlabel = "Time (s)", ylabel = "Displacement (m)")
lines!(ax_1, t, x_free_gα[1, :], label = "x₁ - Generalized-α")
lines!(ax_1, t, x_free_cd[1, :], label = "x₁ - CentralDiff", linestyle = :dash)
lines!(ax_1, t, x_free_rk[1, :], label = "x₁ - RK4", linestyle = :dashdot)
lines!(ax_1, t, x_free_modal[1, :], label = "x₁ - Modal", linestyle = :dot)
axislegend(ax_1, position = :ct,
           backgroundcolor = (:white, 0.5), orientation = :horizontal)
xlims!(ax_1, minimum(t), maximum(t))

lines!(ax_2, t, x_free_gα[2, :], label = "x₂ - Generalized-α")
lines!(ax_2, t, x_free_cd[2, :], label = "x₂ - CentralDiff", linestyle = :dash)
lines!(ax_2, t, x_free_rk[2, :], label = "x₂ - RK4", linestyle = :dashdot)
lines!(ax_2, t, x_free_modal[2, :], label = "x₂ - Modal", linestyle = :dot)
axislegend(ax_2, position = :cb,
           backgroundcolor = (:white, 0.5), orientation = :horizontal)
xlims!(ax_2, minimum(t), maximum(t))

## Forced response
u0 = (zeros(2), zeros(2))

F0 = 10.
tstart = 2.
duration = 5.
haversine = HaverSine(F0, tstart, duration)
F0 = excitation(haversine, t)
F = zeros(2, length(t))
F[1, :] .= F0

prob_forced = DirectTimeProblem(K, M, C, F, u0, t)
x_forced_gα = solve(prob_forced).u
x_forced_cd = solve(prob_forced, CentralDiff()).u
x_forced_rk = solve(prob_forced, RK4()).u

prob_forced_modal = ForcedModalTimeProblem(K, M, ξ, F, u0, t)
x_forced_modal = solve(prob_forced_modal).u

fig_forced = Figure(size = (600, 600))
ax_forced_1 = Axis(fig_forced[1, 1], ylabel = "Displacement (m)", title = "Forced response")
ax_forced_2 = Axis(fig_forced[2, 1], xlabel = "Time (s)", ylabel = "Displacement (m)")
lines!(ax_forced_1, t, x_forced_gα[1, :], label = "x₁ - Generalized-α")
lines!(ax_forced_1, t, x_forced_cd[1, :], label = "x₁ - CentralDiff", linestyle = :dash)
lines!(ax_forced_1, t, x_forced_rk[1, :], label = "x₁ - RK4", linestyle = :dashdot)
lines!(ax_forced_1, t, x_forced_modal[1, :], label = "x₁ - Modal", linestyle = :dot)
axislegend(ax_forced_1, position = :ct,
           backgroundcolor = (:white, 0.5), orientation = :horizontal)
xlims!(ax_forced_1, minimum(t), maximum(t))

lines!(ax_forced_2, t, x_forced_gα[2, :], label = "x₂ - Generalized-α")
lines!(ax_forced_2, t, x_forced_cd[2, :], label = "x₂ - CentralDiff", linestyle = :dash)
lines!(ax_forced_2, t, x_forced_rk[2, :], label = "x₂ - RK4", linestyle = :dashdot)
lines!(ax_forced_2, t, x_forced_modal[2, :], label = "x₂ - Modal", linestyle = :dot)
axislegend(ax_forced_2, position = :cb,
           backgroundcolor = (:white, 0.5), orientation = :horizontal)
xlims!(ax_forced_2, minimum(t), maximum(t))