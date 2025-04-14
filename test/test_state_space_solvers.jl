using StructuralVibration, LinearAlgebra
@usingany GLMakie

# Structural parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# Time vector
t = 0.:1e-2:30.

# State-space model
css = ss_model(K, M, C)

## Free response
# Force
F_free = zeros(2, length(t))

# Initial conditions
x0 = [0.2, 0.1]
v0 = zeros(2)
u0 = [x0; v0]

# Problem definition
prob_free = StateSpaceTimeProblem(css, F_free, u0, t)
x_free_zoh = solve(prob_free).u
x_free_foh = solve(prob_free, :foh).u
x_free_blh = solve(prob_free, :blh).u
# x_free_rk = solve(prob_free, :rk4).u
x_free_rk = solve(prob_free, RK4()).u

prob_free_modal =  FreeModalTimeProblem(K, M, ξ, (x0, v0), t)
x_free_modal = solve(prob_free_modal).u

fig = Figure(size = (600, 600))
ax_1 = Axis(fig[1, 1], ylabel = "Displacement (m)", title = "Free response")
ax_2 = Axis(fig[2, 1], xlabel = "Time (s)", ylabel = "Displacement (m)")
lines!(ax_1, t, x_free_zoh[1, :], label = "x₁ - ZOH")
lines!(ax_1, t, x_free_foh[1, :], label = "x₁ - FOH", linestyle = :dash)
lines!(ax_1, t, x_free_blh[1, :], label = "x₁ - BLH", linestyle = :dashdot)
lines!(ax_1, t, x_free_rk[1, :], label = "x₁ - RK4", linestyle = :dashdotdot)
lines!(ax_1, t, x_free_modal[1, :], label = "x₁ - Modal", linestyle = :dot)
axislegend(ax_1, position = :ct,
           backgroundcolor = (:white, 0.5), orientation = :horizontal)
xlims!(ax_1, minimum(t), maximum(t))

lines!(ax_2, t, x_free_zoh[2, :], label = "x₂ - ZOH")
lines!(ax_2, t, x_free_foh[2, :], label = "x₂ - FOH", linestyle = :dash)
lines!(ax_2, t, x_free_blh[2, :], label = "x₁ - BLH", linestyle = :dashdot)
lines!(ax_2, t, x_free_rk[2, :], label = "x₂ - RK4", linestyle = :dashdotdot)
lines!(ax_2, t, x_free_modal[2, :], label = "x₂ - Modal", linestyle = :dot)
axislegend(ax_2, position = :cb,
           backgroundcolor = (:white, 0.5), orientation = :horizontal)
xlims!(ax_2, minimum(t), maximum(t))

## Forced response
# Force
F0 = 10.
tstart = 2.
duration = 5.
haversine = HaverSine(F0, tstart, duration)
F0 = excitation(haversine, t)
F = zeros(2, length(t))
F[1, :] .= F0

# Initial conditions
x0 = zeros(2)
v0 = zeros(2)
u0 = [x0; v0]

# Problem definition
prob_forced = StateSpaceTimeProblem(css, F, u0, t)
x_forced_zoh = solve(prob_forced).u
x_forced_foh = solve(prob_forced, :foh).u
x_forced_blh = solve(prob_forced, :blh).u
x_forced_rk = solve(prob_forced, :rk4).u

prob_forced_modal =  ForcedModalTimeProblem(K, M, ξ, F, (x0, v0), t)
x_forced_modal = solve(prob_forced_modal).u

fig_forced = Figure(size = (600, 600))
ax_forced_1 = Axis(fig_forced[1, 1], ylabel = "Displacement (m)", title = "Forced response")
ax_forced_2 = Axis(fig_forced[2, 1], xlabel = "Time (s)", ylabel = "Displacement (m)")
lines!(ax_forced_1, t, x_forced_zoh[1, :], label = "x₁ - ZOH")
lines!(ax_forced_1, t, x_forced_foh[1, :], label = "x₁ - FOH", linestyle = :dash)
lines!(ax_forced_1, t, x_forced_blh[1, :], label = "x₁ - BLH", linestyle = :dashdot)
lines!(ax_forced_1, t, x_forced_rk[1, :], label = "x₁ - RK4", linestyle = :dashdotdot)
lines!(ax_forced_1, t, x_forced_modal[1, :], label = "x₁ - Modal", linestyle = :dot)
axislegend(ax_forced_1, position = :ct,
           backgroundcolor = (:white, 0.5), orientation = :horizontal)
xlims!(ax_forced_1, minimum(t), maximum(t))

lines!(ax_forced_2, t, x_forced_zoh[2, :], label = "x₂ - ZOH")
lines!(ax_forced_2, t, x_forced_foh[2, :], label = "x₂ - FOH", linestyle = :dash)
lines!(ax_forced_2, t, x_forced_blh[2, :], label = "x₁ - BLH", linestyle = :dashdot)
lines!(ax_forced_2, t, x_forced_rk[2, :], label = "x₂ - RK4", linestyle = :dashdotdot)
lines!(ax_forced_2, t, x_forced_modal[2, :], label = "x₂ - Modal", linestyle = :dot)
axislegend(ax_forced_2, position = :cb,
           backgroundcolor = (:white, 0.5), orientation = :horizontal)
xlims!(ax_forced_2, minimum(t), maximum(t))

## FRF
freq = 0.01:0.001:1.

# Problem definition - Case 1 - Direct
type = :acc
prob_frf = StateSpaceFRFProblem(css, freq)
H_direct = solve(prob_frf, type = type, ismat = true).u

# Problem definition - Case 2 - Modal
prob_frf_modal = StateSpaceModalFRFProblem(css, freq)
H_modal = solve(prob_frf_modal, type = type, ismat = true).u

fig_frf_ss = Figure()
ax_frf11_ss = Axis(fig_frf_ss[1, 1], ylabel = "FRF (dB)", title = "H₁₁")
ax_frf12_ss = Axis(fig_frf_ss[1, 2], title = "H₁₂")
ax_frf21_ss = Axis(fig_frf_ss[2, 1], xlabel = "Frequency (Hz)", ylabel = "FRF (dB)", title = "H₂₁")
ax_frf22_ss = Axis(fig_frf_ss[2, 2], xlabel = "Frequency (Hz)", title = "H₂₂")

lines!(ax_frf11_ss, freq, 20log10.(abs.(H_direct[1, 1, :])), label = "Direct")
lines!(ax_frf11_ss, freq, 20log10.(abs.(H_modal[1, 1, :])), label = "Modal", linestyle = :dash)
xlims!(ax_frf11_ss, minimum(freq), maximum(freq))
axislegend(ax_frf11_ss, position = :rt, backgroundcolor = (:white, 0.5))

lines!(ax_frf12_ss, freq, 20log10.(abs.(H_direct[1, 2, :])), label = "Direct")
lines!(ax_frf12_ss, freq, 20log10.(abs.(H_modal[1, 2, :])), label = "Modal", linestyle = :dash)
xlims!(ax_frf12_ss, minimum(freq), maximum(freq))

lines!(ax_frf21_ss, freq, 20log10.(abs.(H_direct[2, 1, :])), label = "Direct")
lines!(ax_frf21_ss, freq, 20log10.(abs.(H_modal[2, 1, :])), label = "Modal", linestyle = :dash)
xlims!(ax_frf21_ss, minimum(freq), maximum(freq))

lines!(ax_frf22_ss, freq, 20log10.(abs.(H_direct[2, 2, :])), label = "Direct")
lines!(ax_frf22_ss, freq, 20log10.(abs.(H_modal[2, 2, :])), label = "Modal", linestyle = :dash)
xlims!(ax_frf22_ss, minimum(freq), maximum(freq))

display(GLMakie.Screen(), fig_frf_ss)

## Frequency spectrum
F = zeros(2, length(freq))
F[1, :] .= 10.
type = :acc

# Problem definition - Case 1 - Direct
prob_freq_ss = StateSpaceFreqProblem(css, F, freq)
y_freq_ss = solve(prob_freq_ss, type = type).u

# Problem definition - Case 2 - Modal
prob_freq_modal_ss = StateSpaceModalFreqProblem(css, F, freq)
y_freq_modal_ss = solve(prob_freq_modal_ss, type = type).u

fig_y_ss = Figure()
ax_y1_ss = Axis(fig_y_ss[1, 1], ylabel = "Displacement (dB)")
ax_y2_ss = Axis(fig_y_ss[2, 1], xlabel = "Frequency (Hz)", ylabel = "Displacement (dB)")

lines!(ax_y1_ss, freq, 20log10.(abs.(y_freq_ss[1, :])), label = "y₁ - Direct")
lines!(ax_y1_ss, freq, 20log10.(abs.(y_freq_modal_ss[1, :])), label = "y₁ - Modal", linestyle = :dash)
xlims!(ax_y1_ss, minimum(freq), maximum(freq))
axislegend(ax_y1_ss, position = :rt, backgroundcolor = (:white, 0.5))

lines!(ax_y2_ss, freq, 20log10.(abs.(y_freq_ss[2, :])), label = "y₂ - Direct")
lines!(ax_y2_ss, freq, 20log10.(abs.(y_freq_modal_ss[2, :])), label = "y₂ - Modal", linestyle = :dash)
xlims!(ax_y2_ss, minimum(freq), maximum(freq))
axislegend(ax_y2_ss, position = :rt, backgroundcolor = (:white, 0.5))

display(GLMakie.Screen(), fig_y_ss)