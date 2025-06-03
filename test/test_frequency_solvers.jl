using StructuralVibration, LinearAlgebra
@usingany GLMakie

# Structural parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# Frequency vector
freq = 0.01:0.001:1.

type = :acc
prob_frf = DirectFRFProblem(K, M, C, freq)
H_direct = solve(prob_frf, type = type, ismat = true).u

# Problem definition - Case 2 - Modal
prob_frf_modal = ModalFRFProblem(ω, ξ, freq, Φ, Φ)
H_modal = solve(prob_frf_modal, type = type, ismat = true).u

fig_frf = Figure()
ax_frf11 = Axis(fig_frf[1, 1], ylabel = "FRF (dB)", title = "H₁₁")
ax_frf12 = Axis(fig_frf[1, 2], title = "H₁₂")
ax_frf21 = Axis(fig_frf[2, 1], xlabel = "Frequency (Hz)", ylabel = "FRF (dB)", title = "H₂₁")
ax_frf22 = Axis(fig_frf[2, 2], xlabel = "Frequency (Hz)", title = "H₂₂")

lines!(ax_frf11, freq, 20log10.(abs.(H_direct[1, 1, :])), label = "Direct")
lines!(ax_frf11, freq, 20log10.(abs.(H_modal[1, 1, :])), label = "Modal", linestyle = :dash)
xlims!(ax_frf11, minimum(freq), maximum(freq))
axislegend(ax_frf11, position = :rt, backgroundcolor = (:white, 0.5))

lines!(ax_frf12, freq, 20log10.(abs.(H_direct[1, 2, :])), label = "Direct")
lines!(ax_frf12, freq, 20log10.(abs.(H_modal[1, 2, :])), label = "Modal", linestyle = :dash)
xlims!(ax_frf12, minimum(freq), maximum(freq))

lines!(ax_frf21, freq, 20log10.(abs.(H_direct[2, 1, :])), label = "Direct")
lines!(ax_frf21, freq, 20log10.(abs.(H_modal[2, 1, :])), label = "Modal", linestyle = :dash)
xlims!(ax_frf21, minimum(freq), maximum(freq))

lines!(ax_frf22, freq, 20log10.(abs.(H_direct[2, 2, :])), label = "Direct")
lines!(ax_frf22, freq, 20log10.(abs.(H_modal[2, 2, :])), label = "Modal", linestyle = :dash)
xlims!(ax_frf22, minimum(freq), maximum(freq))

display(GLMakie.Screen(), fig_frf)

## Frequency spectrum
F = zeros(2, length(freq))
F[1, :] .= 10.
type = :acc

# Problem definition - Case 1 - Direct
prob_freq = DirectFreqProblem(K, M, C, F, freq)
y_freq = solve(prob_freq, type = type).u

# Problem definition - Case 2 - Modal
prob_freq_modal = ModalFreqProblem(ω, ξ, Φ'*F, freq, Φ)
y_freq_modal = solve(prob_freq_modal, type = type).u

fig_y = Figure()
ax_y1 = Axis(fig_y[1, 1], ylabel = "Displacement (dB)")
ax_y2 = Axis(fig_y[2, 1], xlabel = "Frequency (Hz)", ylabel = "Displacement (dB)", title = "y₂")

lines!(ax_y1, freq, 20log10.(abs.(y_freq[1, :])), label = "y₁ - Direct")
lines!(ax_y1, freq, 20log10.(abs.(y_freq_modal[1, :])), label = "y₁ - Modal", linestyle = :dash)
xlims!(ax_y1, minimum(freq), maximum(freq))
axislegend(ax_y1, position = :rt, backgroundcolor = (:white, 0.5))

lines!(ax_y2, freq, 20log10.(abs.(y_freq[2, :])), label = "y₂ - Direct")
lines!(ax_y2, freq, 20log10.(abs.(y_freq_modal[2, :])), label = "y₂ - Modal", linestyle = :dash)
xlims!(ax_y2, minimum(freq), maximum(freq))
axislegend(ax_y2, position = :rt, backgroundcolor = (:white, 0.5))

display(GLMakie.Screen(), fig_y)