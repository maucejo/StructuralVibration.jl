using StructuralVibration
@usingany GLMakie

# Structure parameters of the beam
L = 1.        # Length
b = 0.03      # Width
h = 0.01      # Thickness
S = b*h       # Cross-section area
Iz = b*h^3/12 # Moment of inertia

# Material parameters
E = 2.1e11  # Young's modulus
ρ = 7850.   # Density
ξ = 0.01    # Damping ratio

# Mesh
xexc = 0:0.05:L
xm = xexc[2]

# Mode calculation - Simply supported boundary conditions
beam = Beam(L, S, Iz, E, ρ)
fmax = 500

ωn, kn = modefreq(beam, 2fmax)
ϕexc = modeshape(beam, kn, xexc)
ϕm = modeshape(beam, kn, xm)

# FRF calculation
freq = 1.:0.1:fmax
prob = ModalFRFProblem(ωn, ξ, freq, ϕm, ϕexc)
H = solve(prob; ismat = true).u

# Natural frequencies and damping ratios extraction
fn_pp, ξn_pp = freq_extraction(freq, H[1, 2, :], PeakPicking())
fn_cf, ξn_cf = freq_extraction(freq, H[1, 2, :], CircleFit())

# Mode shape extraction
ϕid = modeshape_extraction(freq, H, fn_pp, ξn_pp, [1, 2])

# Plotting
fig_f = Figure()
ax_f1 = Axis(fig_f[1, 1], title = "Natural frequencies", xlabel = "Mode ID", ylabel = "Natural frequency [Hz]")
ax_f2 = Axis(fig_f[1, 2], title = "Damping ratios", xlabel = "Mode ID", ylabel = "Damping ratio  (%)")

scatter!(ax_f1, 1:4, ωn[1:4]/2π, marker = :rect, markersize = 25, label = "Reference")
scatter!(ax_f1, 1:4, fn_pp, markersize = 20, label = "Peak picking")
scatter!(ax_f1, 1:4, fn_cf, marker = :star4, markersize = 15, label = "Circle fit")
axislegend(ax_f1, position = :rb)

scatter!(ax_f2, 1:4, 100ξ*ones(4), marker = :rect, markersize = 25, label = "Reference")
scatter!(ax_f2, 1:4, 100ξn_pp, markersize = 20, label = "Peak picking")
scatter!(ax_f2, 1:4, 100ξn_cf, marker = :star4, markersize = 15, label = "Circle fit")

fig_mode = Figure()
ax_mode1 = Axis(fig_mode[1, 1], title = "Mode shape 1", ylabel = "Value")
ax_mode2 = Axis(fig_mode[1, 2], title = "Mode shape 2")
ax_mode3 = Axis(fig_mode[2, 1], title = "Mode shape 3", xlabel = "Position [m]", ylabel = "Value")
ax_mode4 = Axis(fig_mode[2, 2], title = "Mode shape 4", xlabel = "Position [m]")

lines!(ax_mode1, xexc, ϕexc[:, 1], label = "Reference")
lines!(ax_mode1, xexc, ϕid[:, 1], linestyle = :dash, label = "Estimated")
xlims!(ax_mode1, 0, L)
axislegend(ax_mode1, position = :cb)

lines!(ax_mode2, xexc, ϕexc[:, 2], label = "Reference")
lines!(ax_mode2, xexc, ϕid[:, 2], linestyle = :dash, label = "Estimated")
xlims!(ax_mode2, 0, L)

lines!(ax_mode3, xexc, ϕexc[:, 3], label = "Reference")
lines!(ax_mode3, xexc, ϕid[:, 3], linestyle = :dash, label = "Estimated")
xlims!(ax_mode3, 0, L)

lines!(ax_mode4, xexc, ϕexc[:, 4], label = "Reference")
lines!(ax_mode4, xexc, ϕid[:, 4], linestyle = :dash, label = "Estimated")
xlims!(ax_mode4, 0, L)