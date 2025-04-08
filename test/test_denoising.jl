using StructuralVibration
@usingany GLMakie

# Beam definition
L = 1.
b = 3e-2
h = 1e-2
E = 2.1e11
ρ = 7850.
ξn = 1e-2
S = b*h
Iz = b*h^3/12.

beam = Beam(L, S, Iz, E, ρ)

# Measurement mesh
Δx = 5e-2
Npoint = 20
Xm = LinRange(Δx, L - Δx, Npoint)

# Modes calculation
ωn, kn = modefreq(beam, 2000.)
ϕm = modeshape(beam, kn, Xm)
ϕe = modeshape(beam, kn, Xm[13])

# Modal model
Kn, Mn, Cn = modal_matrices(ωn, ξn)

# Problem definition & solution
tmax = 0.1
nt = 5_000
t = LinRange(0., tmax, nt)

harmo = SineWave(1e4, 0., tmax, 2π*10.)
F = excitation(harmo, t)
Fn = ϕe'*F'

u0 = (zeros(length(ωn)), zeros(length(ωn)))
prob = DirectTimeProblem(Kn, Mn, Cn, Fn, u0, t)
sol = solve(prob)

u = ϕm*sol.u

# Gaussian White Noise
snr_dB = 10.
y = agwn(u, snr_dB)

# Denoising
yc_gcv = denoising(y, GCVDenoising())
yc_lcurve = denoising(y, LCurveDenoising())
yc_kalman = denoising(y, KalmanDenoising())
yc_rts = denoising(y, KalmanDenoising(rts = true))

fig = Figure()
Label(fig[0, 1:2], "Regularization Denoising", fontsize = 22)
ax1 = Axis(fig[1, 1], title = "Node 4", ylabel = "Displacement (m)")
ax2 = Axis(fig[2, 1], title = "Node 13", xlabel = "Time (s)", ylabel = "Displacement (m)")
lines!(ax1, t, y[4, :], color = (:blue, 0.05), label = "Noisy signal", transparency = true)
lines!(ax1, t, u[4, :], label = "Reference")
lines!(ax1, t, yc_gcv[4, :], color = :red, linewidth = 2, linestyle = :dash, label = "GCV")
lines!(ax1, t, yc_lcurve[4, :], color = :green, linewidth = 2, linestyle = :dashdot, label = "LCurve")
xlims!(ax1, 0., tmax)

lines!(ax2, t, y[13, :], color = (:blue, 0.05), label = "Noisy signal", transparency = true)
lines!(ax2, t, u[13, :], label = "Reference")
lines!(ax2, t, yc_gcv[13, :], color = :red, linewidth = 2, linestyle = :dash, label = "GCV")
lines!(ax2, t, yc_lcurve[13, :], color = :green, linewidth = 2, linestyle = :dashdot, label = "LCurve")
xlims!(ax2, 0., tmax)

Legend(fig[1:2, 2], ax1)

fig1 = Figure()
Label(fig1[0, 1:2], "Kalman Denoising", fontsize = 22)
ax3 = Axis(fig1[1, 1], title = "Node 4", ylabel = "Displacement (m)")
ax4 = Axis(fig1[2, 1], title = "Node 13", xlabel = "Time (s)", ylabel = "Displacement (m)")
lines!(ax3, t, y[4, :], color = (:blue, 0.05), label = "Noisy signal", transparency = true)
lines!(ax3, t, u[4, :], label = "Reference")
lines!(ax3, t, yc_kalman[4, :], color = :red, linewidth = 2, linestyle = :dash, label = "Kalman")
lines!(ax3, t, yc_rts[4, :], color = :green, linewidth = 2, linestyle = :dashdot, label = "RTS")
xlims!(ax3, 0., tmax)

lines!(ax4, t, y[13, :], color = (:blue, 0.05), label = "Noisy signal", transparency = true)
lines!(ax4, t, u[13, :], label = "Reference")
lines!(ax4, t, yc_kalman[13, :], color = :red, linewidth = 2, linestyle = :dash, label = "Kalman")
lines!(ax4, t, yc_rts[13, :], color = :green, linewidth = 2, linestyle = :dashdot, label = "RTS")
xlims!(ax4, 0., tmax)

Legend(fig1[1:2, 2], ax3)