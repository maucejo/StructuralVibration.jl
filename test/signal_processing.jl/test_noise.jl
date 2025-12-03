using StructuralVibration, Statistics

## Noise estimation
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

# Excitation
tmax = 0.5
nt = 10_000
t = LinRange(0., tmax, nt)

harmo = SineWave(1e4, 0., tmax, 2π*10.)
F = excitation(harmo, t)
Fn = ϕe'*F'

# Solution calculation
u0 = (zeros(length(ωn)), zeros(length(ωn)))
prob = DirectTimeProblem(Kn, Mn, Cn, Fn, u0, t)
u = ϕm*solve(prob).u

# Signal corruption - Additive Gaussian White Noise
SNR_ref = 25.
y = agwn(u, SNR_ref)

# Variance estimation - Average over all channels
v1 = varest(y, GCVEst())
SNR_est1 = mean(estimated_SNR(y, v1))

v2 = varest(y, LCurveEst())
SNR_est2 = mean(estimated_SNR(y, v2))

v3 = varest(y, DerricoEst())
SNR_est3 = mean(estimated_SNR(y, v3))

## Denoising
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

# Regularization denoising
yc_gcv = denoising(y, GCVDenoising())
yc_lcurve = denoising(y, LCurveDenoising())

# Kalman denoising
yc_kalman = denoising(y, KalmanDenoising())
yc_rts = denoising(y, KalmanDenoising(rts = true))