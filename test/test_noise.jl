using StructuralVibration

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

##Problem definition & solution
tmax = 0.5
nt = 10000
t = LinRange(0., tmax, nt)

harmo = SineWave(1e4, 0., tmax, 2π*10.)
F = excitation(harmo, t)
Fn = ϕe'*F'

u0 = (zeros(length(ωn)), zeros(length(ωn)))
prob = DirectTimeProblem(Kn, Mn, Cn, Fn, u0, t)
sol = solve(prob)

u = ϕm*sol.u

# Gaussian White Noise
y1 = agwn(u, 25.)                # Noisy signal
vary1 = varest(y1)               # Variance estimation
SNRy1 = estimated_SNR(y1, vary1) # Signal to noise ratio estimation

# Multiplicative noise
y2 = mgwn(u, 25.)                # Noisy signal
vary2 = varest(y2)               # Variance estimation
SNRy2 = estimated_SNR(y2, vary2) # Signal to noise ratio estimation

# Mixed noise
y3 = mixed_noise(u, 25.)         # Noisy signal
vary3 = varest(y3)               # Variance estimation
SNRy3 = estimated_SNR(y3, vary3) # Signal to noise ratio estimation

# Colored noise
y4 = acn(u, 25., 1/h, :pink)     # Noisy signal
vary4 = varest(y4)               # Variance estimation
SNRy4 = estimated_SNR(y4, vary4) # Signal to noise ratio estimation