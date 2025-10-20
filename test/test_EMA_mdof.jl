using StructuralVibration

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
Δx = 0.01
xexc = Δx:0.05:(L - Δx)
xm = xexc[2]

# Mode calculation - Simply supported boundary conditions
beam = Beam(L, S, Iz, E, ρ)

# (; freq) = FFTParameters(1024, 4096)
# fmax = freq[end]
fmax = 500.

ωn, kn = modefreq(beam, 2fmax)
ϕexc = modeshape(beam, kn, xexc)
ϕm = modeshape(beam, kn, xm)

# FRF calculation
freq = 1.:0.1:fmax
prob = ModalFRFProblem(ωn, ξ, freq, ϕm, ϕexc)
H = solve(prob; ismat = true).u

p_lsce = lsce(H, freq, 20)
p_lscf = lscf(H, freq, 20)
p_plscf = plscf(H, freq, 20)

# EMA-MDOF pole stability analysis
sol_stab = stabilization(H, freq, 15, LSCE())
