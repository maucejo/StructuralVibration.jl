using StructuralVibration
@usingany CairoMakie

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

# Acquisition parameters
sample_rate = 4096
block_size = 4096
fft_params = FFTParameters(sample_rate, block_size)

freq = fft_params.freq
t = fft_params.t
dt = fft_params.dt

# Mode calculation - Simply supported boundary conditions
beam = Beam(L, S, Iz, E, ρ)
ωn, kn = modefreq(beam, 2freq[end])
ϕexc = modeshape(beam, kn, xexc)
ϕm = modeshape(beam, kn, xm)

# Chirp excitation
F0 = 10.
chirp = SweptSine(F0, t[1], 0.8t[end], freq[1], freq[end], zero_end = true)
force = zeros(length(xexc), length(t))
force[2, :] .= excitation(chirp, t)

# Response calculation
prob = ForcedModalTimeProblem(ωn, ϕexc, ξ*ones(length(kn)), ϕexc'force, (zeros(length(xexc)), zeros(length(xexc))), t, ismodal = true)
y = solve(prob).u

# OMA problem definition
tukeywin(x) = tukey(x, 0.5)
prob_oma = OMAProblem(y, t, sample_rate, block_size, win = tukeywin)

p_lsce = poles_extraction(prob_oma, 30, LSCE())
p_lscf = poles_extraction(prob_oma, 30, LSCF())
p_plscf = poles_extraction(prob_oma, 30, PLSCF())

# EMA-MDOF pole stability analysis
sol_stab = stabilization(prob_oma, 30, LSCF())

# Plot stabilization diagram
stabilization_plot(sol_stab)

res = mode_residues(prob_oma, p_lscf)
ms_id = modeshape_extraction(res, p_lscf, LSCF(), modetype = :oma)[1]