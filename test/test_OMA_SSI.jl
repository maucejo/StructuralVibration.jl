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

# Mode calculation - Simply supported boundary conditions
beam = Beam(L, S, Iz, E, ρ)
fmax = 500.

ωn, kn = modefreq(beam, 2fmax)
ϕexc = modeshape(beam, kn, xexc)
ϕm = modeshape(beam, kn, xm)

# FRF calculation
# Acquisition parameters for one block
sample_rate = 4096
block_size = 4096
fft_params = FFTParameters(sample_rate, block_size)

id_start = 1
id_end = length(fft_params.freq) - 0
freq = fft_params.freq
freq_calc = freq[id_start:id_end]

# Acquisition parameters
nblocks = 1
tb = fft_params.t
dt = fft_params.dt
t = tb[1]:dt:(nblocks*(tb[end] + dt) - dt)

# Chirp excitation
F0 = 10.
chirp = SweptSine(F0, tb[1], 0.8tb[end], freq[1], freq[end], zero_end = true)
x = repeat(excitation(chirp, tb), outer = nblocks)

force = zeros(length(xexc), length(t))
force[2, :] .= x

prob = ForcedModalTimeProblem(ωn, ϕexc, ξ*ones(length(kn)), ϕexc'force, (zeros(length(xexc)), zeros(length(xexc))), t, ismodal = true)
y = solve(prob).u

# OMA problem definition
tukeywin(x) = tukey(x, 0.5)
prob_oma = OMAProblem(y, t, sample_rate, block_size, win = tukeywin)
p_dssi, ms_dssi = modes_extraction(prob_oma, 10, DataSSI())
# p_cssi, ms_cssi = modes_extraction(prob_oma, 10, CovSSI())

p_dssi = poles_extraction(prob_oma, 28, DataSSI(), stabdiag = true)
p_cssi = poles_extraction(prob_oma, 50, CovSSI(), stabdiag = true)

# EMA-MDOF pole stability analysis
sol_stab = stabilization(prob_oma, 50, DataSSI())

# Plot stabilization diagram
stabilization_plot(sol_stab)