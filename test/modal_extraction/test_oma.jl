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
x = 0:0.05:L

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
ms_ref = modeshape(beam, kn, x)

# Chirp excitation
F0 = 10.
chirp = SweptSine(F0, t[1], 0.8t[end], freq[1], freq[end], zero_end = true)
force = zeros(length(x), length(t))
force[2, :] .= excitation(chirp, t)

# Response calculation
prob = ForcedModalTimeProblem(ωn, ms_ref, ξ*ones(length(kn)), ms_ref'force, (zeros(length(x)), zeros(length(x))), t, ismodal = true)
y = solve(prob).u


## EMA-based OMA
# OMA problem definition
tukeywin(x) = tukey(x, 0.5)
prob_oma = OMAProblem(y, t, sample_rate, block_size, win = tukeywin, frange = [1., 1000.])

# EMA-MDOF pole stability analysis
sol_stab = stabilization(prob_oma, 30, LSCF())

# Extraction of stable poles
stab_poles = sol_stab.poles[end][.!isnan.(sol_stab.poles[end])]

# Plot stabilization diagram
stabilization_plot(sol_stab)

# Computation of the mode residues
res = mode_residues(prob_oma, stab_poles)

# Extraction of the mode shapes
ms_ema = modeshape_extraction(res, stab_poles, LSCF(), modetype = :oma)[1]

# Convert to real mode shapes
ms_ema_real = real_normalization(ms_ema)

# Mode shape scaling for comparison using the modal scale factor
scaling = msf(ms_ema_real, ms_ref[:, 1:size(ms_ema_real, 2)])
ms_ema_scaled = ms_ema_real .* scaling'

# Half-spectrum reconstruction - without residuals
Syy = half_csd_reconstruction(res, stab_poles, prob_oma.freq)

# Half-spectrum reconstruction - with residuals
lr, ur = compute_residuals(prob_oma, res, stab_poles)
Syy2 = half_csd_reconstruction(res, stab_poles, prob_oma.freq, lr = lr, ur = ur)


## CovSSI
# Stabilization analysis
sol_stab_cssi = stabilization(prob_oma, 30, CovSSI())

# Plot stabilization diagram
stabilization_plot(sol_stab_cssi)

# Extraction of stable modes
p_cssi, ms_cssi = modes_extraction(prob_oma, 30, CovSSI())
ms_cssi_real = real_normalization(ms_cssi)

# Mode shape scaling for comparison using the modal scale factor
scaling = msf(ms_cssi_real, ms_ref[:, 1:size(ms_cssi_real, 2)])
ms_cssi_scaled = ms_cssi_real .* scaling'


## DataSSI
# Stabilization analysis
sol_stab_dssi = stabilization(prob_oma, 30, DataSSI(), progress = true)

# Plot stabilization diagram
stabilization_plot(sol_stab_dssi)

# Extraction of stable modes
p_dssi, ms_dssi = modes_extraction(prob_oma, 30, DataSSI())
ms_dssi_real = real_normalization(ms_dssi)

# Mode shape scaling for comparison using the modal scale factor
scaling = msf(ms_dssi_real, ms_ref[:, 1:size(ms_dssi_real, 2)])
ms_dssi_scaled = ms_dssi_real .* scaling'