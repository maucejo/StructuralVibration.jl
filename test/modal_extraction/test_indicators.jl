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
xexc = 0:0.05:L
xm = xexc[2]

# Mode calculation - Simply supported boundary conditions
beam = Beam(L, S, Iz, E, ρ)
fmax = 500.

ωn, kn = modefreq(beam, 2fmax)
ms_th = modeshape(beam, kn, xexc)
ms_m = modeshape(beam, kn, xm)

# FRF calculation
freq = 1.:0.1:fmax
prob = ModalFRFProblem(ωn, ξ, freq, ms_m, ms_th)
H = solve(prob).u


# EMA problem
prob_mdof = EMAProblem(H, freq)

# Poles extraction
order = 10 # Model order
p_lscf = poles_extraction(prob_mdof, order, LSCF())

# Driving point indices
dpi = [1, 2]

# Computation of the mode residues
res_lscf = mode_residues(prob_mdof, p_lscf)

# Extraction of the mode shapes
ms, ci = modeshape_extraction(res_lscf, p_lscf, LSCF(), dpi = dpi, modetype = :emar)

# Convert to real mode shapes
ms_real = c2r_modeshape(ms)

# FRF reconstruction - with residuals
lr, ur = compute_residuals(prob_mdof, res_lscf, p_lscf)
H_est = frf_reconstruction(res_lscf, p_lscf, freq, lr = lr, ur = ur)


## Mode complexity indicator
mov_indicator = mov(p_lscf, ms, ci)
mpc_indicator = mpc(ms)
mcf_indicator = mcf(ms)
mpd_indicator = mpd(ms)


## Correlation indicators
mac_indicator = mac(ms[:, 1], ms_th[ :, 1])
mac_matrix = mac(ms, ms_th[:, 1:size(ms, 2)])
comac_indicator = comac(ms, ms_th[:, 1:size(ms, 2)])
ecomac_indicator = ecomac(ms, ms_th[:, 1:size(ms, 2)])
frac_indicator = frac(H_est, H)


## Indicator functions
cmif_indicator = cmif(H)
psif_indicator = psif(H)