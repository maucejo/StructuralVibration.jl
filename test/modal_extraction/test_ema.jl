using StructuralVibration, Peaks, Statistics
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
ms_exc = modeshape(beam, kn, xexc)
ms_m = modeshape(beam, kn, xm)

# FRF calculation
freq = 1.:0.1:fmax
prob = ModalFRFProblem(ωn, ξ, freq, ms_m, ms_exc)
H = solve(prob).u

## Sdof methods
# Natural frequencies and damping ratios extraction
prob_sdof = EMAProblem(H, freq)
poles_pp = poles_extraction(prob_sdof, PeakPicking())
poles_cf = poles_extraction(prob_sdof, CircleFit())
poles_lsf = poles_extraction(prob_sdof, LSFit())

pks = findmaxima(vec(mean(abs, H, dims = 2)))
poles_pp2 = poles_extraction(prob_sdof, PeakPicking(), pks_indices = pks.indices)
poles_cf2 = poles_extraction(prob_sdof, CircleFit(), pks_indices = pks.indices)
poles_lsf2 = poles_extraction(prob_sdof, LSFit(), pks_indices = pks.indices)

# Mode shape extraction
dpi = [1, 2]
ms_id = modeshape_extraction(prob_sdof, poles_pp, PeakPicking(), dpi = dpi)

# Computation of the mode residues
res_pp = mode2residues(ms_id, poles_pp, [2], 1:21)[1]

# FRF reconstruction - without residuals
H_sdof = frf_reconstruction(res_pp, poles_pp, freq)

# FRF reconstruction - with residuals
lr, ur = compute_residuals(prob_sdof, res_pp, poles_pp)
H_sdof2 = frf_reconstruction(res_pp, poles_pp, freq, lr = lr, ur = ur)

# prob_sdof_ema = AutoEMASdofProblem(prob_sdof, PeakPicking(), dpi = dpi, idx_m = [2], idx_e = 1:21)
# sol_ema = solve(prob_sdof_ema)
sol_sdof_ema = solve(prob_sdof, PeakPicking(), dpi = dpi, idx_m = [2], idx_e = 1:21)
fn_sdof_ema, ξn_sdof_ema = poles2modal(sol_sdof_ema.poles)
ms_sdof_ema = sol_sdof_ema.ms


## Mdof methods
# EMA problem
prob_mdof = EMAProblem(H, freq)

# Poles extraction
order = 10 # Model order

p_lsce = poles_extraction(prob_mdof, order, LSCE())
p_lscf = poles_extraction(prob_mdof, order, LSCF())
p_plscf = poles_extraction(prob_mdof, order, pLSCF())

# Stabilization diagram analysis using the LSCF method
stab = stabilization(prob_mdof, order, LSCF())

# Visualization of the stabilization diagram
stabilization_plot(stab)

# Driving point indices
dpi = [1, 2]

# Computation of the mode residues
res_lscf = mode_residues(prob_mdof, p_lscf)

# Extraction of the mode shapes
ms_est = modeshape_extraction(res_lscf, p_lscf, LSCF(), dpi = dpi, modetype = :emar)[1]

# Convert to real mode shapes
ms_est_real = real_normalization(ms_est)

# FRF reconstruction - without residuals
H_mdof = frf_reconstruction(res_lscf, p_lscf, freq)

# FRF reconstruction - with residuals
lr, ur = compute_residuals(prob_mdof, res_lscf, p_lscf)
H_mdof2 = frf_reconstruction(res_lscf, p_lscf, freq, lr = lr, ur = ur)

# prob_mdof_ema = AutoEMAMdofProblem(prob_mdof, 10, LSCF(), dpi = dpi)
sol_mdof_ema = solve(prob_mdof, 10, LSCF(), dpi = dpi)