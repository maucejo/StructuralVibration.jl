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
ϕexc = modeshape(beam, kn, xexc)
ϕm = modeshape(beam, kn, xm)

# FRF calculation
freq = 1.:0.1:fmax
prob = ModalFRFProblem(ωn, ξ, freq, ϕm, ϕexc)
H = solve(prob; ismat = true).u

# Natural frequencies and damping ratios extraction
prob_sdof = EMAProblem(H, freq)
poles_pp = poles_extraction(prob_sdof, PeakPicking())
poles_cf = poles_extraction(prob_sdof, CircleFit())
poles_lsf = poles_extraction(prob_sdof, LSFit())

# Mode shape extraction
dpi = [1, 2]
ms_id = modeshape_extraction(prob_sdof, poles_pp, PeakPicking(), dpi = dpi)

# Automatic EMA
prob_ema = AutoEMASdofProblem(prob_sdof, PeakPicking(), dpi = dpi, idx_m = [2], idx_e = 1:21)
sol_ema = solve(prob_ema)
fn_ema, ξn_ema = poles2modal(sol_ema.poles)
ϕn_ema = sol_ema.ms

# FRF calculation
res_pp = mode2residues(ms_id, poles_pp, [2], 1:21)[1]
H_sdof = frf_reconstruction(res_pp, poles_pp, freq)

# FRF reconstruction - with residuals
lr, ur = compute_residuals(prob_sdof, res_pp, poles_pp)
H_sdof2 = frf_reconstruction(res_pp, poles_pp, freq, lr = lr, ur = ur)