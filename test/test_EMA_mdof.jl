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
Δx = 0.01
xexc = Δx:0.05:(L - Δx)
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

prob_mdof = EMAMdofProblem(H, freq)
p_lsce = poles_extraction(prob_mdof, 20, LSCE())
p_lscf = poles_extraction(prob_mdof, 20, LSCF())
p_plscf = poles_extraction(prob_mdof, 20, PLSCF())

# EMA-MDOF pole stability analysis
sol_stab = stabilization(prob_mdof, 15, LSCF())

# Plot stabilization diagram
stabilization_plot(sol_stab)

# Mode shape extraction
dpi = [1, 2]
res = mode_residues(prob_mdof, p_lsce)
ϕid, ci = modeshape_extraction(res, p_lsce, dpi, type = :real)
ϕr = c2r_modeshape(ϕid)

# Automatic EMA-MDOF procedure
prob_ema = AutoEMAMdofProblem(prob_mdof, 20, dpi, LSCE())
sol_ema = solve(prob_ema)
fn_ema, ξn_ema = poles2modal(sol_ema.poles)
ϕn_ema = sol_ema.ms

# FRF reconstruction - without residuals
H_mdof = frf_reconstruction(res, p_lscf, freq)

# FRF reconstruction - with residuals
lr, ur = compute_residuals(prob_mdof, res, p_lscf)
H_mdof2 = frf_reconstruction(res, p_lscf, freq, lr = lr, ur = ur)