using StructuralVibration

L = 1.
b = 0.03
h = 0.01
S = b*h
I = b*h^3/12
E = 2.1e11
ρ = 7850.

# Computation parameters
fmax = 5000.

# Initialization of the data types
beam = Beam(L, S, I, E, ρ)
om_th, k_th = modefreq(beam, fmax)

# Mesh definition
oned_mesh = OneDMesh(beam, 0., 20, :SS)

# Construction of K and M
Kfe, Mfe = assembly(beam, oned_mesh)

# Application of the BCs
Kbc = apply_bc(Kfe, oned_mesh)
Mbc = apply_bc(Mfe, oned_mesh)
nddl = size(Kbc, 1)

# Calculation of the damping matrix
om_fe, ms_fe = eigenmode(Kbc, Mbc)
Cbc = rayleigh_damping_matrix(Kbc, Mbc, 1e-4, 5e-4)
# Cmodal = modal_damping_matrix(Mbc, om_fe, 0.01, ms_fe)

t = 0.:1e-6:0.01
nt = length(t)

# pull middle node
x0 = 1e-3ones(nddl)
v0 = zeros(nddl)
u0 = (x0, v0)

# No force
Fbc = zeros(nddl, nt)
# Direct time problem
prob = DirectTimeProblem(Kbc, Mbc, Cbc, Fbc, u0, t)
# prob = FreeModalTimeProblem(Kbc, Mbc, .01, u0, t)

# u0m = (ms_fe'u0[1], ms_fe'u0[2])
# prob = DirectTimeProblem(ms_fe'Kbc*ms_fe, ms_fe'Mbc*ms_fe, ms_fe'Cmodal*ms_fe, ms_fe'F_free, u0m, t)

res = solve(prob)
res_cd = solve(prob, CentralDiff())
res_rk4 = solve(prob, RK4())

# State-space model
css = ss_model(Kbc, Mbc, Cbc)
u0 = [x0; v0]
prob_ss = StateSpaceTimeProblem(css, Fbc, u0, t)
res_ss = solve(prob_ss, :zoh)