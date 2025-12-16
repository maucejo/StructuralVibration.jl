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

t = 0.:1e-5:1.
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

# DifferentialEquations.jl
@usingany OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock
import OrdinaryDiffEqTsit5 as ODET
import OrdinaryDiffEqRosenbrock as ODER
function ode_solve!(du, u, p, t)
    A = p[1].Ac

    du .= A*u
end

u0 = [x0; v0]
css = ss_model(Kbc, Mbc, Cbc)
prob_ode = ODEProblem(ode_solve!, u0, (t[1], t[end]), (css,))
sol_ode = ODET.solve(prob_ode, ODET.AutoTsit5(ODER.Rosenbrock23()))
u_ode = sol_ode(t)[1:nddl, :]

function sde_solve!(ddu, du, u, p, t)
    M = p[1]
    K = p[2]
    C = p[3]

    ddu .= M\(-C*du .- K*u)
end
prob_sde = SecondOrderODEProblem(sde_solve!, v0, x0, (t[1], t[end]), (Mbc, Kbc, Cbc))
sol_sde = ODET.solve(prob_sde, ODET.AutoTsit5(ODER.Rosenbrock23()))
u_sde = sol_sde(t)[nddl+1:2*nddl, :]