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

# Direct time problem
prob = DirectTimeProblem(Kbc, Mbc, Cbc, u0, t)


u_gα = solve(prob).u
u_cd = solve(prob, CentralDiff()).u

# DifferentialEquations.jl
import ShareAdd as SA
SA.make_importable("OrdinaryDiffEqTsit5", "OrdinaryDiffEqRosenbrock")
import OrdinaryDiffEqTsit5 as Odet
import OrdinaryDiffEqRosenbrock as Oder
function ode_solve!(du, u, p, t)
    A = p[1].Ac

    du .= A*u
end

u0 = [x0; v0]
css = ss_model(Kbc, Mbc, Cbc)
prob_ode = Odet.ODEProblem(ode_solve!, u0, (t[1], t[end]), (css,))
sol_ode = Odet.solve(prob_ode, Odet.AutoTsit5(Oder.Rosenbrock23()))
u_ode = sol_ode(t)[1:nddl, :]

function sde_solve!(ddu, du, u, p, t)
    M = p[1]
    K = p[2]
    C = p[3]

    ddu .= M\(-C*du .- K*u)
end
prob_sde = Odet.SecondOrderODEProblem(sde_solve!, v0, x0, (t[1], t[end]), (Mbc, Kbc, Cbc))
sol_sde = Odet.solve(prob_sde, Odet.AutoTsit5(Oder.Rosenbrock23()))
u_sde = sol_sde(t)[nddl+1:2*nddl, :]