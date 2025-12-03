using StructuralVibration, LinearAlgebra

## Free response
# System parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# Time vector
t = 0.:1e-2:30.

# Initial conditions
x0 = [0.2, 0.1]
v0 = zeros(2)
u0 = (x0, v0)

# External forces
F_free = zeros(2, length(t))

# Direct time problem
prob_free = DirectTimeProblem(K, M, C, F_free, u0, t)
x_free_gα = solve(prob_free).u
x_free_cd = solve(prob_free, CentralDiff()).u
x_free_rk = solve(prob_free, RK4()).u

# Modal time problem - For comparison
prob_free_modal =  FreeModalTimeProblem(K, M, ξ, u0, t)
x_free_modal = solve(prob_free_modal).u



## Forced response
# System parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# Time vector
t = 0.:1e-2:30.

# Initial conditions
x0 = zeros(2)
v0 = zeros(2)
u0 = (x0, v0)

# External forces
F0 = 10.
tstart = 2.
duration = 5.
haversine = HaverSine(F0, tstart, duration)
F0 = excitation(haversine, t)
F = zeros(2, length(t))
F[1, :] .= F0

# Direct time problem
prob_forced = DirectTimeProblem(K, M, C, F, u0, t)
x_forced_gα = solve(prob_forced).u
x_forced_cd = solve(prob_forced, CentralDiff()).u
x_forced_rk = solve(prob_forced, RK4()).u

# Modal time problem - For comparison
prob_forced_modal = ForcedModalTimeProblem(K, M, ξ, F, u0, t)
x_forced_modal = solve(prob_forced_modal).u