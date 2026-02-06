using StructuralVibration, LinearAlgebra

## Free response
# System parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

# Time vector
t = 0.:1e-2:30.

# Initial conditions
x0 = [0.2, 0.1]
v0 = zeros(2)

# Problem definition - case 1 - Provide the stiffness and mass matrices
u0 = (x0, v0)
prob = FreeModalTimeProblem(K, M, ξ, u0, t)

# Problem definition - case 2 - Provide the squared natural frequencies and mode shapes
ωm, Φm = eigenmode(K, M)
x0m = Φm'*M*x0
v0m = Φm'*M*v0
u0m = (x0m, v0m)
prob_modal = FreeModalTimeProblem(ωm, Φm, ξ, u0m, t)

# Solution
x_free = solve(prob).u
x_free_modal = solve(prob_modal).u


## Forced response - Harmonic excitation
# System parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

# Time vector
t = 0.:1e-2:30.

# Initial conditions
x0 = [0., 1e-4]
v0 = zeros(2)
u0 = (x0, v0)

# Excitation parameters
F = [1., 2.]
freq = 0.5

# Problem definition - case 1 - Provide the stiffness and mass matrices
prob_harmo = HarmonicModalTimeProblem(K, M, ξ, F, 2π*freq, u0, t)

# Problem definition - case 2 - Provide the squared natural frequencies and mode shapes
ωm, Φm = eigenmode(K, M)
x0m = Φm'*M*x0
v0m = Φm'*M*v0
u0m = (x0m, v0m)
Lm = Φm'*F
prob_harmo_modal = HarmonicModalTimeProblem(ωm, Φm, ξ, Lm, 2π*freq, u0m, t)

# Solution
x_harmo = solve(prob_harmo).u
x_harmo_modal = solve(prob_harmo_modal).u


## Forced response - Arbitrary excitation
# System parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

# Time vector
t = 0.:1e-2:30.

# Initial conditions
u0 = (zeros(2), zeros(2))

# Excitation parameters
F0 = 10.
tstart = 2.
duration = 5.
haversine = HaverSine(F0, tstart, duration)
F0 = excitation(haversine, t)
F = zeros(2, length(t))
F[1, :] .= F0

# Problem definition - case 1 - Provide the stiffness and mass matrices
prob_forced = ForcedModalTimeProblem(K, M, ξ, F, u0, t)

# Problem definition - case 2 - Provide the squared natural frequencies and mode shapes
ωm, Φm = eigenmode(K, M)
u0m = (zeros(2), zeros(2))
Lm = Φm'*F
prob_forced_modal = ForcedModalTimeProblem(ωm, Φm, ξ, Lm, u0m, t)

# Solution
x_forced = solve(prob_forced).u
x_forced_modal = solve(prob_forced_modal).u


## Impulse response
# System parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

# Time vector
t = 0.:1e-2:30.

# Impulse response matrix
h = impulse_response(K, M, ξ, t).u