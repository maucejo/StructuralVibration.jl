using StructuralVibration, LinearAlgebra

## Free response
# Structural parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# Time vector
t = 0.:1e-2:30.

# State-space model
css = ss_model(K, M, C)

# Force
F_free = zeros(2, length(t))

# Initial conditions
x0 = [0.2, 0.1]
v0 = zeros(2)
u0 = [x0; v0]

# Problem definition
prob_free = StateSpaceTimeProblem(css, F_free, u0, t)
x_free_zoh = solve(prob_free).u
x_free_foh = solve(prob_free, :foh).u
x_free_blh = solve(prob_free, :blh).u
x_free_rk = solve(prob_free, RK4()).u
# Other possibility
x_free_rk2 = solve(prob_free, :rk4).u

# Modal free response - For comparison
prob_free_modal =  FreeModalTimeProblem(K, M, ξ, (x0, v0), t)
x_free_modal = solve(prob_free_modal).u



## Forced response
# Structural parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# Time vector
t = 0.:1e-2:30.

# State-space model
css = ss_model(K, M, C)

F0 = 10.
tstart = 2.
duration = 5.
haversine = HaverSine(F0, tstart, duration)
F0 = excitation(haversine, t)
F = zeros(2, length(t))
F[1, :] .= F0

# Initial conditions
x0 = zeros(2)
v0 = zeros(2)
u0 = [x0; v0]

# Problem definition
prob_forced = StateSpaceTimeProblem(css, F, u0, t)
x_forced_zoh = solve(prob_forced).u
x_forced_foh = solve(prob_forced, :foh).u
x_forced_blh = solve(prob_forced, :blh).u
x_forced_rk = solve(prob_forced, :rk4).u

# Modal forced response - For comparison
prob_forced_modal = ForcedModalTimeProblem(K, M, ξ, F, (x0, v0), t)
x_forced_modal = solve(prob_forced_modal).u


## Frequency response function
# Structural parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# State-space model
css = ss_model(K, M, C)

# Frequency vector
freq = 0.01:0.001:1.

# Problem definition - Case 1 - Direct
prob_frf = StateSpaceFRFProblem(css, freq)
H_direct = solve(prob_frf).u

# Problem definition - Case 2 - Modal
prob_frf_modal = StateSpaceModalFRFProblem(css, freq)
H_modal = solve(prob_frf_modal).u


## Response spectrum
# Structural parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# State-space model
css = ss_model(K, M, C)

# Frequency vector
freq = 0.01:0.001:1.

# Force matrix
F = zeros(2, length(freq))
F[1, :] .= 10.

# Problem definition - Case 1 - Direct
prob_freq = StateSpaceFreqProblem(css, F, freq)
y_freq = solve(prob_freq).u

# Problem definition - Case 2 - Modal
prob_freq_modal = StateSpaceModalFreqProblem(css, F, freq)
y_freq_modal = solve(prob_freq_modal).u