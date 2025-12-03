using StructuralVibration

## Frequency response function
# Structural parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# Frequency vector
freq = 0.01:0.001:1.

prob_frf = DirectFRFProblem(K, M, C, freq)
H_direct = solve(prob_frf).u

# Problem definition - Case 2 - Modal
prob_frf_modal = ModalFRFProblem(ω, ξ, freq, Φ, Φ)
H_modal = solve(prob_frf_modal).u

## Response spectrum
# Structural parameters
M = Diagonal([2., 1.])
K = [6. -2.; -2. 4.]
ξ = 0.05

ω, Φ = eigenmode(K, M)
C = modal_damping_matrix(M, ω, ξ, Φ)

# Frequency vector
freq = 0.01:0.001:1.

# Force matrix
F = zeros(2, length(freq))
F[1, :] .= 10.

# Problem definition - Case 1 - Direct
prob_freq = DirectFreqProblem(K, M, C, F, freq)
y_freq = solve(prob_freq).u

# Problem definition - Case 2 - Modal
prob_freq_modal = ModalFreqProblem(ω, ξ, Φ'F, freq, Φ)
y_freq_modal = solve(prob_freq_modal).u