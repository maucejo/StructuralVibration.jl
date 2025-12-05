using StructuralVibration, LinearAlgebra

# System matrices
m_ss = Diagonal([2., 1.])
k_ss = [6. -2.; -2. 4.]
c_ss = [0.67 -0.11; -0.11 0.39]

# Continuous-time state space from system matrices
css = ss_model(k_ss, m_ss, c_ss)
λ, Ψ = eigenmode(css.Ac)
ω, ξ = modal_parameters(λ)
Ψr = real_normalization(Ψ[1:2, 2:2:end])

# Continuous-time state space from modal information
ωn, ϕn = eigenmode(k_ss, m_ss)
css_modal = ss_modal_model(ωn, 0.01, ϕn)

# Discrete-time state space
dss = c2d(css, 0.01, :foh)