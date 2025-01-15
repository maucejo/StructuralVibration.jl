using Parameters, DSP, LinearAlgebra, ProgressMeter, Interpolations
@usingany GLMakie

includet("../src/models/excitation.jl")
includet("../src/solvers/modal_time_solvers.jl")
includet("../src/utils/calculus.jl")

## Geometrical parameters
r = 0.06
h = 0.6
S = π*r^2

## Material parameters
E = 70e9
ρ = 2700.
ξ = 1e-3

## Stiffness and mass matrices
K = (E*S/h).*[(1. + √2/4.) -√2/4.; -√2/4. √2/4.]
M = (ρ*S*h).*[(1. + √2/2.)/3. -√2/6.; -√2/6. √2/6.];

## Solver
tmax = 0.5
nt = 10000
t = LinRange(0., tmax, nt)

## Free response
# u0 = ([0.2, 0.1], zeros(2))
# prob = FreeModalTimeProblem(K, M, ξ, u0, t)

## Harmonic response
# u0 = ([0., 1e-4], zeros(2))
# prob = HarmonicModalTimeProblem(K, M, ξ, u0, t, [1e4, 0.], 2π*1000.)

## Forced response
u0 = ([0., 1e-4], zeros(2))
harmo = SineWave(1e4, 0., tmax, 2π*1000., π/2.)
F₀ = excitation(harmo, t)
F = zeros(2, nt)
F[1, :] .= F₀
prob = ForcedModalTimeProblem(K, M, ξ, u0, t, F)

sol = solve(prob)