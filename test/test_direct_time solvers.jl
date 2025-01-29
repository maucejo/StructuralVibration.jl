using LazyGrids
using StructuralVibration

## Plate definition
plate = Plate(0.6, 0.4, 5e-3, 2.1e11, 7800., 0.3)

## Computation of the resonance frequencies
fmax = 1e3
ωₙ, kₙ = modefreq(plate, 10fmax)
Nmodes = length(ωₙ)

## Construction of the modale model
Kₙ, Mₙ, Cₙ = modal_matrices(ωₙ, 1e-2)

## Calculation of the modal force vector
Δt = 1e-6 # Pas de temps
tf = 0.07 # instant final
t = 0.:Δt:tf
loc = [0.1, 0.2]

hammer = Hammer(1., 8e-3, 9.7, 6e-4)
F = excitation(hammer, t)
ϕₑ = modeshape(plate, kₙ, loc[1], loc[2])
Fₙ = (F*ϕₑ)'

## Listening points mesh
Nx = 40
Ny = 40
x, y = ndgrid(LinRange(0., plate.L, Nx), LinRange(0., plate.b, Ny))
X = x[:]
Y = y[:]
ϕₒ = modeshape(plate, kₙ, X, Y)

# Computation of the modal coordinates
u0 = (zeros(Nmodes), zeros(Nmodes))
prob = DiscreteTimeProblem(Kₙ, Mₙ, Cₙ, u0, Δt, Fₙ)

sol = solve(prob)
(; D, V, A) = sol

## Computation of the displacement, velocity and acceleration at the listening points
Disp = ϕₒ*D
Vel = ϕₒ*V
Acc = ϕₒ*A