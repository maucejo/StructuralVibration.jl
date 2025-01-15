using Parameters, ProgressMeter, LinearAlgebra, LazyGrids

includet("../src/models/plate.jl")
includet("../src/models/model.jl")
includet("../src/solvers/frequency_solvers.jl")

plate = Plate(0.6, 0.4, 5e-3, 2.1e11, 7800., 0.3)

## Définition des fréquences jusqu'à fmax
fmax = 1e3
ωₙ, kₙ = eigval(plate, fmax)
Nmodes = length(ωₙ)

## Définition du modèle modal
Kₙ, Mₙ, Cₙ = modal_model(ωₙ, 1e-2)

# Définition de la déformée modale
Nx = 5
Ny = 5
x, y = ndgrid(LinRange(0., plate.L, Nx), LinRange(0., plate.b, Ny))
X = x[:]
Y = y[:]
ϕₑ = eigmode(plate, kₙ, X, Y)

# Calcul des coordonnées généralisées
freq = 1:500
prob = ModalFRFProblem(ωₙ, 1e-2, ϕₑ, ϕₑ)

FRF = solve(prob, freq, :acc)