using StructuralVibration

plaq = Plate(0.6, 0.4, 5e-3, 2.1e11, 7800., 0.3)

## Définition des fréquences jusqu'à fmax
fmax = 1e3
ωₙ, kₙ = eigval(plaq, fmax)
Nmodes = length(ωₙ)

## Définition du modèle modal
Kₙ, Mₙ, Cₙ = modal_model(ωₙ, 1e-2)

# Définition de la déformée modale
loc = [0.1, 0.1]
ϕₑ = eigmode(plaq, kₙ, loc[1], loc[2])

# Calcul des coordonnées généralisées
freq = 1:500
prob = ModalFRF(ωₙ, 1e-2, ϕₑ, ϕₑ, freq)

FRF = frf(prob, :acc)