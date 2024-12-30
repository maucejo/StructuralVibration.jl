"""
    modal_model(ωₙ, ηₙ)

Computes the modal mass, stiffness, and damping matrices

# Parameters
    * p : Structure containing data related to the plate
    * ωₙ : Vector of natural frequencies
    * ξₙ : Modal damping

# Outputs
    * Kₙ : Generalized stiffness matrix
    * Mₙ : Generalized mass matrix (identity matrix, due to mass normalization)
    * Cₙ : Generalized damping matrix
"""
function modal_model(ωₙ, ξₙ)
    return Diagonal(ωₙ.^2), I(length(ωₙ)), Diagonal(2ξₙ.*ωₙ)
end