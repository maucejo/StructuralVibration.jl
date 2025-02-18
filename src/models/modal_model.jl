"""
    modal_matrices(ωₙ, ξₙ)

Computes the modal mass, stiffness, and damping matrices

**Inputs**
* `ωₙ`: Vector of natural frequencies
* `ξₙ`: Modal damping

**Outputs**
* `Kₙ`: Generalized stiffness matrix
* `Mₙ`: Generalized mass matrix (identity matrix, due to mass normalization)
* `Cₙ`: Generalized damping matrix
"""
function modal_matrices(ωₙ, ξₙ)
    return Diagonal(ωₙ.^2), I(length(ωₙ)), Diagonal(2ξₙ.*ωₙ)
end

"""
    modal_effective_mass(M, ϕₙ, r)

Computes the effective mass of a mode

**Inputs**
* `M`: Mass matrix
* `ϕₙ`: Mode shape
* `r`: Influence vector (rigid body mode)

**Outputs**
* `meff`: Modal effective mass

*Note: The modeshapes are supposed to be mass-normalized*
"""
function modal_effective_mass(M, ϕₙ, r)
    Lₙ = ϕₙ'*M*r

    return Lₙ'*Lₙ
end