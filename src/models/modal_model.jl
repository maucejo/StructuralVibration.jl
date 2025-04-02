"""
    modal_matrices(ωn, ξn)

Computes the modal mass, stiffness, and damping matrices

**Inputs**
* `ωn`: Vector of natural frequencies
* `ξn`: Modal damping

**Outputs**
* `Kn`: Generalized stiffness matrix
* `Mn`: Generalized mass matrix (identity matrix, due to mass normalization)
* `Cn`: Generalized damping matrix
"""
function modal_matrices(ωn, ξn)
    return Diagonal(ωn.^2), I(length(ωn)), Diagonal(2ξn.*ωn)
end

"""
    modal_effective_mass(M, ϕn, r)

Computes the effective mass of a mode

**Inputs**
* `M`: Mass matrix
* `ϕn`: Mode shape
* `r`: Influence vector (rigid body mode)

**Outputs**
* `meff`: Modal effective mass

*Note: The modeshapes are supposed to be mass-normalized*
"""
function modal_effective_mass(M, ϕn, r)
    Ln = ϕn'*M*r

    return Ln'*Ln
end