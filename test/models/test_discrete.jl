using StructuralVibration

## Sdof
m = 1.
f₀ = 10.
ξ = 0.01
sdof = Sdof(m, f₀, ξ)

## Mdof
# Definition of the structural parameters
k_mdof = [1., 1.]
m_mdof = ones(3)
c_mdof = [0.1, 0.1]

# Initialization of Mdof
mdof = Mdof(k_mdof, m_mdof, c_mdof)

# Definition of a MdofMesh
mdof_mesh = MdofMesh(mdof, :CF)

# System assembly
K_mdof, M_mdof, C_mdof = assembly(mdof)

# Apply boundary conditions (if any)
K_bc = apply_bc(K_mdof, mdof_mesh)
M_bc = apply_bc(M_mdof, mdof_mesh)
C_bc = apply_bc(C_mdof, mdof_mesh)

# Compute the eigenmodes of the systems
ωn, Φn = eigenmode(K_bc, M_bc)

# Computation of the modal matrices
Kmodal, Mmodal, Cmodal = modal_matrices(ωn, 0.01)

# Computation of modal effective mass
meff = modal_effective_mass(M_bc, Φn, ones(2))

## FE model
# Dimensions
L = 1.
d = 3e-2

# Section features
S = π*d^2/4
Iz = π*d^4/64

# Material
E = 2.1e11
ρ = 7800.

# Computation parameters
fmax = 2000.

# Initialization of the data types
beam = Beam(L, S, Iz, E, ρ)

# Mesh definition
oned_mesh = OneDMesh(beam, 0., 20, :SS)

# Construction of K and M
Kfe, Mfe = assembly(beam, oned_mesh)

# Application of the BCs
Kbc = apply_bc(Kfe, oned_mesh)
Mbc = apply_bc(Mfe, oned_mesh)

# Computation ofthe eigenmodes of the structure
ωfe, Φfe = eigenmode(Kbc, Mbc)

# Calculation of the damping matrix
Cray = rayleigh_damping_matrix(Kbc, Mbc, 1., 1.)
Cmodal = modal_damping_matrix(Mbc, ωfe, 0.01, Φfe)