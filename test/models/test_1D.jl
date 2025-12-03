## 1D models testing
using StructuralVibration

# Dimensions
L = 1.
d = 3e-2

# Section features
S = π*d^2/4
Iz = π*d^4/64
IG = 2Iz
J = IG

# Tension for string
T = 5000.

# Material
E = 2.1e11
ν = 0.33
G = E/(1 - 2*ν)
ρ = 7800.

# Computation parameters
x = [0.1, 0.9]

# Initialization of the data types
bar = Bar(L, S, E, ρ)
rod = Rod(L, IG, J, G, ρ)
strings = Strings(L, S, T, ρ)
beam = Beam(L, S, Iz, E, ρ)

# Computation of the natural frequencies
om_bar, k_bar = modefreq(bar, 10_000.)
om_rod, k_rod = modefreq(rod, 10_000.)
om_strings, k_strings = modefreq(strings, 100.)
om_beam, k_beam = modefreq(beam, 2000.)

# Computation of the corresponding mode shapes
phi_bar = modeshape(bar, k_bar, x, :CC)
phi_rod = modeshape(rod, k_rod, x, :FF)
phi_strings = modeshape(strings, k_strings, x, :CC)
phi_beam = modeshape(beam, k_beam, x, :SS)