#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
#| output: false
using StructuralVibration, ShareAdd, LinearAlgebra
@usingany CairoMakie
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
@doc  freq_extraction
```
#
#
#
#
#
#
#
#| echo: false
@doc  modeshape_extraction
```
#
#
#
#
#| output: false
# Structure parameters of the beam
L = 1.        # Length
b = 0.03      # Width
h = 0.01      # Thickness
S = b*h       # Cross-section area
Iz = b*h^3/12 # Moment of inertia

# Material parameters
E = 2.1e11  # Young's modulus
ρ = 7850.   # Density
ξ = 0.01    # Damping ratio

# Mesh
xexc = 0:0.1:L
xm = xexc[2]

# Mode calculation - Simply supported boundary conditions
beam = Beam(L, S, Iz, E, ρ)
fmax = 500

ωn, kn = modefreq(beam, 2fmax)
ϕexc = modeshape(beam, kn, xexc)
ϕm = modeshape(beam, kn, xm)

# FRF calculation
freq = 1.:0.1:fmax
prob = ModalFRFProblem(ωn, ξ, freq, ϕm, ϕexc)
H = solve(prob; ismat = true).u

# Natural frequencies and damping ratios extraction
fn_pp, ξn_pp = freq_extraction(freq, H[1, 2, :], PeakPicking())
fn_cf, ξn_cf = freq_extraction(freq, H[1, 2, :], CircleFit())

# Mode shape extraction
ϕid = modeshape_extraction(freq, H, fn, ξn, [1, 2])
#
#
#
#
#
