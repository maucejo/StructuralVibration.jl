using StructuralVibration

# Dimensions
Lp = 0.6
bp = 0.4
hp = 1e-3

# Material parameters
E = 2.1e11
ρ = 7800.
ν = 0.33

# Computation parameters
fmax = 500.
xp = [0.1, 0.5]
yp = [0.1, 0.3]

# Initialization of the data types
plate = Plate(Lp, bp, hp, E, ρ, ν)
membrane = Membrane(Lp, bp, ρ*hp, 100.)

# Computation of the natural frequencies
om_plate, k_plate = modefreq(plate, fmax)[1:2]
om_membrane, k_membrane = modefreq(membrane, fmax)[1:2]

# Computation of the corresponding mode shapes
phi_plate = modeshape(plate, k_plate, xp, yp)
phi_membrane = modeshape(membrane, k_membrane, xp, yp)

# Computation of the wave parameters
om_plate, c_plate, k_plate, λ_plate = wave_parameters(plate, fmax)
om_membrane, c_membrane, k_membrane, λ_membrane = wave_parameters(membrane, fmax)