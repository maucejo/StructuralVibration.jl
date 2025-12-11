using StructuralVibration

## ODS from response spectra
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
x = 0:0.05:L
xexc = x[2]

# Mode calculation - Simply supported boundary conditions
beam = Beam(L, S, Iz, E, ρ)
fmax = 500.

ωn, kn = modefreq(beam, 2fmax)
ms_exc = modeshape(beam, kn, xexc)
ms_m = modeshape(beam, kn, x)

# Response spectrum calculation
freq = 1.:0.1:fmax
Fn = repeat(ms_exc', 1, length(freq)) # Frequency-independent force vector
prob = ModalFreqProblem(ωn, ξ, Fn, freq, ms_m)
y = solve(prob).u

pks_indices = [225, 929, 2101, 3743]
ods_complex, freq_peaks = ods(y, freq, pks_indices)
ods_real = real_normalization(ods_complex)

scaling = msf(ods_real, ms_m[:, 1:size(ods_real, 2)])
ods_scaled = ods_real .* scaling'

## ODS from CSD matrices
# Acquisition parameters
sample_rate = 4096
block_size = 4096
fft_params = FFTParameters(sample_rate, block_size)

freq = fft_params.freq
t = fft_params.t
dt = fft_params.dt

# Mode calculation - Simply supported boundary conditions
ωn, kn = modefreq(beam, 2freq[end])
ms_ref = modeshape(beam, kn, x)

# Chirp excitation
F0 = 10.
chirp = SweptSine(F0, t[1], 0.8t[end], freq[1], freq[end], zero_end = true)
force = zeros(length(x), length(t))
force[2, :] .= excitation(chirp, t)

# Response calculation
prob = ForcedModalTimeProblem(ωn, ms_ref, ξ*ones(length(kn)), ms_ref'force, (zeros(length(x)), zeros(length(x))), t, ismodal = true)
y = solve(prob).u

tukeywin(x) = tukey(x, 0.5)
Gyy, freq = csd(y, y, block_size, tukeywin, fs = sample_rate)

frange = [0., 500.]
freq_idx = freq .>= frange[1] .&& freq .<= frange[2]
freq = freq[freq_idx]
Gyy = Gyy[:, :, freq_idx]

pks_indices = [23, 94, 211, 375]
ods_complex, freq_peaks = ods(Gyy, freq, pks_indices, 0.9)
ods_real = real_normalization(ods_complex)

scaling = msf(ods_real, ms_m[:, 1:size(ods_real, 2)])
ods_scaled = ods_real .* scaling'