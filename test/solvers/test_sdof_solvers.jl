using StructuralVibration

## Free response
m = 1.
f0 = 1.

# Time vector
t = 0.:0.01:10.

# Initial conditions
u0 = [1., -2.]

# Undamped system
sdof_nd = Sdof(m, f0, 0.)
prob_nd = SdofFreeTimeProblem(sdof_nd, u0, t)
x_nd = solve(prob_nd).u

# Underdamped system
sdof_ud = Sdof(m, f0, 0.1)
prob_ud = SdofFreeTimeProblem(sdof_ud, u0, t)
x_ud = solve(prob_ud).u

# Critically damped system
sdof_cd = Sdof(m, f0, 1.)
prob_cd = SdofFreeTimeProblem(sdof_cd, u0, t)
x_cd = solve(prob_cd).u

# Overdamped system
sdof_od = Sdof(m, f0, 1.2)
prob_od = SdofFreeTimeProblem(sdof_od, u0, t)
x_od = solve(prob_od).u



## Harmonic response
# Structural parameters
m = 1.
f0 = 1.
ξ = 0.1

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

# Harmonic excitation - Force amplitude + excitation frequency
F0 = 10.
f = 2.

# Instantiation of the problem
t = 0.:0.01:10.
u0 = [1., -2.]
prob = SdofHarmonicTimeProblem(sdof, u0, t, F0, f, :force)

# Solve the problem
x_harmo = solve(prob).u

## Forced response
# Structural parameters
m = 1.
f0 = 1.
ξ = 0.1

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

# Haversine excitation signal
F0 = 10.
tstart = 0.5
duration = 2.
haversine = HaverSine(F0, tstart, duration)

t = 0.:0.01:10.
F = excitation(haversine, t)

# Instantiation of the problem
u0 = [1., -2.]
prob = SdofForcedTimeProblem(sdof, u0, t, F, :force)

# Solve the problem
x_arb_filt = solve(prob, method = :filt).u
x_arb_interp = solve(prob, method = :interp).u
x_arb_conv = solve(prob, method = :conv).u

## SRS calculation
# Definition of the base acceleration
t = 0.:1e-5:0.5
base_acc = HalfSine(10., 0., 1e-3)

# Compute the SRS
f = 50.:25:1e3
srs_basic = srs(base_acc, f, t)
srs_rec_int = srs(base_acc, f, t, alg = :RecursiveInt)
srs_rec_filt = srs(base_acc, f, t, alg = :RecursiveFilt)
srs_smallwood = srs(base_acc, f, t, alg = :Smallwood)


## Impulse response
# Structural parameters
m = 1.
f0 = 1.
ξ = 0.1

t = 0.:0.01:10.

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

h = impulse_response(sdof, t)


## Frequency response function
# Structural parameters
m = 1.
f0 = 10.
ξ = 0.01

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

# Definition of the frequency range
f = 1.:0.1:30.

# Instantiation of the problem - Force excitation
prob = SdofFRFProblem(sdof, f)

# Solve the problem
H = solve(prob).u



## Response spectrum
# Structural parameters
m = 1.
f0 = 10.
ξ = 0.01

# Instantiation of the Sdof object
sdof = Sdof(m, f0, ξ)

# Definition of the frequency range
f = 1.:0.1:30.

# Instantiation of the problem - Force excitation
F = fill(10., length(f))
prob = SdofFrequencyProblem(sdof, f, F)

# Solve the problem
y = solve(prob).u