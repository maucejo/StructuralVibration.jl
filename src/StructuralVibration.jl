module StructuralVibration

using DSP, FastGaussQuadrature, FFTW, Interpolations, LinearAlgebra, Optim, Parameters, Peaks, PrecompileTools, ProgressMeter, Random, SpecialFunctions, Statistics

# Structs - Models
export Sdof, Mdof, Bar, Rod, Strings, Beam, Plate, Membrane,
       ContinuousStateSpace, DiscreteStateSpace

# Structs - FE and discrete models
export OneDMesh, MdofMesh

# Structs - Excitations
export Rectangle, Triangle, Hammer, SmoothRect, SineWave,
       HalfSine, HaverSine,SweptSine, GaussianPulse, ColoredNoise

# Structs - Problems
export SdofFreeTimeProblem, SdofHarmonicTimeProblem, SdofForcedTimeProblem,
       SdofFRFProblem, SdofFrequencyProblem, StateSpaceTimeProblem, StateSpaceFRFProblem, StateSpaceFreqProblem,
       FreeModalTimeProblem, ForcedModalTimeProblem, ModalFRFProblem, DirectFRFProblem, ModalFreqProblem, DirectFreqProblem, DirectTimeProblem

# Structs - Time solvers
export CentralDiff, RK4, FoxGoodwin, LinearAcceleration,
       Newmark, HHT, WBZ, GeneralizedAlpha, MidPoint

# Structs - Noise estimation
export BayesianEst, GCVEst, LcurveEst, DerricoEst

# Structs - Noise denoising
export RegDenoising, KalmanDenoising

# Structs - Modal extraction
export BodeExtract, NyquistExtract

# Structs - Signal processing
export FFTParameters

# Functions
export excitation, solve, impulse_response, srs

export modefreq, modeshape, eigenmode, modal_matrices, modal_effective_mass,
       assembly, apply_bc, selection_matrix, rayleigh_damping_matrix, modal_damping_matrix

export agwn, acn, mult_noise, mixed_noise, varest, estimated_SNR, denoising

export c2d, ss_model, ss_modal_model, modal_parameters, c2r_modeshape

export freq_extraction, modeshape_extraction

export tfestimation, welch, spectrum, exponential, force, flattop, nutall,
       blackman_nutall, parzen, planck

export gradient, detrend

# Functions for plotting
export sv_plot, bode_plot, nyquist_plot, waterfall_plot

# Utility functions
export undefs

# Include files - Utils
include("utils/undefs.jl")

# Include files - Estimation
include("signal_processing/gradient.jl")
include("signal_processing/detrend.jl")
include("signal_processing/noise_estimation.jl")
include("signal_processing/denoising.jl")
include("signal_processing/modal_extraction.jl")
include("signal_processing/windows.jl")
include("signal_processing/signal.jl")

# Include files - Sdof
include("models/sdof_mdof.jl")
include("solvers/sdof_solvers.jl")

# Include files - Continuous structures
include("models/oned_structure.jl")
include("models/twod_structure.jl")

# Include files - Discrete structures
include("models/FEmodel.jl")
include("models/modal_model.jl")
include("solvers/direct_time_solvers.jl")
include("solvers/frequency_solvers.jl")
include("solvers/modal_time_solvers.jl")

# Include files - State space
include("models/state_space.jl")
include("solvers/state_space_solvers.jl")

# Include files - Excitation models
include("models/excitation.jl")

# Include files - Noise models
include("models/noise.jl")

# Include files - Visualization
include("utils/visualization.jl")

# Include files - Precompilation
# include("precompilation/precompilation.jl")
end