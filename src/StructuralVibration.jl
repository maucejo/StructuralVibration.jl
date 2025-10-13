module StructuralVibration

using FastGaussQuadrature, FFTW, Interpolations, LinearAlgebra, Optim, Peaks, Polynomials, PrecompileTools, ProgressMeter, Random, SpecialFunctions, Statistics, ToeplitzMatrices

using DSP: conv, filt, remez, filtfilt, rms

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
       SdofFRFProblem, SdofFrequencyProblem, StateSpaceTimeProblem, StateSpaceFRFProblem, StateSpaceModalFRFProblem, StateSpaceFreqProblem, StateSpaceModalFreqProblem, FreeModalTimeProblem,HarmonicModalTimeProblem, ForcedModalTimeProblem, ModalFRFProblem, DirectFRFProblem, ModalFreqProblem, DirectFreqProblem, DirectTimeProblem

# Structs - Time solvers
export CentralDiff, RK4, FoxGoodwin, LinearAcceleration,
       Newmark, HHT, WBZ, GeneralizedAlpha, MidPoint

# Structs - Noise estimation
export BayesEst, GCVEst, LCurveEst, DerricoEst

# Structs - Noise denoising
export BayesDenoising, GCVDenoising, LCurveDenoising, KalmanDenoising

# Structs - Modal extraction
export PeakPicking, CircleFit, LSFit, AutoEMASdofProblem, EMASdofSolution

# Structs - Signal processing
export FFTParameters

# Functions
export excitation, solve, impulse_response, srs

export modefreq, modeshape, eigenmode, modal_matrices, modal_effective_mass,
       assembly, apply_bc, selection_matrix, rayleigh_damping_matrix, modal_damping_matrix

export agwn, acn, mgwn, mixed_noise, varest, estimated_SNR, denoising,
       anti_alias

export c2d, ss_model, ss_modal_model, modal_parameters, c2r_modeshape

export freq_extraction, modeshape_extraction

export mof, mov, mpc, mcf, mpd, cmif, psif, msf, comac, ecomac, mac, frac, modal2poles, poles2modal, impulse_response, lsce, lscf, plscf

export tfestimate, welch, spectrum, rect, hanning, hamming, tukey, cosine,
       lanczos, triang, bartlett, gaussian, bartlett_hann, blackman, kaiser, dpss, exponential, force, flattop, nuttall, blackman_nuttall, blackman_harris, parzen, planck

export gradient, detrend

# Functions for plotting
export sv_plot, bode_plot, nyquist_plot, waterfall_plot, theme_choice

# Include files - Utils
include("utils/macro_utils.jl")

# Include files - Excitation models
include("models/excitation.jl")

# Include files - Noise models
include("models/noise.jl")

# Include files - Estimation
include("signal_processing/gradient.jl")
include("signal_processing/detrend.jl")
include("signal_processing/noise_estimation.jl")
include("signal_processing/denoising.jl")
include("signal_processing/windows.jl")
include("signal_processing/signal.jl")

# Include files - Modal extraction
include("modal_extraction/ema_sdof.jl")
include("modal_extraction/ema_indicators.jl")
include("modal_extraction/ema_mdof_poles.jl")
include("modal_extraction/ema_mdof_utils.jl")

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

# Include files - Visualization
include("utils/visualization.jl")

# Include files - Precompilation
# include("precompilation/precompilation.jl")
end