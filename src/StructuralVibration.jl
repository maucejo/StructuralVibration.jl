module StructuralVibration

using Parameters, ProgressMeter, LinearAlgebra, Statistics,
      DSP, FFTW, Interpolations, Optim, SpecialFunctions, PrecompileTools

# Structs - Models
export Sdof, Bar, Rod, Strings, Beam, Plate, RectMembrane,
       ContinuousStateSpace, DiscreteStateSpace

# Structs - FE and discrete models
export Mesh

# Structs - Excitations
export Rectangle, Triangle, Hammer, SmoothRect,
       SineWave, SweptSine, GaussianPulse, ColoredNoise

# Structs - Problems
export SdofFreeTimeProblem, SdofHarmonicTimeProblem, SdofForcedTimeProblem,
       SdofFRFProblem, SdofFrequencyProblem, StateSpaceTimeProblem, StateSpaceFRFProblem, StateSpaceFreqProblem,
       FreeModalTimeProblem, ForcedModalTimeProblem, ModalFRFProblem, DirectFRFProblem, ModalFreqProblem, DirectFreqProblem, DirectTimeProblem

# Structs - Time solvers
export CentralDiff, RK4, FoxGoodwin, LinearAcceleration,
       Newmark, HHT, WBZ, GeneralizedAlpha, MidPoint

# Structs - Noise estimation
export BayesianEst, GCVEst, LcurveEst, DerricoEst

# Functions
export excitation, modefreq, modeshape, eigenmode, modal_matrices, solve,
       assembly, selection_matrix, agwn, acn, mult_noise, mix_noise, varest, estimated_SNR, c2d, ss_model

# Functions for plotting
export plot, bode_plot, nyquist_plot, waterfall_plot

# Include files - Sdof
include("models/sdof.jl")
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

# Include files - Utils
include("utils/calculus.jl")

# Include files - Visualization
include("utils/visualization.jl")

# Include files - Precompilation
# include("precompilation/precompilation.jl")
end