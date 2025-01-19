module StructuralVibration

using Parameters, ProgressMeter, LinearAlgebra, Statistics,
      DSP, FFTW, Interpolations, PrecompileTools

# # Structs - Models
export Sdof, Bar, Rod, Beam, Plate,
       ContinuousStateSpace, DiscreteStateSpace

# Structs - FE and discrete models
export Mesh

# Structs - Excitations
export Rectangle, Triangle, Hammer, SmoothRect,
       SineWave, SweptSine, GaussianPulse, ColoredNoise

# # Structs - Problems
export SdofFreeTimeProblem, SdofHarmonicTimeProblem, SdofForcedTimeProblem,
       SdofFRFProblem, SdofFrequencyProblem, StateSpaceTimeProblem, StateSpaceFRFProblem, StateSpaceFreqProblem,
       FreeModalTimeProblem, ForcedModalTimeProblem, ModalFRFProblem, DirectFRFProblem, ModalFreqProblem, DirectFreqProblem, DiscreteTimeProblem


# Structs - Time solvers
export CentralDiff, RK4, FoxGoodwin, LinearAcceleration,
       Newmark, HHT, WBZ, GeneralizedAlpha, MidPoint

# # Functions
export excitation, modefreq, modeshape, eigenmode, modal_matrices, solve,
       assembly, selection_matrix, agwn, acn, mult_noise, mix_noise, varest, estimated_SNR, c2d, ss_model

# Include files - Models
include("models/sdof.jl")
include("models/oned_structure.jl")
include("models/plate.jl")
include("models/modal_model.jl")
include("models/state_space.jl")
include("models/FEmodel.jl")
include("models/excitation.jl")
include("models/noise.jl")

# Include files - Solvers
include("solvers/sdof_solvers.jl")
include("solvers/state_space_solvers.jl")
include("solvers/frequency_solvers.jl")
include("solvers/direct_time_solvers.jl")
include("solvers/modal_time_solvers.jl")

# Include files - Utils
include("utils/calculus.jl")


# Include files - Precompilation
# include("precompilation/precompilation.jl")
end