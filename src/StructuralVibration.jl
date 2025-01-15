module StructuralVibration

using Parameters, ProgressMeter, LinearAlgebra, Statistics,
      DSP, FFTW, Interpolations, PrecompileTools

# Structs - Models
export Sdof, Bar, Rod, Beam, Plate,
       ContinuousStateSpace, DiscreteStateSpace

# Structs - FE and discrete models
export Mesh

# Structs - Excitations
export Rectangle, Triangle, Hammer, SmoothRect,
       SineWave, SweptSine, GaussianPulse, ColoredNoise

# Structs - Problems
export SdofTimeProblem, SdofFrequencyProblem, FreeModalTimeProblem,
       ForcedModalTimeProblem, ModalFRFProblem, DirectFRFProblem, ModalFreqProblem, DirectFreqProblem, DiscreteTimeProblem

# Structs - Time solvers
export CentralDiff, RK4, FoxGoodwin, LinearAcceleration,
       Newmark, HHT, WBZ, GeneralizedAlpha, MidPoint

# Functions
export excitation, modefreq, modeshape, eigenmode, modal_matrices, solve,
       assembly, selection_matrix, agwn, varest, estimated_SNR, c2d, ss_model

# Include files - Solvers
include("solvers/direct_time_solvers.jl")
include("solvers/modal_time_solvers.jl")
include("solvers/frequency_solvers.jl")

# Include files - Utils
include("utils/calculus.jl")

# Include files - Models
include("models/excitation.jl")
include("models/noise.jl")
include("models/sdof.jl")
include("models/oned_structure.jl")
include("models/plate.jl")
include("models/modal_model.jl")
include("models/FEmodel.jl")
include("models/state_space.jl")

# Include files - Precompilation
# include("precompilation/precompilation.jl")
end