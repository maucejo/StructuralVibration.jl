module VibrationData

using Parameters, ProgressMeter, LinearAlgebra,
      DSP, Interpolations, PrecompileTools

# Structs - Models
export Plate, Bar, Rod, Beam, SDOF

# Structs - Excitations
export Rectangle, Triangle, RandomExc, Hammer, SmoothRect

# Structs - Problems
export LinearTimeProblem, ModalFRF, DirectFRF

# Functions
export excitation, eigval, eigmode, modal_model, solve, frf

# Time solvers
export CentralDiff, RK4, FoxGoodwin, LinearAcceleration,
       Newmark, HHT, WBZ, GeneralizedAlpha, MidPoint

# Noise utils
export agwn, varest, estimated_SNR

# Include files - Models
include("models/sdof.jl")
include("models/plate.jl")
include("models/bar_rod.jl")
include("models/beam.jl")
include("models/model.jl")

# Include files - Solvers
include("solvers/time_solvers.jl")
include("solvers/frequency_solvers.jl")

# Include files - Utils
include("utils/excitation.jl")
include("utils/noise.jl")

# Include files - Precompilation
# include("precompilation/precompilation.jl")
end