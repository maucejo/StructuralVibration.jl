module StructuralVibration
       using DispatchDoctor
       @stable begin

       using FastGaussQuadrature, FFTW, Interpolations, LinearAlgebra, Optim,
             Peaks, Polynomials, ProgressMeter, Random, SpecialFunctions, Statistics, SkipNan, ToeplitzMatrices

       using DSP: conv, filt, filtfilt, remez, rms

       # Structs - Models
       export Bar, Beam, ContinuousStateSpace, DiscreteStateSpace, Mdof,
              Membrane, Plate, Rod, Sdof, Strings

       # Structs - FE and discrete models
       export MdofMesh, OneDMesh

       # Structs - Excitations
       export ColoredNoise, GaussianPulse, HalfSine, Hammer, HaverSine,
              Rectangle, SineWave, SmoothRect, SweptSine, Triangle

       # Structs - Problems
       export DirectFRFProblem, DirectFreqProblem, DirectTimeProblem,
              ForcedModalTimeProblem, FreeModalTimeProblem, HarmonicModalTimeProblem,ModalFRFProblem, ModalFreqProblem, SdofForcedTimeProblem, SdofFreeTimeProblem, SdofFrequencyProblem, SdofFRFProblem, SdofHarmonicTimeProblem, StateSpaceFreqProblem, StateSpaceFRFProblem, StateSpaceModalFRFProblem, StateSpaceModalFreqProblem, StateSpaceTimeProblem

       # Structs - Time solvers
       export CentralDiff, FoxGoodwin, GeneralizedAlpha, HHT,
              LinearAcceleration, MidPoint, Newmark, RK4, WBZ

       # Structs - Noise estimation
       export BayesianEst, DerricoEst, GCVEst, LCurveEst

       # Structs - Noise denoising
       export BayesDenoising, GCVDenoising, KalmanDenoising, LCurveDenoising

       # Structs - Modal extraction
       export AutoEMAMdofProblem, AutoEMASdofProblem, EMAProblem,
              EMASolution, OMAProblem, StabilizationAnalysis

       export CircleFit, CovSSI, DataSSI, LSCE, LSCF, LSFit, PeakPicking, pLSCF

       # Structs - Signal processing
       export FFTParameters

       # Functions
       export excitation, impulse_response, solve, srs

       export apply_bc, assembly, eigenmode, modefreq, modeshape,
              modal_damping_matrix, modal_effective_mass, modal_matrices, rayleigh_damping_matrix, selection_matrix

       export acn, agwn, anti_alias, denoising, estimated_SNR, mgwn,
              mixed_noise, varest

       export c2d, c2r_modeshape, modal_parameters, ss_model, ss_modal_model

       # Modal extraction
       export cmif, comac, compute_residuals, convert_csd, ecomac, frac,
              frf_reconstruction, half_csd_reconstruction, impulse_response, mac, mcf, modal2poles, mode_residues, mode2residues, modes_extraction, modeshape_extraction, mof, mpc, mpd, msf, mov, poles_extraction, poles2modal, psif, stabilization, csd_from_tf, half_csd, xcorr

       # Signal processing
       export bartlett, bartlett_hann, blackman, blackman_harris,
              blackman_nuttall, cosine, csd, dpss, exponential, flattop, flattri, force, gaussian, hamming, hanning, kaiser, lanczos, nuttall, parzen, planck, rect, spectrum, tfestimate, triang, tukey, welch

       export detrend, gradient

       # Functions for plotting
       export bode_plot, nyquist_plot, peaksplot, stabilization_plot, sv_plot,
              theme_choice, waterfall_plot

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
       include("modal_extraction/modal_analysis_types.jl")
       include("modal_extraction/ema_sdof.jl")
       include("modal_extraction/ema_indicators.jl")
       include("modal_extraction/ema_mdof.jl")
       include("modal_extraction/ema_frf.jl")
       include("modal_extraction/ema_utils.jl")
       include("modal_extraction/oma_utils.jl")
       include("modal_extraction/oma_mdof.jl")
       include("modal_extraction/oma_halfspec.jl")

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

end