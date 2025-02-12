@setup_workload begin
    # Excitation
    Δt = 1e-3
    t = 0.:Δt:0.035
    F₀ = 1.
    tstart = Δt
    duration = 0.03
    include("precompile_excitation.jl")

    # Sdof
    m = 1.
    f₀ = 10.
    ξ = 0.01
    freqs = 1.:2.:20.
    include("precompile_sdof.jl")

    # OneD structure
    L = 1.
    b = 3e-2
    h = 1e-2
    S = b*h
    Iz = b*h^3/12.
    J = Iz
    E = 2.1e11
    ν = 0.33
    G = E/(1 - 2*ν)
    ρ = 7800.
    fmax = freqs[end]
    x = [0.1, 0.9]
    include("precompile_oned_structure.jl")

    # TwoD structure
    Lp = 0.6
    bp = 0.4
    hp = 1e-3
    xp = [0.1, 0.5]
    yp = [0.1, 0.3]
    include("precompile_twod_structure.jl")

    # Mdof
    k_mdof = [1., 1.]
    m_mdof = ones(3)
    c_mdof = [0.1, 0.1]
    include("precompile_mdof.jl")

    # Modal Time solvers
    m_ss = Diagonal([2., 1.])
    k_ss = [6. -2.; -2. 4.]
    c_ss = [0.67 -0.11; -0.11 0.39]
    sine_exc = SineWave(F₀, tstart, duration, 10.)
    Fexc = excitation(sine_exc, t)
    include("precompile_modal_time_solvers.jl")

    # Direct Time solvers
    include("precompile_direct_time_solvers.jl")

    # State space
    include("precompile_state_space.jl")

    # State space solvers
    css = ss_model(k_ss, m_ss, c_ss)
    include("precompile_state_space_solvers.jl")

    # FEmodel
    include("precompile_FEmodel.jl")

    # Noise model
    x = sin.(2π*10*t)
    fs = 1/Δt
    snr_dB = 25.
    include("precompile_noise.jl")

    # Signal processing - Noise estimation and denoising
    y = agwn(x, snr_dB)
    include("precompile_noise_estimation.jl")

    # Signal processing - Signal
    ts = 0:Δt:3-Δt
    sig = cos.(2π*10*ts)
    ref = randn(length(ts))
    include("precompile_signal.jl")
end