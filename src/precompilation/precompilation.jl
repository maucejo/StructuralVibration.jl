@setup_workload begin
    # Excitation
    Δt = 1e-3
    t = 0.:Δt:10.
    F₀ = 1.
    tstart = 5Δt
    duration = 5.
    include("precompile_excitation.jl")

    # Sdof
    m = 1.
    f₀ = 10.
    ξ = 0.01
    freqs = 1.:0.01:500.
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

    # State space
    m_ss = Diagonal([2., 1.])
    k_ss = [6. -2.; -2. 4.]
    c_ss = [0.67 -0.11; -0.11 0.39]
    include("precompile_state_space.jl")

    # FEmodel
    include("precompile_FEmodel.jl")
end