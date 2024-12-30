@setup_workload begin
    # Material properties
    E = 2.1e11
    ν = 0.33
    ρ = 7800.
    G = E/2(1 + ν)

    # Beam/bar dimensions
    L = 1.          # Length
    b = 3e-2        # Width
    h = 1e-2        # Height
    d = 1e-2        # Diameter
    S = b*h         # Section
    Iz = b*h^3/12.   # Second moment of area

    # Rod dimensions
    I₀ = π*d^4/32.  # Polar moment of inertia
    J = I₀          # Torsional constant

    # Plate dimensions
    a = 0.6         # Length
    b = 0.4         # Width
    hₚ = 5e-3       # Thickness

    # SDOF system
    m = 1.
    ω₀ = 2π*10.
    ξ = 0.01

    @compile_workload begin
        # Discrete model - SDOF
        sdof = SDOF(m, ω₀, ξ)

        # Continuous models - 1D
        bar = Bar(L, S, E, ρ)
        rod = Rod(L, I₀, J, G, ρ)
        beam = Beam(L, S, Iz, E, ρ)

        # Continuous models - 2D
        plate = Plate(a, b, hₚ, E, ρ, ν)

        # Eigenvalues and eigenvectors
        fmax = 3e3
        x = range(0., L, length = 3)

        ωbar, kbar = eigval(bar, fmax)
        ϕbar = eigmode(bar, kbar, x)

        ωrod, krod = eigval(rod, fmax)
        ϕrod = eigmode(rod, krod, x)

        ωbeam, kbeam =  eigval(beam, fmax)
        ϕbeam = eigmode(beam, kbeam, x)

        loc = [0.1, 0.2]
        ωplate, kplate = eigval(plate, fmax)
        ϕplate = eigmode(plate, kplate, loc[1], loc[2])

        # Excitations
        Δt = 1e-2
        t = 0:Δt:2Δt

        rect = Rectangle(1., 8e-3, 1e-2)
        F_rect = excitation(rect, t)

        hammer = Hammer(1., 8e-3, 9.7, 6e-4)
        F_hammer = excitation(hammer, t)

        triangle = Triangle(1., 8e-3, 5e-2)
        F_triangle = excitation(triangle, t)

        srect = SmoothRect(1., 8e-3, 5e-3, 5e-2)
        F_srect = excitation(srect, t)

        randexc = RandomExc(1., 8e-3, 5e-2, 0.1)
        F_rand = excitation(randexc, t)

        # Time solvers
        Kₙ, Mₙ, Cₙ = modal_model(ωplate, ξ)
        Fₙ = (F_hammer*ϕplate)'

        prob = LinearTimeProblem(Kₙ, Mₙ, Cₙ, Fₙ, t)
        CI = (D₀ = zeros(length(ωplate)), V₀ = zeros(length(ωplate)))

        solga = solve(prob, CI, GeneralizedAlpha())
        solcd = solve(prob, CI, CentralDiff())
        solhht = solve(prob, CI, HHT())
        solfg = solve(prob, CI, FoxGoodwin())
        solla = solve(prob, CI, LinearAcceleration())
        solnk = solve(prob, CI, Newmark())
        solwbz = solve(prob, CI, WBZ())
        solmp = solve(prob, CI, MidPoint())
        solrk = solve(prob, CI, RK4())

        # Frequency solvers
        freq = 70:71
        prob_freq = ModalFRF(ωplate, ξ, ϕplate, ϕplate, freq)
        solfrf = frf(prob_freq, :acc)
    end
end