bar = Bar(L, S, E, ρ)
rod = Rod(L, Iz, J, G, ρ)
strings = Strings(L, S, 100., ρ)
beam = Beam(L, S, Iz, E, ρ)

@compile_workload begin
    ωbar, kbar = modefreq(bar, fmax)
    ωrod, krod = modefreq(rod, fmax)
    ωstrings, kstrings = modefreq(strings, fmax)
    ωbeam, kbeam = modefreq(beam, fmax)

    ϕbar = modeshape(bar, kbar, x)
    ϕrod = modeshape(rod, krod, x)
    ϕstrings = modeshape(strings, kstrings, x)
    ϕbeam = modeshape(beam, kbeam, x)
end