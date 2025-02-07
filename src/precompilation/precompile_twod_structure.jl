plate = Plate(Lp, bp, h, E, ρ, ν)
membrane = Membrane(Lp, bp, ρ*h, 1e3)

@compile_workload begin
    ωp, kp = modefreq(plate, fmax)
    ωm, km = modefreq(membrane, fmax)

    ϕp = modeshape(plate, kp, xp, yp)
    ϕm = modeshape(membrane, km, xp, yp)
end
