@compile_workload begin
    oned_mesh = OneDMesh(beam, 0., 5)
    Kfe, Mfe = assembly(beam, oned_mesh)
    Kbc = apply_bc(Kfe, oned_mesh)
    Mbc = apply_bc(Mfe, oned_mesh)

    ωfe, Φfe = eigenmode(Kbc, Mbc)

    Cray = rayleigh_damping_matrix(Kbc, Mbc, 1., 1.)
    Cmodal = modal_damping_matrix(Mbc, ωfe, 0.01ones(length(ωfe)), Φfe)
end