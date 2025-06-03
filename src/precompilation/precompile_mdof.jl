mdof = Mdof(k_mdof, m_mdof, c_mdof)

@compile_workload begin
    mdof_mesh = MdofMesh(mdof)
    K_mdof, M_mdof, C_mdof = assembly(mdof)
end
