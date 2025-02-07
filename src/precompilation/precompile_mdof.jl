mdof = Mdof(k_mdof, m_mdof, c_mdof)

@compile_workload begin
    mdof_mesh = MdofMesh(mdof, 2)
    K_mdof, M_dof, C_dof = assembly(mdof)
end
