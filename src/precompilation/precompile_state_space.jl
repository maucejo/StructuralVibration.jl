@compile_workload begin
    css = ss_model(k_ss, m_ss, c_ss)
    dss = c2d(css, 1e-3, :zoh)

    λss, Ψss = eigenmode(css.Ac)
    ωss, ξss = modal_parameters(λss)

    ωₙ, ϕₙ = eigenmode(k_ss, m_ss)
    css_modal = ss_modal_model(ωₙ, 0.01, ϕₙ)
end