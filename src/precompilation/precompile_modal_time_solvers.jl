@compile_workload begin
    ωn, ϕn = eigenmode(k_ss, m_ss)
    Kn, Mn, Cn = modal_matrices(ωn, 0.01)
    meff = modal_effective_mass(m_ss, ϕn, ones(2))

    prob_modal_free = FreeModalTimeProblem(k_ss, m_ss, 0.01, ([1., 0.], zeros(2)), t)
    prob_modal_harmo = HarmonicModalTimeProblem(k_ss, m_ss, 0.01, (zeros(2), zeros(2)), t, [1e4, 0.], 10.)
    prob_modal_forced = ForcedModalTimeProblem(k_ss, m_ss, 0.01, (zeros(2), zeros(2)), t, [zeros(length(t)) Fexc]')

    sol_modal_free = solve(prob_modal_free)
    sol_modal_harmo = solve(prob_modal_harmo)
    sol_modal_forced = solve(prob_modal_forced)
    sol_impulse = impulse_response(k_ss, m_ss, 0.01, t)
end