@compile_workload begin
    prob = DirectTimeProblem(k_ss, m_ss, [zeros(length(t)) Fexc]', zeros(2, 2), (zeros(2), zeros(2)), t)

    sol_cd = solve(prob, CentralDiff(), progress = false)
    sol_rk4 = solve(prob, RK4(), progress = false)
    sol_newmarkfamily = solve(prob, progress = false)
end