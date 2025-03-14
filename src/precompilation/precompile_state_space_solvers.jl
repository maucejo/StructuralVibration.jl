@compile_workload begin
    prob_ss_time = StateSpaceTimeProblem(css, zeros(4), t, [zeros(length(t))  Fexc]')
    prob_ss_frf = StateSpaceFRFProblem(css, freqs)
    prob_ss_modalfrf = StateSpaceModalFRFProblem(css, freqs)
    prob_ss_freq = StateSpaceFreqProblem(css, freqs, [zeros(length(freqs)) ones(length(freqs))]')
    prob_ss_modalfreq = StateSpaceModalFreqProblem(css, freqs, [zeros(length(freqs)) ones(length(freqs))]')

    sol_ss_time = solve(prob_ss_time, :zoh, progress = false)
    sol_ss_time_rk4 = solve(prob_ss_time, RK4(), progress = false)
    sol_ss_frf = solve(prob_ss_frf, progress = false)
    sol_ss_modalfrf = solve(prob_ss_modalfrf, progress = false)
    sol_ss_freq = solve(prob_ss_freq, progress = false)
    sol_ss_modalfreq = solve(prob_ss_modalfreq, progress = false)
end