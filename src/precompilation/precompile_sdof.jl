sdof = Sdof(m, f₀, ξ)

prob_sdof_free = SdofFreeTimeProblem(sdof, [1., 0.], t)
prob_sdof_harmo = SdofHarmonicTimeProblem(sdof, [0., 0.], t, 1., 10.)
prob_sdof_forced = SdofForcedTimeProblem(sdof, [0., 0.], t, excitation(hammer, t))
prob_sdof_frf = SdofFRFProblem(sdof, freqs)
prob_sdof_freq = SdofFrequencyProblem(sdof, freqs, ones(length(freqs)))

#solvers
@compile_workload begin
    sol_free = solve(prob_sdof_free)
    sol_harmo = solve(prob_sdof_harmo)
    sol_forced = solve(prob_sdof_forced)
    sol_frf = solve(prob_sdof_frf)
    sol_freq = solve(prob_sdof_freq)
    h = impulse_response(sdof, t)
end