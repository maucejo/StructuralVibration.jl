@compile_workload begin
    tf = tfestimation(ref, sig, 256, fs = fs, overlap = 0.5)
    pxx = welch(sig, 256, fs = fs, overlap = 0.5)
    spec = spectrum(sig, 256, fs = fs, overlap = 0.5)
end