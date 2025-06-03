@compile_workload begin
    y = agwn(x, snr_dB)
    y .= acn(x, snr_dB, fs, band_freq = [100., 300.])
    y .= mgwn(x, snr_dB)
    y .= mixed_noise(x, snr_dB)
end