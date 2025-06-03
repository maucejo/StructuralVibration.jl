@compile_workload begin
    noisevar = varest(y)
    noisevar = varest(y, GCVEst())
    noisevar = varest(y, LcurveEst())
    noisevar = varest(y, DerricoEst())

    y_clean = denoising(y, GCVDenoising())
    y_clean .= denoising(y, LCurveDenoising())
    y_clean .= denoising(y, KalmanDenoising(rts = true))
end