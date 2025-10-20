"""
    residues(frf, freq, poles; frange = [freq[1], freq[end]], type = :dis)

Compute residues and lower and upper residuals of a frequency response function (FRF) given its poles.

**Inputs**
- `frf`: Frequency response function (vector)
- `freq`: Frequency vector (Hz)
- `poles`: Vector of complex poles
**Keyword Arguments**
- `frange`: Frequency range for residue calculation (default: full range)
- `type`: Type of FRF
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

**Outputs**
- `res`: Residues corresponding to each pole
- `lr`: Lower residual
- `ur`: Upper residual
"""
function residues(frf, freq, poles; frange = [freq[1], freq[end]], type = :dis)
    # FRF post-processing - Frequency range reduction
    fidx = @. frange[1] ≤ freq ≤ frange[2]
    frf_red = frf[:, :, fidx]
    freq_red = freq[fidx]

    
    return res, lr, ur
end