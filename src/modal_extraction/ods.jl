"""
    ods(freq, y::Matrix{T}, pks_indices)

Compute the operational deflection shapes (ODS) from response spectrum.

**Inputs**
- `freq`: Frequency vector (Hz)
- `y`: Response spectrum matrix (no measurement points x no frequency lines)
- `pks_indices`: Indices of the peaks in the frequency vector

**Outputs**
- `freq_pks`: Frequencies at the peaks (Hz)
- `ods_amplitude`: ODS amplitudes (normalized to unit amplitude)
- `ods_phase`: ODS phases (radians)
"""
function ods(freq, y::Matrix{T}, pks_indices) where T<:Complex

    ods_amplitude = abs.(y[:, pks_indices])

    # Both approches are equivalent

    # For lightly damped structures, the response at resonance is approximately 90° out of phase with the excitation (i.e. the response phasor is mostly quadrature). Taking the imaginary part (or equivalently rotating the complex vector by ±90° and taking the real part) removes the arbitrary global phase (up to a sign) if the response is mostly quadrature.

    # ods_sign = sign.(imag(y[:, pks_indices]))
    # ods_mode = (ods_amplitude .* ods_sign) ./ maximum(ods_amplitude; dims = 1)

    # ods_phase = angle.(y[:, pks_indices]) .+ π/2
    # ods_mode = real(ods_amplitude .* cis.(ods_phase) ./ maximum(ods_amplitude; dims = 1))

    ods_phase = angle.(y[:, pks_indices])
    ods_mode = c2r_modeshape(ods_amplitude .* cis.(ods_phase) ./ maximum(ods_amplitude; dims = 1))

    return freq[pks_indices], ods_mode
end