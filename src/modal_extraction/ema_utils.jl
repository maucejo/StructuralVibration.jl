"""
    modal2poles(fn, ξn)

Convert natural frequencies and damping ratios to complex poles.

**Inputs**
- `fn`: Vector of natural frequencies (Hz)
- `ξn`: Vector of damping ratios

**Output**
- `poles`: Vector of complex poles
"""
function modal2poles(fn, ξn)
    ωn = 2π*fn

    return @. -ξn*ωn + 1im*ωn*√(1 - ξn^2)
end

"""
    poles2modal(poles)

Convert complex poles to natural frequencies and damping ratios.

**Input**
- `poles`: Vector of complex poles

**Outputs**
- `fn`: Vector of natural frequencies (Hz)
- `ξn`: Vector of damping ratios
"""
function poles2modal(poles)
    ωn = abs.(poles)
    fn = ωn/(2π)
    ξn = @. -real(poles)/ωn
    return fn, ξn
end

"""
    impulse_response(H, freq, fs)

Compute the impulse response from a frequency response function (FRF) using IFFT.

**Inputs**
- `H`: Frequency response function (can be a vector or a matrix with FRFs in rows)
- `freq`: Frequency vector (Hz)
- `fs`: Sampling frequency (Hz)

**Output**
- `h`: Impulse response

**Notes**
- The function enforces conjugate symmetry to ensure a real-valued impulse response.
"""
function impulse_response(H, freq, fs)
    # N points IFFT
    df = freq[2] - freq[1]
    N = Int(round(fs/df))

    # Data preparation
    nd = ndims(H)
    H_padded = zpad(H, N, :inv)

    return irfft(H_padded, N, nd)
end
