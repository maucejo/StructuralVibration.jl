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
    H_padded = zpad(H, N, :ifft)

    return irfft(H_padded, N, nd)
end

"""
    poles_validity(poles, freq, stabdiag)

Check the validity of extracted poles and retain only those within the specified frequency range. Poles that do not appear as pairs of complex conjugate numbers with a positive real part or that are purely real are suppressed. For complex conjugate poles, only the pole with a positive imaginary part is retained.

**Inputs**
- `raw_poles::Vector{Complex}`: Extracted poles
- `order::Int`: Model order of the system
- `freq::Vector{Float64}`: Frequency vector
- `stabdiag::Bool`: Whether to compute stabilization diagram

**Output**
- `poles::Vector{Complex}`: Validated poles within the frequency range
"""
function poles_validity(raw_poles, order, freq, stabdiag)
    # Poles that do not appear as pairs of complex conjugate numbers with a positive real part or that are purely real are suppressed. For complex conjugate poles, only the pole with a positive imaginary part is retained.
    valid_poles = raw_poles[@. imag(raw_poles) < 0. && real(raw_poles) ≤ 0. && !isreal(raw_poles)]
    p_valid = intersect(raw_poles, conj.(valid_poles))
    sort!(p_valid, by = abs)

    # To avoid modifying the stability of the dynamic system in case of frequency offset, we set real(poles_new) = real(poles_old). This means that dn_new = dn_old * fn_old / fn_new.
    fn, ξn = poles2modal(p_valid)
    if freq[1] != 0.
        @. ξn *= fn / (fn + freq[1])
        fn .+= freq[1]
    end

    # Keep only the poles within frange
    fidx = @. freq[1] ≤ fn ≤ freq[end]
    poles = modal2poles(fn[fidx], ξn[fidx])

    # If Stability diagram
    return stabdiag ? [poles; fill(complex(NaN, NaN), order - length(poles))] : poles
end
