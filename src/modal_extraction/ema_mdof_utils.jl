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
    if nd > 1
        np, nf = size(H)
        h = similar(freq, np, N)
    else
        nf = length(H)
        h = similar(freq, N)
    end

    if N > nf
        if nd == 1
            zero_pad = zeros(eltype(H), N - nf)
            H_padded = [H; zero_pad]
        else
            zero_pad = zeros(eltype(H), np, N - nf)
            H_padded = [H zero_pad]
        end
    else
        if nd == 1
            H_padded = H[1:N]
        else
            H_padded = H[:, 1:N]
        end
    end

    # IFFT - Enforce conjugate symmetry
    if N%2 == 0
        if nd == 1
            Hsym = [H_padded[1]; H_padded[2:Int(N/2 + 1) ]; conj.(reverse(H_padded[2:Int(N/2)]))]
        else
            Hsym = [H_padded[:, 1] H_padded[:, 2:Int(N/2 + 1)] conj.(reverse(H_padded[:, 2:Int(N/2)], dims = 2))]
        end
    else
        if nd == 1
            Hsym = [H_padded[1]; H_padded[2:Int((N+1)/2)]; conj.(reverse(H_padded[2:Int((N+1)/2)]))]
        else
            Hsym = [H_padded[:, 1] H_padded[:, 2:Int((N+1)/2)] conj.(reverse(H_padded[:, 2:Int((N+1)/2)], dims = 2))]
        end
    end

    return real(ifft(Hsym, nd))
end