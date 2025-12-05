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
    real_normalization(Ψ)

Converts the complex modes to real modes

**Input**
* `Ψ`: Complex modes

**Output**
* `ϕn`: Real modes

**Reference**

[1] E. Hiremaglur. "Real-Normalization of Experimental Complex Modal Vectors with Modal Vector Contamination". MS Thesis. University of Cincinnati, 2014.
"""
function real_normalization(Ψ)

    # Initialization
    m = size(Ψ, 1)
    ϕ = similar(real(Ψ))
    x = similar(ϕ, m)
    y = similar(x)
    p = similar(x, 2)

    # Real mode shape calculation
    for (i, Ψi) in enumerate(eachcol(Ψ))
        x .= real(Ψi)
        y .= imag(Ψi)

        # Fit a first order line to the data
        p .= polyfit(x, y, 1)

        # Angle of maximum correlation line
        θ = atan(p[1])

        ϕ[:, i] .= real(Ψi*cis(-θ))
    end

    return ϕ
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
    p = modal2poles(fn[fidx], ξn[fidx])

    poles = length(p) ≤ order ? p : p[1:order]

    # If Stability diagram
    return stabdiag ? [poles; fill(complex(NaN, NaN), order - length(poles))] : poles
end

"""
    peak_proms!(pks; min = 0., max = Inf)

Compute the prominence of peaks and filter them based on minimum and maximum prominence.

**Inputs**
- `pks`: A NamedTuple containing peak information with fields:
    - `indices`: Indices of the peaks
    - `heights`: Heights of the peaks
    - `data`: The original data array from which peaks were identified
- `min`: Minimum prominence threshold (default: 0.)
- `max`: Maximum prominence threshold (default: Inf)

**Output**
- `pks`: Filtered `pks` to include only peaks with prominence within the specified range. The structure will also include a new field `proms` containing the computed prominences.
"""
function peak_proms!(pks; min = 0., max = Inf)
    data = pks.data

    n = length(data)
    npreaks = length(pks.indices)
    prominences = similar(data, npreaks)

    for (i, (idx, peak_height)) in enumerate(zip(pks.indices, pks.heights))
        # Find left base
        left_min = peak_height
        for j in idx-1:-1:1
            left_min = minimum([left_min, data[j]])
            if j > 1 && data[j] < data[j-1] && data[j] < data[j+1]
                break
            end
        end

        # Find right base
        right_min = peak_height
        for j in idx+1:n
            right_min = minimum([right_min, data[j]])
            if j < n && data[j] < data[j-1] && data[j] < data[j+1]
                break
            end
        end

        # Prominence is the height above the highest of the two bases
        base = maximum([left_min, right_min])
        prominences[i] = peak_height - base
    end

    # Filter peaks by prominence
    valid_mask = (prominences .≥ min) .& (prominences .≤ max)
    # Update pks in place
    filtered_pks = (
        indices = pks.indices[valid_mask],
        heights = pks.heights[valid_mask],
        proms = prominences[valid_mask],
        data = data
    )

    return merge(pks, filtered_pks)
end

"""
    peak_widths(pks; relheight = 0.5, min = 0., max = Inf)

Compute the widths of peaks at a specified relative height and filter them based on minimum and maximum width.

**Inputs**
- `pks`: A NamedTuple containing peak information with fields:
    - `indices`: Indices of the peaks
    - `heights`: Heights of the peaks
    - `proms`: Prominences of the peaks
    - `data`: The original data array from which peaks were identified
- `relheight`: Relative height at which to measure the width (default: 0.5)
- `min`: Minimum width threshold (default: 0.)
- `max`: Maximum width threshold (default: Inf)

**Output**
- `pks`: Filtered `pks` to include only peaks with widths within the
specified range. The structure will also include new fields `widths` and `edges` containing the computed widths and their corresponding left and right edges.
"""
function peak_widths!(pks; relheight = 0.5, min = 0., max = Inf)
    data = pks.data
    proms = pks.proms
    n = length(data)
    npeaks = length(pks.indices)

    widths = similar(proms, npeaks)
    left_edges = similar(proms, npeaks)
    right_edges = similar(proms, npeaks)

    for (i, (idx, peak_height, prom)) in enumerate(zip(pks.indices, pks.heights, proms))
        # Reference height based on prominence and relative height
        ref_height = peak_height - prom * relheight

        # Find left edge (where signal crosses reference height)
        left_idx = idx
        for j in idx-1:-1:1
            if data[j] ≤ ref_height
                # Linear interpolation for sub-sample precision
                if j < idx - 1
                    # Interpolate between j and j+1
                    t = (ref_height - data[j]) / (data[j+1] - data[j])
                    left_idx = j + t
                else
                    left_idx = j
                end
                break
            end
            if j == 1
                left_idx = 1.
            end
        end

        # Find right edge (where signal crosses reference height)
        right_idx = idx
        for j in idx+1:n
            if data[j] ≤ ref_height
                # Linear interpolation for sub-sample precision
                if j > idx + 1
                    # Interpolate between j-1 and j
                    t = (ref_height - data[j]) / (data[j-1] - data[j])
                    right_idx = j - 1 + t
                else
                    right_idx = j
                end
                break
            end
            if j == n
                right_idx = float(n)
            end
        end

        widths[i] = right_idx - left_idx
        left_edges[i] = left_idx
        right_edges[i] = right_idx
    end

    # Filter peaks by width
    valid_mask = (widths .≥ min) .& (widths .≤ max)

    filtered_edges = [(left_edges[i], right_edges[i]) for i in findall(valid_mask)]

    # Update pks in place
    filtered_pks = (
        indices = pks.indices[valid_mask],
        heights = pks.heights[valid_mask],
        proms = pks.proms[valid_mask],
        widths = widths[valid_mask],
        edges = filtered_edges,
        data = data
    )

    return merge(pks, filtered_pks)
end