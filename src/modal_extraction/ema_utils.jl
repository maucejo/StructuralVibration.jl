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
    peak_proms!(pks; wlen = -1, min = 0., max = Inf)

Compute the prominence of peaks and filter them based on minimum and maximum prominence.

**Inputs**
- `pks`: A NamedTuple containing peak information with fields:
    - `indices`: Indices of the peaks
    - `heights`: Heights of the peaks
    - `data`: The original data array from which peaks were identified
- `wlen`: Window length in samples for evaluating prominence. If < 2, the entire signal is used (default: -1)
- `min`: Minimum prominence threshold (default: 0.)
- `max`: Maximum prominence threshold (default: Inf)
- `strict`: If true, only strict local maxima are considered peaks (default: false)

**Output**
- `pks`: Filtered `pks` to include only peaks with prominence within the specified range. New fields added:
    - `proms`: Prominences of the peaks
    - `left_bases`: Indices of the left bases of the peaks
    - `right_bases`: Indices of the right bases of the peaks

**Notes**
- This implementation follows the scipy.signal.peak_prominences algorithm, but implements additional heuristics to handle plateaus and saddle points in case of non-strict local maxima.
"""
function peak_proms!(pks; wlen = -1, min = 0., max = Inf, strict = false)
    data = pks.data
    n = length(data)
    npeaks = length(pks.indices)

    prominences = similar(data, npeaks)
    left_bases = similar(pks.indices, npeaks)
    right_bases = similar(pks.indices, npeaks)

    for (peak_nr, (peak_idx, peak_height)) in enumerate(zip(pks.indices, pks.heights))
        # Initialize window bounds
        i_min = 1
        i_max = n

        # Adjust window if wlen is specified
        if wlen >= 2
            # If wlen is even, the resulting window length is implicitly
            # rounded to next odd integer
            i_min = maximum([peak_idx - wlen ÷ 2, i_min])
            i_max = minimum([peak_idx + wlen ÷ 2, i_max])
        end

        # Check if peak is a strict local maximum
        peak_idx_adjusted = peak_idx
        is_strict_maximum = true

        if !strict
            # Check if left neighbor is equal to peak height
            if peak_idx > 1 && data[peak_idx - 1] ≥ peak_height * (1 - eps(eltype(data)) * 10)
                is_strict_maximum = false
            end

            # Check if right neighbor is equal to peak height
            if peak_idx < n && data[peak_idx + 1] ≥ peak_height * (1 - eps(eltype(data)) * 10)
                is_strict_maximum = false
            end

            # If not a strict maximum, find the plateau extent and use its center
            if !is_strict_maximum
                left_plateau = peak_idx
                right_plateau = peak_idx

                # Extend left while on plateau (within tolerance)
                while left_plateau > 1 && abs(data[left_plateau - 1] - peak_height) ≤ eps(eltype(data)) * maximum([abs(peak_height), 1.0]) * 10
                    left_plateau -= 1
                end

                # Extend right while on plateau (within tolerance)
                while right_plateau < n && abs(data[right_plateau + 1] - peak_height) ≤ eps(eltype(data)) * maximum([abs(peak_height), 1.0]) * 10
                    right_plateau += 1
                end

                # Use the center of the plateau
                peak_idx_adjusted = (left_plateau + right_plateau) ÷ 2
            end
        end

        # Find the left base in interval [i_min, peak]
        i = peak_idx_adjusted
        left_bases[peak_nr] = peak_idx_adjusted
        left_min = peak_height

        while i >= i_min && data[i] ≤ peak_height
            if data[i] < left_min
                left_min = data[i]
                left_bases[peak_nr] = i
            end
            i -= 1
        end

        # Find the right base in interval [peak, i_max]
        i = peak_idx_adjusted
        right_bases[peak_nr] = peak_idx_adjusted
        right_min = peak_height

        while i <= i_max && data[i] ≤ peak_height
            if data[i] < right_min
                right_min = data[i]
                right_bases[peak_nr] = i
            end
            i += 1
        end

        # Compute prominence
        base_height = maximum([left_min, right_min])
        prom = peak_height - base_height

        # Heuristic for zero or near-zero prominence (plateau or saddle point)
        if !strict && prom < eps(eltype(data)) * maximum([abs(peak_height), 1.0])
            # Try to find a reasonable descent in a neighborhood
            # Look for the minimum value in a window around the peak
            search_width = minimum([10, (i_max - i_min) ÷ 4])  # Adaptive search width

            left_search = maximum([peak_idx_adjusted - search_width, i_min])
            right_search = minimum([peak_idx_adjusted + search_width, i_max])

            local_min = minimum(data[left_search:right_search])

            # Estimate prominence based on local statistics
            if local_min < peak_height
                prom = peak_height - local_min
            else
                # Last resort: use a small fraction of the peak height
                # This represents the uncertainty in the peak location
                prom = 0.01peak_height  # 1% of peak height
            end

            # Update bases to reflect the search region
            left_bases[peak_nr] = left_search
            right_bases[peak_nr] = right_search
        end

        prominences[peak_nr] = prom
    end

    # Filter peaks by prominence
    valid_mask = (prominences .≥ min) .& (prominences .≤ max)

    # Update pks in place
    filtered_pks = (
        indices = pks.indices[valid_mask],
        heights = pks.heights[valid_mask],
        proms = prominences[valid_mask],
        left_bases = left_bases[valid_mask],
        right_bases = right_bases[valid_mask],
        data = data
    )

    return merge(pks, filtered_pks)
end

"""
    peak_widths!(pks; relheight = 0.5, min = 0., max = Inf)

Compute the widths of peaks at a specified relative height and filter them based on minimum and maximum width.

**Inputs**
- `pks`: A NamedTuple containing peak information with fields:
    - `indices`: Indices of the peaks
    - `heights`: Heights of the peaks
    - `proms`: Prominences of the peaks
    - `left_bases`: Left bases of the peaks (from `peak_proms!`)
    - `right_bases`: Right bases of the peaks (from `peak_proms!`)
    - `data`: The original data array from which peaks were identified
- `relheight`: Relative height at which to measure the width (default: 0.5)
- `min`: Minimum width threshold (default: 0.)
- `max`: Maximum width threshold (default: Inf)

**Output**
- `pks`: Filtered `pks` to include only peaks with widths within the
specified range.
New fields added:
    - `widths`: Widths of the peaks at the specified relative height
    - `width_heights`: Heights at which widths were measured
    - `edges`: Tuples of (left_intercept, right_intercept) for each peak

**Notes**
- This implementation follows the scipy.signal.peak_widths algorithm
- Requires `left_bases` and `right_bases` from `peak_proms!`
"""
function peak_widths!(pks; relheight = 0.5, min = 0., max = Inf)
    # Validation
    if relheight < 0
        throw(ArgumentError("`relheight` must be greater or equal to 0.0"))
    end

    data = pks.data
    proms = pks.proms
    left_bases = pks.left_bases
    right_bases = pks.right_bases
    n = length(data)
    npeaks = length(pks.indices)

    widths = similar(proms, npeaks)
    width_heights = similar(proms, npeaks)
    left_ips = similar(proms, npeaks)
    right_ips = similar(proms, npeaks)

    for (p, (peak_idx, peak_height, prom, i_min, i_max)) in enumerate(zip(pks.indices, pks.heights, proms, left_bases, right_bases))
        # Validate bounds and order
        if !(1 <= i_min <= peak_idx <= i_max <= n)
            throw(ArgumentError("prominence data is invalid for peak at index $peak_idx"))
        end

        # # Heuristic for zero or near-zero prominence peaks: half-power bandwidth method
        # if min_width_heuristic && prom < eps(eltype(data)) * maximum([abs(peak_height), 1.0])
        #     # Half-power height: peak_height / √2
        #     height = peak_height*(1.  - 1/√2)
        #     width_heights[p] = height

        #     # Find left intersection point where data crosses height
        #     i = peak_idx
        #     while data[i] > height && i > 1
        #         i -= 1
        #     end

        #     left_ip = float(i)
        #     # Interpolate for sub-sample precision
        #     if i < peak_idx && data[i] ≤ height && data[i + 1] > height
        #         left_ip += (height - data[i]) / (data[i + 1] - data[i])
        #     end

        #     # Find right intersection point where data crosses height
        #     i = peak_idx
        #     while data[i] > height && i < n
        #         i += 1
        #     end

        #     right_ip = float(i)
        #     # Interpolate for sub-sample precision
        #     if i > peak_idx && data[i] ≤ height && data[i - 1] > height
        #         right_ip -= (height - data[i]) / (data[i - 1] - data[i])
        #     end

        #     widths[p] = right_ip - left_ip
        #     left_ips[p] = left_ip
        #     right_ips[p] = right_ip
        #     continue
        # end

        # Height of the contour line at which width is measured
        height = peak_height - prom * relheight
        width_heights[p] = height

        # Find intersection point on left side
        # Start at peak and move left until we find where signal crosses reference height
        i = peak_idx
        while i_min < i && height < data[i]
            i -= 1
        end

        left_ip = float(i)
        if data[i] < height
            # Interpolate if true intersection height is between samples
            left_ip += (height - data[i]) / (data[i + 1] - data[i])
        end

        # Find intersection point on right side
        # Start at peak and move right until we find where signal crosses reference height
        i = peak_idx
        while i < i_max && height < data[i]
            i += 1
        end

        right_ip = float(i)
        if data[i] < height
            # Interpolate if true intersection height is between samples
            right_ip -= (height - data[i]) / (data[i - 1] - data[i])
        end

        widths[p] = right_ip - left_ip
        left_ips[p] = left_ip
        right_ips[p] = right_ip
    end

    # Filter peaks by width
    valid_mask = (widths .≥ min) .& (widths .≤ max)

    filtered_edges = [(left_ips[i], right_ips[i]) for i in findall(valid_mask)]

    # Update pks in place
    filtered_pks = (
        indices = pks.indices[valid_mask],
        heights = pks.heights[valid_mask],
        proms = pks.proms[valid_mask],
        left_bases = left_bases[valid_mask],
        right_bases = right_bases[valid_mask],
        widths = widths[valid_mask],
        width_heights = width_heights[valid_mask],
        edges = filtered_edges,
        data = data
    )

    return merge(pks, filtered_pks)
end