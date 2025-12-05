"""
    ods(y, freq, pks_indices)
    ods(Gyy, freq, pks_indices)

Compute the operational deflection shapes (ODS) from response spectrum.

**Inputs**
- `y::Matrix{Complex}`: Response spectrum matrix (n_output x n_frequency)
- `Gyy::Array{Complex, 3}`: Cross-spectral density matrix (n_ouput x n_input x n_frequency)
- `freq`: Frequency vector (Hz)
- `pks_indices`: Indices of the peaks in the frequency vector
- `min_mac` (optional): Minimum MAC value to consider a candidate ODS for averaging

**Outputs**
- `ods_mode`: Operational deflection shape
- `freq_pks`: Frequencies at the peaks (Hz)

**Note**

- The ODS are normalized such that the maximum absolute value of each mode shape is 1.
"""
function ods(y::Matrix{T}, freq, pks_indices) where {T <: Complex}

    # Both approches are equivalent
    # For lightly damped structures, the response at resonance is approximately 90° out of phase (in quadrature) with the excitation. Taking the imaginary part (or equivalently rotating the complex vector by ±90° and taking the real part) removes the arbitrary global phase (up to a sign) if the response is mostly quadrature.

    ods_mode = y[:, pks_indices]
    ods_mode ./= maximum(abs.(ods_mode); dims = 1)

    # ods_amplitude = abs.(y[:, pks_indices])
    # ods_phase = angle.(y[:, pks_indices]) .+ π/2
    # ods_mode = real(ods_amplitude .* cis.(ods_phase))

    # ods_sign = sign.(imag(y[:, pks_indices]))
    # ods_mode = (ods_amplitude .* ods_sign) ./ maximum(ods_amplitude; dims = 1)

    return ods_mode, freq[pks_indices]
end

function ods(y::Matrix{T}, freq, pks_indices, min_mac) where {T <: Complex}

    ny, nf = size(y)
    ods_mode_ref, freq_pks = ods(y, freq, pks_indices)

    if pks_indices isa Int
        pks_indices = [pks_indices]
    end

    ods_ref = similar(ods_mode_ref, ny)
    ods_avg = similar(ods_mode_ref)
    candidate = similar(ods_ref)
    ods_cand = zeros(eltype(ods_ref), ny)
    for (p, idx_peak) in enumerate(pks_indices)
        ods_ref .= ods_mode_ref[:, p]

        # Search for all the candidates in a frequency band defined by MAC values
        lp = idx_peak - 1
        mac_val = 1.
        navg = 1
        while lp ≥ 1 && mac_val ≥ min_mac
            # Define a candidate
            candidate .= ods(y, freq, lp)[1]

            # Compute the MAC value
            mac_val = mac(candidate, ods_ref)

            if mac_val ≥ min_mac
                navg += 1
                c = msf(candidate, ods_ref)
                # push!(ods_cand, candidate .* c)
                ods_cand .+= candidate .* c^2
            end
            lp -= 1
        end

        rp = idx_peak + 1
        mac_val = 1.
        while rp ≤ nf && mac_val ≥ min_mac
            # Define a candidate
            candidate .= ods(y, freq, rp)[1]

            # Compute the MAC value
            mac_val = mac(candidate, ods_ref)

            if mac_val ≥ min_mac
                navg += 1
                c = msf(candidate, ods_ref)
                ods_cand .+= candidate .* c
            end
            rp += 1
        end

        # Average all the candidates
        ods_avg[:, p] .= ods_cand/navg
        ods_avg[:, p] ./= maximum(abs.(ods_avg[:, p]))
        fill!(ods_cand, eltype(ods_ref)(0.))
    end

    return ods_avg, freq_pks
end

function ods(Gyy::Array{T, 3}, freq, pks_indices) where {T <: Complex}
    npeak = length(pks_indices)
    no, ni = size(Gyy)[1:2]

    if no ≤ ni
        Syy = permutedims(Gyy, (2, 1, 3))
    else
        Syy = Gyy
    end

    ods_complex = similar(Syy, max(no, ni), npeak)
    for (p, idx_peak) in enumerate(pks_indices)
        ods_complex[:, p] .= svd(Syy[:, :, idx_peak]).U[:, 1]
    end

    return ods_complex ./ maximum(abs.(ods_complex); dims = 1), freq[pks_indices]
end

function ods(Gyy::Array{T, 3}, freq, pks_indices, min_mac) where {T <: Complex}

    ods_mode_ref, freq_pks = ods(Gyy, freq, pks_indices)
    ny = size(ods_mode_ref, 1)
    nf = size(Gyy, 3)

    if pks_indices isa Int
        pks_indices = [pks_indices]
    end

    ods_ref = similar(ods_mode_ref, ny)
    ods_avg = similar(ods_mode_ref)
    candidate = similar(ods_ref)
    ods_cand = zeros(eltype(ods_ref), ny)
    for (p, idx_peak) in enumerate(pks_indices)
        ods_ref .= ods_mode_ref[:, p]

        # Search for all the candidates in a frequency band defined by MAC values
        navg = 1
        lp = idx_peak - 1
        mac_val = 1.
        while lp ≥ 1 && mac_val ≥ min_mac
            # Define a candidate
            candidate .= ods(Gyy, freq, lp)[1]

            # Compute the MAC value
            mac_val = mac(candidate, ods_ref)

            if mac_val ≥ min_mac
                navg += 1
                c = msf(candidate, ods_ref)
                ods_cand .+= candidate .* c
            end
            lp -= 1
        end

        rp = idx_peak + 1
        mac_val = 1.
        while rp ≤ nf && mac_val ≥ min_mac
            # Define a candidate
            candidate .= ods(Gyy, freq, rp)[1]

            # Compute the MAC value
            mac_val = mac(candidate, ods_ref)

            if mac_val ≥ min_mac
                navg += 1
                c = msf(candidate, ods_ref)
                ods_cand .+= candidate .* c
            end
            rp += 1
        end

        # Average all the candidates
        ods_avg[:, p] .= ods_cand/navg
        ods_avg[:, p] ./= maximum(abs.(ods_avg[:, p]))
        fill!(ods_cand, eltype(ods_ref)(0.))
    end

    return ods_avg, freq_pks

end