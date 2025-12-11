function modes_extraction(prob::OMAProblem, alg::FSDD; width::Int = 1, min_prom = 0., max_prom = Inf, pks_indices = Int[], min_mac = 1.)

    # Extraction of full CSD spectrum and freq from problem
    (; fullspec, freq)

    # Compute first the singular value and left singular vector w.r.t frequency
    no, ni = size(fullspec)[1:2]

    if ni ≠ no
        throw(ArgumentError("The number of outputs must be equal to the number of intputs"))
    end
    no = size(Gyy, 1)

    sv, U = compute_sv(Gyy)

    # Estimate the peaks properties in the singular value spectrum
    if isempty(pks_indices)
        # Find peaks in the FRF
        pks = findmaxima(log10.(sv), width)

        # Define a peak prominence threshold to filter spurious peaks
        pks = peakproms!(pks, min = min_prom, max = max_prom) |> peakwidths!
        pks_indices = pks.indices
    end

    # Compute modes
    ms = similar(U, no, length(pks_indices))
    for (p, idx_peak) in enumerate(pks_indices)
        # Extract mode shapes
        ms[:, p], id_bell = estimate_modeshape_fsdd(sv, U, idx_peak, min_mac)

        # Compute enhanced CSD
        eGyy = compute_ecsd(Gyy, ms[:, p], id_bell)

        # Compute Least-squares matrices
        
    end
end

function compute_sv(Gyy)
    no, _, nf = size(Gyy)

    sv = zeros(nf)
    U = similar(Gyy, no, nf)
    for f in 1:nf
        F = svd(Gyy[:, :, f])
        sv[f] = F.S[1]
        U[:, f] .= F.U[:, 1]
    end

    return sv, U
end

function estimate_modeshape_fsdd(sv, U, idx_peak, min_mac)
    nf = length(sv)
    ms_ref = U[:, idx_peak]

    if min_mac == 1.
        return ms_ref
    end

    # Search for all the candidates in a frequency band defined by MAC values
    id_bell = [idx_peak]

    # Search left
    lp = idx_peak - 1
    mac_val = 1.
    while lp ≥ 1 && mac_val ≥ min_mac
        # Compute the MAC value
        mac_val = mac(U[:, lp], ms_ref)
        push!(id_bell, lp)
        lp -= 1
    end

    # Search right
    rp = idx_peak + 1
    mac_val = 1.
    while rp ≤ nf && mac_val ≥ min_mac
        # Compute the MAC value
        mac_val = mac(U[:, rp], ms_ref)
        push!(id_bell, rp)
        rp += 1
    end

    # Sort indices
    sort!(id_bell)

    # Estimate weighted averaged mode shapes
    ms_avg = (U[:, id_bell]*sv[id_bell])/sum(sv[id_bell])

    return ms_avg, id_bell
end

function compute_ecsd(Gyy, ms, id_bell)

    nf = length(id_bell)
    eGyy = similar(Gyy, nf)
    for (f, idb) in enumerate(id_bell)
        eGyy[f] = ms'Gyy[:, :, idb]*ms
    end

    return eGyy
end