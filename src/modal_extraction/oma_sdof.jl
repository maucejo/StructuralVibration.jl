function modes_extraction(prob::OMAProblem, alg::FSDD; width::Int = 1, min_prom = 0., max_prom = Inf, pks_indices = Int[], min_mac = 0.9, avg_alg = :sv)

    # Extraction of full CSD spectrum and freq from problem
    (; fullspec, freq) = prob

    # Compute first the singular value and left singular vector w.r.t frequency
    no, ni = size(fullspec)[1:2]

    if ni ≠ no
        throw(ArgumentError("The number of outputs must be equal to the number of intputs"))
    end

    sv, U = compute_sv(fullspec)

    # Estimate the peaks properties in the singular value spectrum
    if isempty(pks_indices)
        # Find peaks in the FRF
        pks = findmaxima(log10.(sv), width)

        # Define a peak prominence threshold to filter spurious peaks
        pks = peakproms!(pks, min = min_prom, max = max_prom) |> peakwidths!
        pks_indices = pks.indices
    end

    # Compute modes
    npeak = length(pks_indices)
    ms = similar(U, no, npeak)
    poles = similar(U, npeak)
    for (p, idx_peak) in enumerate(pks_indices)
        # Extract mode shapes
        ms[:, p], id_bell = estimate_modeshape_fsdd(sv, U, idx_peak, min_mac, avg_alg)

        # Compute enhanced CSD
        eGyy = compute_ecsd(fullspec, ms[:, p], id_bell)
        ωb = 2π * freq[id_bell]

        # Compute Least-squares matrices
        A = [eGyy -(2ωb.*eGyy) -one.(eGyy)]
        b = -ωb.^2 .* eGyy

        coeff = A\b

        ωn = sqrt(coeff[1])
        Ωn = coeff[2]
        # ξn = sqrt(1. - (Ωn/ωn)^2)

        poles[p] = -sqrt(ωn^2 - Ωn^2) + 1im*Ωn
    end

    return poles, ms ./ maximum(abs, ms; dims = 1)
end

"""
    compute_sv(Gyy)

Compute the first singular value and corresponding left singular vector of the CSD matrix `Gyy` for each frequency line.

**Input**
- `Gyy::Array{ComplexF64, 3}`: CSD matrix of size (no, no, nf), where `no` is the number of outputs and `nf` is the number of frequency lines.

**Outputs**
- `sv::Vector{Float64}`: Vector of first singular values for each frequency line.
- `U::Array{ComplexF64, 2}`: Matrix of left singular vectors corresponding to the first singular value for each frequency line.
"""
function compute_sv(Gyy::Array{T, 3}) where T <: Complex
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

"""
    estimate_modeshape_fsdd(sv, U, idx_peak, min_mac, avg_alg)

Estimate the mode shape corresponding to a peak in the singular value spectrum using the FSDD method.

**Inputs**
- `sv::Vector{Float64}`: Vector of singular values.
- `U::Array{ComplexF64, 2}`: Matrix of left singular vectors
- `idx_peak::Int`: Index of the peak in the singular value spectrum.
- `min_mac::Float64`: Minimum MAC value to consider neighboring frequencies.
- `avg_alg::Symbol`: Averaging algorithm to use
    - `:sv`: Weighted averaging using singular values
    - `:lin`: Linear averaging

**Outputs**
- `ms_avg::Vector{ComplexF64}`: Estimated mode shape.
- `id_bell::Vector{Int}`: Indices of frequency lines used in the
    estimation.
"""
function estimate_modeshape_fsdd(sv, U, idx_peak, min_mac, avg_alg)
    nf = length(sv)
    ms_ref = U[:, idx_peak]

    # Search for all the candidates in a frequency band defined by MAC values
    id_bell = [idx_peak]
    if min_mac == 1.
        return ms_ref, id_bell
    end

    

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
    if avg_alg == :sv
        ms_avg = (U[:, id_bell]*sv[id_bell])/sum(sv[id_bell])
    else
        ms_avg = mean(U[:, id_bell], dims = 2)
    end

    return ms_avg, id_bell
end

"""
    compute_ecsd(Gyy, ms, id_bell)

Compute the enhanced CSD spectrum using the estimated mode shape.

**Inputs**
- `Gyy::Array{ComplexF64, 3}`: CSD matrix of size (no, no, nf), where `no` is the number of outputs and `nf` is the number of frequency lines.
- `ms::Vector{ComplexF64}`: Estimated mode shape.
- `id_bell::Vector{Int}`: Indices of frequency lines used in the estimation.

**Output**
- `eGyy::Vector{ComplexF64}`: Enhanced CSD spectrum at the
    frequency lines defined by `id_bell`.
"""
function compute_ecsd(Gyy, ms, id_bell)

    nf = length(id_bell)
    eGyy = similar(Gyy, nf)
    for (f, idb) in enumerate(id_bell)
        eGyy[f] = ms'Gyy[:, :, idb]*ms
    end

    return eGyy
end