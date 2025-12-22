"""
    modes_extraction(prob, alg::FSDD; width, min_prom, max_prom , pks_indices, min_mac, avg_alg)

Extract modes from an Operational Modal Analysis problem using the Full Spectrum Decomposition and Decay (FSDD) method

**Inputs**
- `prob::OMAProblem`: OMA problem containing the CSD matrix and frequency lines.
- `alg::SdofOMA`: OMA method to use for pole extraction
    - `FSDD`: Frequency-spatial Domain Decomposition)
- `width::Int = 1`: Width parameter for peak detection in the singular value spectrum
- `min_prom::Float64 = 0.`: Minimum peak prominence for peak detection
- `max_prom::Float64 = Inf`: Maximum peak prominence for peak detection
- `pks_indices::Vector{Int} = Int[]`: Predefined indices of peaks
- `min_mac::Float64 = 0.9`: Minimum MAC value to consider neighboring frequencies for mode shape estimation
- `avg_alg::Symbol = :sv`: Averaging algorithm to use for mode shape estimation
    - `:sv`: Weighted averaging using singular values
    - `:lin`: Linear averaging

**Outputs**
- `poles::Vector{ComplexF64}`: Extracted poles
- `ms::Array{ComplexF64, 2}`: Extracted mode shapes
"""
function modes_extraction(prob::OMAProblem, alg::FSDD; width::Int = 1, min_prom = 0., max_prom = Inf, pks_indices = Int[], min_mac = 0.9, avg_alg = :sv)

    # Extraction of full CSD spectrum and freq from problem
    (; fullspec, freq) = prob

    # Compute first the singular value and left singular vector w.r.t frequency
    no, ni = size(fullspec)[1:2]

    if no < ni
        Gyy = permutedims(fullspec, (2, 1, 3))
    else
        Gyy = fullspec
    end
    no, ni, nf = size(Gyy)

    sv, U, V = compute_sv(Gyy)

    # Estimate the peaks properties in the singular value spectrum
    if isempty(pks_indices)
        # Find peaks in the FRF
        pks = findmaxima(log10.(sv), width)

        # Define a peak prominence threshold to filter spurious peaks
        pks = peakproms!(pks, min = min_prom, max = max_prom) |> peakwidths!
        pks_indices = pks.indices
    end

    # Initialization for EFDD if needed
    npoints = 2(nf - 1)
    dt = 1/2.56freq[end]
    t_sdof = 0:dt:(npoints÷2 - 1)*dt
    sdof_bell = zeros(eltype(fullspec), nf)
    corr_sdof = similar(freq, npoints÷2)

    # Compute modes
    npeak = length(pks_indices)
    ms = similar(U, no, npeak)
    poles = similar(U, npeak)
    for (p, idx_peak) in enumerate(pks_indices)
        # Extract mode shapes
        up, vp, id_bell = estimate_bell(sv, U, V, idx_peak, min_mac, avg_alg)

        ms[:, p] .= up

        # Compute enhanced CSD
        eGyy = compute_ecsd(Gyy, up, vp, id_bell)
        ωb = 2π * freq[id_bell]

        # Compute Least-squares matrices

        A = [eGyy -(2ωb.*eGyy) -one.(eGyy)]
        b = -ωb.^2 .* eGyy

        # Solve the system
        # Tips: Solving the complex system directly can lead to numerical issues
        # so we separate the real and imaginary parts to solve a real system of double size using
        Sa = [real(A); imag(A)]
        Sb = [real(b); imag(b)]
        res = (Sa'Sa)\(Sa'Sb)

        ωn = √res[1]
        Ωn = res[2]

        if ωn < Ωn
            @warn "Damping ratio is estimated to be greater than 100% for the mode at $(freq[idx_peak]) Hz. Using EFDD estimation of poles instead."
            sdof_bell[id_bell] .= eGyy

            corr_sdof .= irfft(sdof_bell, npoints)[1:npoints÷2]
            pks_sdof = findmaxima(corr_sdof)
            env = abs.(hilbert(corr_sdof))[pks_sdof.indices[1]:pks_sdof.indices[end]]
            t_env = t_sdof[pks_sdof.indices[1]:pks_sdof.indices[end]]

            A = [t_env one.(t_env)]
            b = log.(env)

            res = (A'A)\(A'b)

            poles[p] = res[1] + 1im*Ωn
        else
            poles[p] = -√(ωn^2 - Ωn^2) + 1im*Ωn
        end
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
    no, ni, nf = size(Gyy)

    sv = zeros(nf)
    U = similar(Gyy, no, nf)
    V = similar(Gyy, ni, nf)
    for f in 1:nf
        F = svd(Gyy[:, :, f])
        sv[f] = F.S[1]
        U[:, f] .= F.U[:, 1]
        V[:, f] .= F.V[:, 1]
    end

    return sv, U, V
end

"""
    estimate_bell(sv, U, idx_peak, min_mac, avg_alg)

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
function estimate_bell(sv, U, V, idx_peak, min_mac, avg_alg)
    nf = length(sv)
    u_ref = U[:, idx_peak]

    # Search for all the candidates in a frequency band defined by MAC values
    id_bell = [idx_peak]
    if min_mac == 1.
        return u_ref, id_bell
    end

    # Search left
    lp = idx_peak - 1
    mac_val = 1.
    while lp ≥ 1 && mac_val ≥ min_mac
        # Compute the MAC value
        mac_val = mac(U[:, lp], u_ref)
        push!(id_bell, lp)
        lp -= 1
    end

    # Search right
    rp = idx_peak + 1
    mac_val = 1.
    while rp ≤ nf && mac_val ≥ min_mac
        # Compute the MAC value
        mac_val = mac(U[:, rp], u_ref)
        push!(id_bell, rp)
        rp += 1
    end

    # Sort indices
    sort!(id_bell)

    # Estimate weighted averaged mode shapes
    if avg_alg == :sv
        u_avg = (U[:, id_bell]*sv[id_bell])/sum(sv[id_bell])
        v_avg = (V[:, id_bell]*sv[id_bell])/sum(sv[id_bell])
    else
        u_avg = mean(U[:, id_bell], dims = 2)
        v_avg = mean(V[:, id_bell], dims = 2)
    end

    return u_avg, v_avg, id_bell
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
function compute_ecsd(Gyy, u, v, id_bell)

    nf = length(id_bell)
    eGyy = similar(Gyy, nf)
    for (f, idb) in enumerate(id_bell)
        eGyy[f] = u'Gyy[:, :, idb]*v
    end

    return real(eGyy)
end