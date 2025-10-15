## Structure definition for method selection
struct LSCE end
struct LSCF end
struct PLSCF end

"""
    EMAMdofStabilization(poles, modefn, mode_stabfn, mode_stabdr)

Data structure to hold the results of the stabilization analysis.

**Fields**
- `poles::Vector{Vector{Complex}}`: Vector of vectors containing extracted poles at each model order
- `modefn::Matrix{Real}`: Matrix containing the natural frequencies (useful for plotting)
- `mode_stabfn::Matrix{Bool}`: Matrix indicating the stability of natural frequencies
- `mode_stabdr::Matrix{Bool}`: Matrix indicating the stability of damping ratios
"""
@show_data struct EMAMdofStabilization{T <: Complex}
    poles::Vector{Vector{T}}   # Extracted poles at each model order
    modefn::Matrix{T}          # Natural frequencies (used for plotting)
    mode_stabfn::Matrix{Bool}  # Stability of natural frequencies
    mode_stabdr::Matrix{Bool}  # Stability of damping ratios
end

## Function for modal extraction
"""
    lsce(frf, freq, order, fs; frange, stabdiag = false)

Perform Least Squares Complex Exponential (LSCE) method to extract complex poles from Frequency Response Function (FRF) data.

**Inputs**
- `frf`: 3D array of Frequency Response Functions (FRF) (array nm x ne x nf)
- `freq`: Vector of frequency values (Hz)
- `order::Int`: Model order (number of poles to extract)
- `frange`: Frequency range for analysis (default: [freq[1], freq[end]])
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)

**Output**
- `poles`: Vector of extracted complex poles
"""
function lsce(frf, freq, order::Int; frange = [freq[1], freq[end]], stabdiag = false)
    # FRF post-processing - Frequency range reduction
    fidx = @. frange[1] ≤ freq ≤ frange[2]
    frf_red = frf[:, :, fidx]
    freq_red = freq[fidx]
    fsred = 2.56*(freq_red[end]-freq_red[1]) # Approximate sampling frequency of the truncated signal

    # Impulse function calculation
    nm, ne, nf = size(frf_red)
    FRF = reshape(frf_red, (nm*ne, nf))
    Hk = impulse_response(FRF, freq_red, fsred)

    # Model order and number of samples
    n = 2order
    nsamples = 10n # Number of samples used for the identification

    # Preallocation
    H0 = similar(Hk, nsamples, n+1)
    A = similar(Hk, nsamples*nm*ne, n)
    b = similar(Hk, nsamples*nm*ne)
    for (k, hk) in enumerate(eachrow(Hk))
        # Hankel matrix construction
        H0 .= Hankel(hk[1:nsamples], hk[nsamples:nsamples+n])

        # Construct A and b
        A[(1:nsamples) .+ (k-1)*nsamples, :] .= H0[:, 1:end-1]
        b[(1:nsamples) .+ (k-1)*nsamples] .= -H0[:, end]
    end

    # Compute the coefficients of the characteristic polynomial
    α = [A\b; 1.]

    # Compute the poles
    V = roots(Polynomial(α))
    p = log.(V)*fsred

    # Poles that do not appear as pairs of complex conjugate numbers with a positive real part or that are purely real are suppressed. For complex conjugate poles, only the pole with a positive imaginary part is retained.
    valid_poles = p[@. imag(p) < 0. && real(p) ≤ 0. && !isreal(p)]
    p_valid = intersect(p, conj.(valid_poles))
    sort!(p_valid, by = abs)

    # To avoid modifying the stability of the dynamic system in case of frequency offset, we set real(poles_new) = real(poles_old). This means that dn_new = dn_old * fn_old / fn_new.
    fn, ξn = poles2modal(p_valid)
    if freq_red[1] != 0.
        @. ξn *= fn / (fn + freq_red[1])
        fn .+= freq_red[1]
    end

    # Keep only the poles within frange
    fidx = @. frange[1] ≤ fn ≤ frange[2]
    poles = modal2poles(fn[fidx], ξn[fidx])

    # If Stability diagram
    return stabdiag ? [poles; fill(NaN, Nmodes - length(poles))] : poles
end

"""
    lscf(frf, freq, order; frange, stabdiag)

Perform Least Squares Complex Frequency (LSCF) method to extract complex poles from Frequency Response Function (FRF) data.

**Inputs**
- `frf`: 3D array of Frequency Response Functions (FRF) (array nm x ne x nf)
- `freq`: Vector of frequency values (Hz)
- `order::Int`: Model order (number of poles to extract)
- `frange`: Frequency range for analysis (default: [freq[1], freq[end]])
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)
- `weighting`: Boolean to indicate if the weighting based on the variance of each FRF is applied (default: true)

**Output**
- `poles`: Vector of extracted complex poles

**Reference**

[1] El-Kafafy M., Guillaume P., Peeters B., Marra F., Coppotelli G. (2012).Advanced Frequency-Domain Modal Analysis for Dealing with Measurement Noise and Parameter Uncertainty. In: Allemang R., De Clerck J., Niezrecki C., Blough J. (eds) Topics in Modal Analysis I, Volume 5. Conference Proceedings of the Society for Experimental Mechanics Series. Springer, New York, NY
"""
function lscf(frf, freq, order::Int; frange = [freq[1], freq[end]], stabdiag = false, weighting = true)
    # FRF post-processing - Frequency range reduction
    fidx = @. frange[1] ≤ freq ≤ frange[2]
    frf_red = frf[:, :, fidx]
    freq_red = freq[fidx]

    # FRF post-processing - Reshape FRF matrix
    (nm, ne, nf) = size(frf_red)
    FRF = reshape(frf_red, (nm*ne, nf))  # each row = one FRF across frequencies

    # LSCF computation - stable numerics
    ω = 2(freq_red .- freq_red[1])                # Reduced angular frequency
    Δt = 1/(2(freq_red[end] - freq_red[1]))    # Sampling interval
    modelOrder = 0:order
    nmodel = order + 1

    # Build basis matrix X0 (nf x nmodel), complex
    # Use the same basis as original: exp.(-1im*ω*modelOrder'*Δt)
    X0 = cispi.(-ω*modelOrder'*Δt)  # nf x nmodel

    # Stable inversion of Rk = X0' * X0 using SVD of X0 (X0 = U * Σ * V')
    Ux, Σx, Vx = svd(real(X0'X0))
    # Filter small singular values for numerical stability
    tol = maximum(Σx) * eps(eltype(Σx)) * length(Σx)
    idx = Σx .> tol
    if !any(idx)
        error("LSCF: basis matrix X0 is rank-deficient (all singular values below tolerance).")
    end
    # Build filtered inverse of Rk = V * diag(S.^2) * V'
    Rk_inv = Vx[:, idx] * Diagonal(1.0 ./ (Σx[idx])) * Ux[:, idx]'

    # Preallocation
    M = similar(freq, nmodel, nmodel)
    Yk = similar(X0)
    Sk = similar(M)
    Tk = similar(M)

    # Process each FRF (each row of FRF is length nf)
    for Hk in eachrow(FRF)
        # ensure H is a column-like vector of length nf
        # Compute weighting scalar for this FRF (keep original behavior)
        vk = weighting ? only(varest(Hk)) : 1.

        @. Yk = -Hk*X0
        Sk .= real(X0'Yk)
        Tk .= real(Yk'Yk)

        # Accumulate robust contribution to M
        # Use real parts since final M is real symmetric
        M .+= (Tk .- Sk'*(Rk_inv*Sk))/vk
    end

    # Extract submatrices for denominator solve
    A = -M[1:order, 1:order]
    b = M[1:order, nmodel]

    # Solve A * α = b in stable way using SVD and filtering
    UA, ΣA, VA = svd(A)
    tolA = maximum(ΣA) * eps(eltype(ΣA)) * length(ΣA)
    idxA = ΣA .> tolA
    if !any(idxA)
        error("LSCF: matrix A is (nearly) singular.")
    end
    α_sub = VA[:, idxA] * Diagonal(1.0 ./ ΣA[idxA]) * (UA[:, idxA]'b)
    α = [α_sub; 1.] # characteristic polynomial coefficients

    # Compute poles from polynomial roots
    V = roots(Polynomial(α))
    p = -log.(V) / Δt

    # Poles filtering (keep negative real part and negative imag to pick one of conjugates)
    valid_poles = p[@. imag(p) < 0. && real(p) ≤ 0. && !isreal(p)]
    poles = intersect(p, conj.(valid_poles))
    sort!(poles, by = abs)

    # Frequency offset correction (restore original frequency offset)
    if freq_red[1] != 0.0 && !isempty(poles)
        fn, ξn = poles2modal(poles)
        poles .= modal2poles(fn .+ freq_red[1], @. ξn * fn / (fn + freq_red[1]))
    end

    # Return format based on stabdiag flag
    return stabdiag ? [poles; fill(NaN, Nmodes - length(poles))] : poles
end

"""
    plscf(frf, freq, order; frange, stabdiag)

Perform Polyreference Least Squares Complex Frequency (pÒLSCF) method to extract complex poles from Frequency Response Function (FRF) data.

**Inputs**
- `frf`: 3D array of Frequency Response Functions (FRF) (array nm x ne x nf)
- `freq`: Vector of frequency values (Hz)
- `order::Int`: Model order (number of poles to extract)
- `frange`: Frequency range for analysis (default: [freq[1], freq[end]])
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)
- `weighting`: Boolean to indicate if the weighting based on the variance of each FRF is applied (default: true)

**Output**
- `poles`: Vector of extracted complex poles

**Reference**

[1] El-Kafafy M., Guillaume P., Peeters B., Marra F., Coppotelli G. (2012).Advanced Frequency-Domain Modal Analysis for Dealing with Measurement Noise and Parameter Uncertainty. In: Allemang R., De Clerck J., Niezrecki C., Blough J. (eds) Topics in Modal Analysis I, Volume 5. Conference Proceedings of the Society for Experimental Mechanics Series. Springer, New York, NY
"""
function plscf(frf, freq, order::Int; frange = [freq[1], freq[end]], stabdiag = false, weighting = true)
    # FRF post-processing - Frequency range reduction
    fidx = @. frange[1] ≤ freq ≤ frange[2]
    frf_red = frf[:, :, fidx]
    freq_red = freq[fidx]

    # FRF post-processing - Reshape FRF matrix
    (nm, ne, nf) = size(frf_red)
    if ne > nm
        frf_red = permutedims(frf_red, (2, 1, 3))
        (nm, ne, nf) = size(frf_red)
    end

    # pLSCF computation
    ω = 2(freq_red .- freq_red[1]) # Reduced angular frequency
    Δt = 1/(2(freq_red[end] - freq_red[1])) # Sampling interval
    modelOrder = 0:order
    nmodel = order + 1

    # Reduced normal equation computation
    Ω0 = cispi.(-ω*modelOrder'*Δt)
    # Precompute R0
    X0 = similar(Ω0, nf, ne*nmodel)
    for f in 1:nf
        X0[f, :] .= kron(Ω0[f, :], I(ne))
    end
    # R0 = real(X0'X0)
    # Stable inversion of Rk = X0' * X0 using SVD of X0 (X0 = U * Σ * V')
    Ux, Σx, Vx = svd(real(X0'X0))
    # Filter small singular values for numerical stability
    tol = maximum(Σx)*eps(eltype(Σx))*length(Σx)
    idx = Σx .> tol
    if !any(idx)
        error("LSCF: basis matrix X0 is rank-deficient (all singular values below tolerance).")
    end
    # Build filtered inverse of R0
    R0_inv = Vx[:, idx] * Diagonal(1.0 ./ (Σx[idx])) * Ux[:, idx]'

    # Preallocation
    M = similar(frange, ne*nmodel, ne*nmodel)
    Y0 = similar(Ω0, nf, ne*nmodel)
    H0 = similar(Ω0, ne, nf)
    S0 = similar(M)
    T0 = similar(M)
    vk = ones(ne)

    for i in 1:nm
        H0 .= frf_red[i, :, :]
        if weighting
            vk .= varest(H0)
        end

        for f in 1:nf
            Y0[f, :] .= -kron(Ω0[f, :], H0[:, f])
        end

        # Build R0, S0, T0
        S0 .= real(X0'Y0)
        T0 .= real(Y0'Y0)

        # Accumulation of M
        # M .+= (T0 .- S0'*(R0\S0))./vk
        M .+= (T0 .- S0'*(R0_inv*S0))./vk
    end

    # Computation of the coefficients of the denominator
    A = -M[1:order*ne, 1:order*ne]
    b = M[1:order*ne, (order*ne + 1):nmodel*ne]

    # Solve A * α = b in stable way using SVD and filtering
    UA, ΣA, VA = svd(A)
    tolA = maximum(ΣA) * eps(eltype(ΣA)) * length(ΣA)
    idxA = ΣA .> tolA
    if !any(idxA)
        error("LSCF: matrix A is (nearly) singular.")
    end
    # α = A\b
    α = VA[:, idxA] * Diagonal(1.0 ./ ΣA[idxA]) * (UA[:, idxA]'b)

    # Construct the companion matrix
    Id = I(ne*(order - 1))
    Ze = zeros(ne*(order - 1), ne)
    C = vcat(hcat(Ze, Id), -α')

    # Pole calculation
    p = -log.(eigvals(C))/Δt

    # Poles that do not appear as pairs of complex conjugate numbers with a positive real part or that are purely real are suppressed. For complex conjugate poles, only the pole with a positive imaginary part is retained.
    valid_poles = p[@. imag(p) < 0. && real(p) ≤ 0. && !isreal(p)]
    poles = intersect(p, conj.(valid_poles))

    # To avoid modifying the stability of the dynamic system in case of frequency offset, we set real(poles_new) = real(poles_old). This means that dn_new = dn_old * fn_old / fn_new.
    if freq_red[1] != 0.
        fn, ξn = poles2modal(poles)
        poles .= modal2poles(fn .+ freq_red[1], @. ξn*fn/(fn + freq_red[1]))
    end
    sort!(poles, by = abs)

    # If Stability diagram
    return stabdiag ? [poles; fill(NaN, Nmodes - length(poles))] : poles
end

## Function for stabilization diagram analysis
"""
    stabilization(frf, freq, max_order::Int, method = LSCF(); frange = [freq[1], freq[end]], weighting = true, stabcrit = [0.01, 0.05])

Perform stabilization diagram analysis using the specified modal extraction method (LSCE, LSCF, or PLSCF).

**Inputs**
- `frf`: 3D array of Frequency Response Functions (FRF) (array nm x ne x nf)
- `freq`: Vector of frequency values (Hz)
- `max_order::Int`: Maximum model order for the stabilization analysis
- `method`: Modal extraction method to use (default: LSCF)
- `frange`: Frequency range for analysis (default: [freq[1], freq[end]])
- `weighting`: Boolean to indicate if the weighting based on the variance of each FRF is applied (default: true)
- `stabcrit`: Vector containing the stability criteria for natural frequencies and damping ratios (default: [0.01, 0.05])

**Output**
- `sol`: Data structure containing the results of the stabilization analysis
"""
function stabilization(frf, freq, max_order::Int, method = LSCF(); frange = [freq[1], freq[end]], weighting = true, stabcrit = [0.01, 0.05])

    # Initialization
    max_order += 1 # For having stability information from order 1 to max_order
    poles = [similar(freq, i) for i in 1:max_order]
    fn = [similar(freq, i) for i in 1:max_order]
    dr = [similar(freq, i) for i in 1:max_order]
    modefn = fill(NaN, max_order, max_order)
    mode_stabfn = fill(false, max_order, max_order)
    mode_stabdr = fill(false, max_order, max_order)

    for order in 1:max_order
        if method isa LSCE
            p = lsce(frf, freq, order, frange = frange, stabdiag = true)
        elseif method isa LSCF
            p = lscf(frf, freq, order, frange = frange, stabdiag = true, weighting = weighting)
        elseif method isa PLSCF
            p = plscf(frf, freq, order, frange = frange, stabdiag = true, weighting = weighting)
        else
            throw(ArgumentError("Unknown modal extraction method"))
        end

        fne, dre = poles2modal(p)
        poles[order] .= p
        fn[order] .= fne
        dr[order] .= dre
        modefn[1:order, order] .= fn[order]
        # Nmodes = length(fne)
        # poles[order] .= p[1:Nmodes]
        # fn[order] .= fne[1:Nmodes]
        # dr[order] .= dre[1:Nmodes]
        # modefn[1:Nmodes, order] .= fn[order]

        if order > 1
            Nm = length(fn[order-1])
            mode_stabfn[1: Nm, order-1], mode_stabdr[1:Nm, order-1] = check_stability(fn[order], fn[order-1], dr[order], dr[order-1], stabcrit)
        end
    end

    # Remove last order (not used for stability check)
    maxorder -= 1
    return EMAMdofStabilization(poles, modefn, mode_stabfn[:, 1:maxorder, end-1], mode_stabdr[:, 1:maxorder, end-1])
end