"""
    poles_extraction(prob::EMAProblem, alg)
    poles_extraction(prob::MdofProblem, order, alg; stabdiag)

Extract poles from the Bode diagram fitting method

**Inputs**
* `prob::EMAProblem`: EMA problem containing FRF data and frequency vector
* `order`: Order of the model (only for Mdof methods)
* `alg`: Algorithm to extract the poles
    * Sdof methods:
        * `PeakPicking`: Peak picking method (default for Sdof methods)
        * `CircleFit`: Circle fitting method
        * `LSFit`: Least squares fitting method
    * Mdof methods:
        * `LSCF`: Least Squares Complex Frequency method (default for Mdof methods)
        * `PLSCF`: Polyreference Least Squares Complex Frequency method
        * `LSCE`: Least Squares Complex Exponential method (only for Mdof methods)
* `stabdiag::Bool`: Boolean to indicate the function is used to build a stability diagram (Only for Mdof methods, default: false)

**Outputs**
* `poles`: Vector of extracted complex poles

**Note**
- For Sdof methods, the natural frequencies and damping ratios are extracted from each FRF (each row of the matrix) and then averaged. The number of FRF used for averaging are those having the maximum (and same) number of peaks detected.
"""
function poles_extraction(prob::MdofProblem, order::Int, alg::MdofModalExtraction = LSCF(); stabdiag = false)

    return compute_poles(prob, order, alg, stabdiag = stabdiag)
end

"""
    compute_poles(prob, order, alg = LSCE(); stabdiag, weighting)

Perform Least Squares Complex Exponential (LSCE) method to extract complex poles from Frequency Response Function (FRF) data.

**Inputs**
- `prob`: EMAMdofProblem containing FRF data and frequency vector
- `order::Int`: Model order (number of poles to extract)
- `alg::LSCE`: Modal extraction method
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)

**Output**
- `poles`: Vector of extracted complex poles
"""
@views function compute_poles(prob::MdofProblem, order::Int, alg::LSCE; stabdiag = false)
    # Extract FRF and frequency from problem
    if prob isa EMAProblem
        (; frf, freq) = prob
    elseif prob isa OMAProblem
        frf = prob.halfspec
        freq = prob.freq
    else
        error("Unsupported problem type for LSCE method.")
    end

    # Approximate sampling frequency of the truncated signal
    fsred = 2.56*(freq[end]-freq[1])

    # Impulse function calculation
    no, ni, nf = size(frf)
    FRF = reshape(frf, (no*ni, nf))
    Hk = impulse_response(FRF, freq, fsred)

    # Model order and number of samples
    n = 2order
    nsamples = 10n # Number of samples used for the identification

    # Preallocation
    H0 = similar(Hk, nsamples, n+1)
    A = similar(Hk, nsamples*no*ni, n)
    b = similar(Hk, nsamples*no*ni)
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

    return poles_validity(p, order, freq, stabdiag)
end

"""
    compute_poles(prob, order, alg = LSCF(); stabdiag)

Perform Least Squares Complex Frequency (LSCF) method to extract complex poles from Frequency Response Function (FRF) data.

**Inputs**
- `prob`: EMAMdofProblem containing FRF data and frequency vector
- `order::Int`: Model order (number of poles to extract)
- `alg::LSCF`: Modal extraction method
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)

**Output**
- `poles`: Vector of extracted complex poles

**Reference**

[1] El-Kafafy M., Guillaume P., Peeters B., Marra F., Coppotelli G. (2012).Advanced Frequency-Domain Modal Analysis for Dealing with Measurement Noise and Parameter Uncertainty. In: Allemang R., De Clerck J., Niezrecki C., Blough J. (eds) Topics in Modal Analysis I, Volume 5. Conference Proceedings of the Society for Experimental Mechanics Series. Springer, New York, NY
"""
@views function compute_poles(prob::MdofProblem, order::Int, alg::LSCF; stabdiag = false)
    # Extract FRF and frequency from problem
    if prob isa EMAProblem
        (; frf, freq) = prob
    elseif prob isa OMAProblem
        frf = prob.halfspec
        freq = prob.freq
    else
        error("Unsupported problem type for LSCF method.")
    end

    # FRF post-processing - Reshape FRF matrix
    (no, ni, nf) = size(frf)
    FRF = reshape(frf, (no*ni, nf))  # each row = one FRF across frequencies

    # LSCF computation - stable numerics
    ω = 2(freq .- freq[1])                # Reduced angular frequency
    Δt = 1/(2(freq[end] - freq[1]))    # Sampling interval
    modelOrder = 0:order
    nmodel = order + 1

    # Build basis matrix X0 (nf x nmodel), complex
    # Use the same basis as original: exp.(-1im*ω*modelOrder'*Δt)
    X0 = cispi.(-ω*modelOrder'*Δt)  # nf x nmodel
    Rk = real(X0'X0)

    # Preallocation
    M = zeros(eltype(freq), nmodel, nmodel)
    Yk = similar(X0)
    Sk = similar(M)
    Tk = similar(M)

    # Process each FRF (each row of FRF is length nf)
    for Hk in eachrow(FRF)
        # Build Yk, Sk, Tk
        @. Yk = -Hk*X0
        Sk .= real(X0'Yk)
        Tk .= real(Yk'Yk)

        # Accumulate robust contribution to M
        M .+= (Tk .- Sk'*(Rk\Sk))
    end

    # Extract submatrices for denominator solve
    A = -M[1:order, 1:order]
    b = M[1:order, nmodel]
    α = [A\b; 1.]

    # Compute poles from polynomial roots
    V = roots(Polynomial(α))
    p = similar(frf, order)
    try
        p .= -log.(V)/Δt
    catch e
        return fill(complex(NaN, NaN), order)
    end

    return poles_validity(p, order, freq, stabdiag)
end

"""
    compute_poles(prob, order, alg = PLSCF(); stabdiag)

Perform Polyreference Least Squares Complex Frequency (pLSCF) method to extract complex poles from Frequency Response Function (FRF) data.

**Inputs**
- `prob`: EMAMdofProblem containing FRF data and frequency vector
- `order::Int`: Model order (number of poles to extract)
- `alg::PLSCF`: Modal extraction method
- `frange`: Frequency range for analysis (default: [freq[1], freq[end]])
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)

**Output**
- `poles`: Vector of extracted complex poles

**Reference**

[1] El-Kafafy M., Guillaume P., Peeters B., Marra F., Coppotelli G. (2012).Advanced Frequency-Domain Modal Analysis for Dealing with Measurement Noise and Parameter Uncertainty. In: Allemang R., De Clerck J., Niezrecki C., Blough J. (eds) Topics in Modal Analysis I, Volume 5. Conference Proceedings of the Society for Experimental Mechanics Series. Springer, New York, NY
"""
@views function compute_poles(prob::MdofProblem, order::Int, alg::PLSCF; stabdiag = false)
    # Extract FRF and frequency from problem
    if prob isa EMAProblem
        (; frf, freq) = prob
    elseif prob isa OMAProblem
        frf = prob.halfspec
        freq = prob.freq
    else
        error("Unsupported problem type for PLSCF method.")
    end

    # FRF post-processing - Reshape FRF matrix
    (no, ni, nf) = size(frf)
    if ni > no
        frf = permutedims(frf, (2, 1, 3))
        (no, ni, nf) = size(frf)
    end

    # pLSCF computation
    ω = 2(freq .- freq[1]) # Reduced angular frequency
    Δt = 1/(2(freq[end] - freq[1])) # Sampling interval
    modelOrder = 0:order
    nmodel = order + 1

    # Reduced normal equation computation
    X0 = cispi.(-ω*modelOrder'*Δt)
    R0 = real(X0'X0)

    # Preallocation
    M = zeros(eltype(freq), ni*nmodel, ni*nmodel)
    Y0 = similar(X0, nf, ni*nmodel)
    H0 = similar(X0, ni, nf)
    S0 = similar(X0, nmodel, ni*nmodel)
    T0 = similar(M)

    for i in 1:no
        H0 .= frf[i, :, :]

        for f in 1:nf
            Y0[f, :] .= -kron(X0[f, :], H0[:, f])
        end

        # Build R0, S0, T0
        S0 .= real(X0'Y0)
        T0 .= real(Y0'Y0)

        # Accumulation of M
        M .+= (T0 .- S0'*(R0\S0))
    end

    # Computation of the coefficients of the denominator
    A = -M[1:order*ni, 1:order*ni]
    b = M[1:order*ni, (order*ni + 1):nmodel*ni]

    # Check condition number for better numerical stability
    if cond(A) > 1e12
        # Use SVD for better stability
        F = svd(A)
        # Filter out very small singular values
        tol = maximum(F.S) * eps(eltype(F.S)) * max(size(A)...)
        inv_S = [s > tol ? 1/s : 0.0 for s in F.S]
        α = F.V * Diagonal(inv_S) * F.U' * b
    else
        α = A\b
    end

    # Construct the companion matrix
    Id = I(ni*(order - 1))
    Ze = zeros(ni*(order - 1), ni)
    C = vcat(hcat(Ze, Id), -α')

    # Pole calculation
    p = similar(frf, ni*order)
    try
        p .= -log.(eigvals(C))/Δt
    catch e
        return fill(complex(NaN, NaN), order)
    end

    return poles_validity(p, order, freq, stabdiag)
end

## Function for stabilization diagram analysis
"""
    stabilization(prob, max_order, method; stabcrit)

Perform stabilization diagram analysis using the specified modal extraction method.

**Inputs**
- `prob::MdofProblem`: EMA-MDOF problem containing FRF data and frequency vector
- `max_order::Int`: Maximum model order for the stabilization analysis
- `alg::Union{MdofModalExtraction, OMAModalExtraction}`: Modal extraction algorithm
    - EMA algorithms:
        - `LSCE()`: Least Squares Complex Exponential method
        - `LSCF()`: Least Squares Complex Frequency method (default)
        - `PLSCF()`: Polyreference Least Squares Complex Frequency method
    - OMA algorithms:
        - `CovSSI()`: Covariance-based SSI method
        - `DataSSI()`: Data-based SSI method
- `frange`: Frequency range for analysis (default: [freq[1], freq[end]])
- `stabcrit`: Vector containing the stability criteria for natural frequencies and damping ratios (default: [0.01, 0.05])

**Output**
- `sol::EMAMdofStabilization`: Data structure containing the results of the stabilization analysis
"""
function stabilization(prob::MdofProblem, max_order::Int, alg::Union{MdofModalExtraction, OMAModalExtraction} = LSCF(); stabcrit = [0.01, 0.05])

    # Extract FRF and frequency from problem
    if prob isa EMAProblem
        frf = prob.frf
    elseif prob isa OMAProblem
        frf = prob.halfspec
    else
        error("Unsupported problem type for LSCF method.")
    end

    # Initialization
    max_order += 1 # For having stability information from order 1 to max_order
    poles = [similar(frf, i) for i in 1:max_order]
    fn = [similar(prob.freq, i) for i in 1:max_order]
    dr = [similar(prob.freq, i) for i in 1:max_order]
    modefn = fill(NaN, max_order, max_order)
    mode_stabfn = falses(max_order, max_order)
    mode_stabdr = falses(max_order, max_order)

    for order in 1:max_order
        try
            poles[order] .= poles_extraction(prob, order, alg, stabdiag = true)
        catch e
            poles[order] .= fill(complex(NaN, NaN), order)
        end

        fne, dre = poles2modal(poles[order])
        fn[order] .= fne
        dr[order] .= dre
        modefn[1:order, order] .= fn[order]

        if order > 1
            Nm = length(fn[order-1])
            mode_stabfn[1:Nm, order-1], mode_stabdr[1:Nm, order-1] = check_stability(fn[order], fn[order-1], dr[order], dr[order-1], stabcrit)
        end
    end

    # Remove last order (not used for stability check)
    max_order -= 1
    return StabilizationAnalysis(prob, poles[1:max_order], modefn[1:max_order, 1:max_order], mode_stabfn[1:max_order, 1:max_order], mode_stabdr[1:max_order, 1:max_order])
end

"""
    check_stability(fn_new, fn_old, dr_new, dr_old, stabcrit)

Check the stability of natural frequencies and damping ratios between two consecutive model orders.

**Inputs**
- `fn_new`: Vector of natural frequencies at the new model order
- `fn_old`: Vector of natural frequencies at the old model order
- `dr_new`: Vector of damping ratios at the new model order
- `dr_old`: Vector of damping ratios at the old model order
- `stabcrit`: Vector containing the stability criteria for natural frequencies and damping ratios

**Outputs**
- `stabfn`: Boolean vector indicating the stability of natural frequencies
- `stabdr`: Boolean vector indicating the stability of damping ratios
"""
function check_stability(fn_new, fn_old, dr_new, dr_old, stabcrit)
    Nmodes_old = length(fn_old)
    stabfn = falses(Nmodes_old)
    stabdr = falses(Nmodes_old)

    for i in 1:Nmodes_old
        if isnan(fn_old[i])
            continue
        end

        # Find closest mode in new order
        idx = argmin(skipnan(abs.(fn_new .- fn_old[i])))

        # Compute differences
        dfn = abs(fn_new[idx] - fn_old[i])
        drn = abs(dr_new[idx] - dr_old[i])

        # Check frequency stability
        stabfn[i] = dfn ≤ stabcrit[1]*fn_old[i]

        # Check damping ratio stability
        stabdr[i] = drn ≤ stabcrit[2]*dr_old[i]
    end

    return stabfn, stabdr
end

## Mode extraction functions
"""
    mode_residues(prob, poles)

Compute residues of a frequency response function (FRF) given its poles.

**Inputs**
- `prob`: EMAMdofProblem containing FRF data and frequency vector
- `poles`: Vector of complex poles

**Output**
- `res`: Residues corresponding to each pole
"""
function mode_residues(prob, poles)
    # Extract FRF and frequency from problem
    if prob isa EMAProblem
        (; frf, freq) = prob
    elseif prob isa OMAProblem
        frf = prob.halfspec
        freq = prob.freq
    else
        error("Unsupported problem type for mode residues calculation.")
    end

    # Correct frange to avoid division by zero
    if freq[1] < 1.
        freq_red = freq[2:end]
        frf_red = permutedims(frf[:, :, 2:end], (3, 1, 2)) # Put the last dimension first for subsequent calculations
    else
        freq_red = freq
        frf_red = permutedims(frf, (3, 1, 2)) # Put the last dimension first for subsequent calculations
    end
    nm, ne = size(frf_red)[2:3]
    ω = 2π*freq_red

    # Residue calculation
    valid_poles = .!isnan.(poles)
    p = [poles[valid_poles]; conj.(poles[valid_poles])]
    np = length(p)

    # The residues are computing separately for the sake of accuracy
    res = similar(frf_red, np, nm, ne)
    P = (1 ./(1im*ω .- transpose(p)))
    for i in 1:ne
        res[:, :, i] .= P\frf_red[:, :, i]
    end

    # Return residues
    return res[1:(np ÷ 2), : , :]
end

"""
    modeshape_extraction(prob, poles, alg; dpi)
    modeshape_extraction(residues, poles, alg; dpi, modetype)

Extract mode shapes using Sdof approximation

**Inputs**
* `prob::EMAProblem`: EMA problem containing FRF data and frequency vector
* `poles`: Vector of complex poles
* `alg::SdofModalExtraction` or `alg::MdofModalExtraction`: Modal extraction algorithm
* `dpi`: Driving point indices - default = [1, 1]
    * `dpi[1]`: Driving point index on the measurement mesh
    * `dpi[2]`: Driving point index on the excitation mesh
* `modetype`: Type of mode shape (only for Mdof methods)
    - `:emac`: EMA - Complex mode shapes (default)
    - `:emar`: EMA - Real mode shapes
    - `:oma`: OMA

**Output**
* `ms`: Mode shapes
* `ci`: Scaling factors vector (only for Mdof methods)

**Note**
- If the number of measurement points is less than the number of excitation points, the mode shapes are estimated at the excitation points (roving hammer test). Otherwise, the mode shapes are estimated at the measurement points (roving accelerometer test).
- The `alg` argument is here for performing multiple dispatch but is not used in the function.
- If `modetype` is set to `:oma`, the mode shapes correspond to the first left singular vectors of the residues matrices (see Ref. [2]).

**References**
[1] M. Géradin and D. J. Rixen. "Mechanical Vibrations: Theory and Application to Structural Dynamics". 3rd. Edition, Wiley, 2015.

[2] C. Rainieri and G. Fabbrocino. "Operational Modal Analysis of Civil Engineering Structures: An Introduction and Guide for Applications". Springer, 2014.

[3] P. Verboven. "Frequency-domain system identification for modal analysis". PhD thesis, Katholieke Universiteit Leuven, 2002.
"""
function modeshape_extraction(residues, poles::Vector{T}, alg::MdofModalExtraction; dpi = [1, 1], modetype = :emac) where {T <: Complex}

    # Data preparation
    np, nm, ne = size(residues)
    if nm < ne
        # Roving hammer test
        Res = permutedims(residues, (3, 1, 2))
        dpi = [dpi[2], dpi[1]]
        nm = ne
    else
        # Roving accelerometer test
        Res = permutedims(residues, (2, 1, 3))
    end

    # Modal constant initialization
    ci = ones(eltype(Res), np)
    ms = similar(Res, nm, np)
    if modetype == :emar || modetype == :emac
        # Scaling factor calculation
        sqR = sqrt.(Res[dpi[1], :, dpi[2]])
        sqR[abs.(sqR) .≤ eps()] .= one(poles[1]) # Avoid division by zero

        # Mode shape extraction
        ms .= Res[:, :, dpi[2]] ./ transpose(sqR)

        # Scaling mode shapes to unit modal mass if real mode shapes are requested
        if modetype == :emar
            fn, ξn = poles2modal(poles)
            ωn = 2π*fn
            @. ci /= 2im*ωn*√(1. - ξn^2)
            ms ./= sqrt.(transpose(ci))
        end
    else # OMA - See Ref. [2, 3]
        for i in 1:np
            # Compute SVD of the residue matrix
            F = svd(Res[:, i, :])

            # Filter significant singular values (relative threshold)
            S_norm = F.S / maximum(F.S)
            significant_idx = findall(S_norm .≥ 0.01)  # Adjustalble threshold

            # Initialize best candidate
            best_idx = 1

            # If first mode, use energy criterion
            # Otherwise, use orthogonality with previous modes
            if i == 1
                # For first mode, select based on energy distribution
                # (mode shapes should not be too localized)
                for idx in significant_idx
                    candidate = F.U[:, idx]
                    energy_spread = std(abs.(candidate)) / mean(abs.(candidate))

                    if energy_spread > 0.3  # Adjustable threshold
                        best_idx = idx
                        break
                    end
                end
            else
                # For subsequent modes, check orthogonality with previous modes
                for idx in significant_idx
                    candidate = F.U[:, idx]

                    # Compute MAC with all previous modes
                    mac_prev = [mac(candidate, ms[:, j]) for j in 1:(i-1)]

                    # Select mode with the lowest MAC to previous modes
                    if maximum(mac_prev) < 0.3  # Orthogonality
                        best_idx = idx
                        break
                    end
                end
            end

            ms[:, i] .= F.U[:, best_idx]
        end
    end

    return ms, ci
end

"""
    solve(prob_ema::AutoEMASdofProblem)
    solve(prob_ema::AutoEMAMdofProblem, order::Int; stabdiag, weighting)

Solve automatically experimental modal analysis problem using Sdof or Mdof methods

**Inputs**
* `prob_ema`: Structure containing the input data for automatic experimental modal analysis using Sdof methods

**Outputs**
* `sol::EMASolution`: Structure containing the solution of the automatic experimental modal analysis using Sdof methods
"""
function solve(prob_ema::AutoEMAMdofProblem)
    # Unpack problem data
    (; prob, order, dpi, alg) = prob_ema

    # Extraction of natural frequencies and damping ratios
    poles = poles_extraction(prob, order, alg)

    # Extraction of mode shapes
    res = mode_residues(prob, poles)

    # Compute the residuals separately for accuracy
    lr, ur = compute_residuals(prob, res, poles)

    phi, ci = modeshape_extraction(res, poles, alg, dpi = dpi, modetype = :emar)

    return EMASolution(poles, c2r_modeshape(phi), ci, res, lr, ur)
end