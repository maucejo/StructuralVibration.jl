## Structure definition for method selection
abstract type MdofModalExtraction end
struct LSCE <: MdofModalExtraction end
struct LSCF <: MdofModalExtraction end
struct PLSCF <: MdofModalExtraction end

"""
    EMAMdofProblem(frf, freq; frange, type_frf)

Data structure defining the inputs for EMA-MDOF modal extraction methods.

**Constructor parameters**
- `frf::Array{Complex, 3}`: 3D FRF matrix (array nm x ne x nf)
- `freq::AbstractArray{Real}`: Vector of frequency values (Hz)
- `frange::Vector{Real}`: Frequency range for analysis (default: [freq[1], freq[end]])
- `type_frf::Symbol`: Type of FRF used in the analysis
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

**Fields**
- `frf::Array{Complex, 3}`: 3D admittance matrix (array nm x ne x nf)
- `freq::AbstractArray{Real}`: Vector of frequency values (Hz)
- `type_frf::Symbol`: Type of FRF used in the analysis

**Note**
The FRF is internally converted to admittance if needed.
"""
@show_data struct EMAMdofProblem{C <: Complex, R <: Real}
    frf::Array{C, 3}
    freq::AbstractArray{R}
    type_frf::Symbol

    function EMAMdofProblem(frf::Array{C, 3}, freq::AbstractArray{R}; frange = [freq[1], freq[end]], type_frf = :dis) where {C <: Complex, R <: Real}

        # FRF post-processing - Frequency range reduction
        fidx = @. frange[1] ≤ freq ≤ frange[2]
        ω = 2π*freq[fidx]
        frf_red = frf[:, :, fidx]

        # Conversion to admittance
        if type_frf == :vel
            for (f, ωf) in enumerate(ω)
                frf_red[:, :, f] ./= (1im*ωf)
            end
        elseif type_frf == :acc
            for (f, ωf) in enumerate(ω)
                frf_red[:, :, f] ./= -ωf^2
            end
        end

        return new{C, R}(frf[:, :, fidx], freq[fidx], type_frf)
    end
end

"""
    AutoEMAMdofProblem(prob, dpi, method; modetype)

Structure containing the input data for automatic experimental modal analysis using Mdof methods

**Fields**
* `prob::EMAMdofProblem`: EMA-MDOF problem containing FRF data and frequency vector
* `order::Int`: Model order (number of poles to extract)
* `dpi::Vector{Int}`: Driving point indices - default = [1, 1]
    * `dpi[1]`: Driving point index on the measurement mesh
    * `dpi[2]`: Driving point index on the excitation mesh
* `method::MdofModalExtraction`: Method to extract the poles
    * `LSCE`: Least Squares Complex Exponential method
    * `LSCF``: Least Squares Complex Frequency method (default)
    * `PLSCF`: Polyreference Least Squares Complex Frequency method
* `modetype::Symbol`: Type of mode shapes to extract
    * `:real`: Real mode shapes (default)
    * `:complex`: Complex mode shapes
"""
@show_data struct AutoEMAMdofProblem
    prob::EMAMdofProblem
    order::Int
    dpi:: Vector{Int}
    method::MdofModalExtraction

    AutoEMAMdofProblem(prob::EMAMdofProblem, order::Int, dpi::Vector{Int} = [1, 1], method::MdofModalExtraction = LSCF()) = new(prob, order, dpi, method)
end

"""
    EMAMdofSolution(poles, ms, ci, res, lr, ur)

Structure containing the solution of the automatic experimental modal analysis using Mdof methods

**Fields**
* `poles::Vector{Complex}`: Extracted poles
* `ms::AbstractArray{Real}`: Mode shapes
* `ci::Vector{Complex}`: Scaling constants for each mode
* `res::AbstractArray{Complex}`: Residues for each mode
* `lr::Matrix{Complex}`: Lower residual
* `ur::Matrix{Complex}`: Upper residual
"""
@show_data struct EMAMdofSolution{Tp <: Complex, Tm <: Real}
    poles::Vector{Tp}
    ms::AbstractArray{Tm}
    ci::Vector{Tp}
    res::AbstractArray{Tp}
    lr::Matrix{Tp}
    ur::Matrix{Tp}
end

"""
    EMAMdofStabilization(prob, poles, modefn, mode_stabfn, mode_stabdr)

Data structure summarizing the results of the stabilization analysis.

**Fields**
- `prob::EMAMdofProblem`: EMA-MDOF problem containing FRF data and frequency vector
- `frange::Vector{Real}`: Frequency range used for the stabilization analysis
- `poles::Vector{Vector{Complex}}`: Vector of vectors containing extracted poles at each model order
- `modefn::Matrix{Real}`: Matrix containing the natural frequencies (useful for plotting)
- `mode_stabfn::Matrix{Bool}`: Matrix indicating the stability of natural frequencies
- `mode_stabdr::Matrix{Bool}`: Matrix indicating the stability of damping ratios

**Note**

This structure is returned by the `stabilization` function after performing a stabilization diagram analysis and used by `stabilization_plot` for visualization.
"""
@show_data struct EMAMdofStabilization{Tp <: Complex, Tf <: Real}
    prob::EMAMdofProblem        # EMA-MDOF problem containing FRF data and frequency vector
    poles::Vector{Vector{Tp}}   # Extracted poles at each model order
    modefn::Matrix{Tf}          # Natural frequencies (used for plotting)
    mode_stabfn::BitMatrix      # Stability of natural frequencies
    mode_stabdr::BitMatrix      # Stability of damping ratios
end

## Functions for modal extraction
"""
    poles_extraction(prob, order, method; stabdiag, weighting)

Extract complex poles from Frequency Response Function (FRF) data using the specified modal extraction method.

**Inputs**
- `prob::EMAMdofProblem`: EMA-MDOF problem containing FRF data and frequency vector
- `order::Int`: Model order (number of poles to extract)
- `method::MdofModalExtraction`: Modal extraction method
    - `LSCE`: Least Squares Complex Exponential method
    - `LSCF`: Least Squares Complex Frequency method (default)
    - `PLSCF`: Polyreference Least Squares Complex Frequency method
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)
- `weighting`: Boolean to indicate whether to apply weighting of the FRF (default: true)

**Output**
- `poles`: Vector of extracted complex poles

**Note**

The `weighting` parameter is only applicable for the LSCF and PLSCF methods.
"""
function poles_extraction(prob::EMAMdofProblem, order::Int, method::MdofModalExtraction = LSCF(); stabdiag = false, weighting = true)
    est_method = let
        if method isa LSCE
            lsce
        elseif method isa LSCF
            lscf
        elseif method isa PLSCF
            plscf
        else
            error("Unknown modal extraction method.")
        end
    end

    # Retry mechanism for modal extraction, since it can fail due to numerical issues
    estimation_success = false
    stop = 0
    local poles
    while !estimation_success && (stop ≤ 10)
        try
            # Modal extraction
             poles = est_method(prob, order, stabdiag = stabdiag, weighting = weighting)

            estimation_success = true
        catch e # Estimation failed
            @warn "Modal extraction failed at order $order with error: $e. Retrying..."
            stop += 1
        end
    end

    if stop == 10
        error("Modal extraction failed after 10 attempts at order $order.")
    end

    return poles
end


"""
    lsce(prob, order; stabdiag, weighting)

Perform Least Squares Complex Exponential (LSCE) method to extract complex poles from Frequency Response Function (FRF) data.

**Inputs**
- `prob`: EMAMdofProblem containing FRF data and frequency vector
- `order::Int`: Model order (number of poles to extract)
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)
- `weighting`: Boolean to indicate whether to apply weighting of the FRF (default: true)

**Output**
- `poles`: Vector of extracted complex poles

**Note**
The weighting parameter is not used in the LSCE method but is included for consistency with other methods.
"""
function lsce(prob::EMAMdofProblem, order::Int; stabdiag = false, weighting = true)
    # Extract FRF and frequency from problem
    (; frf, freq) = prob

    # Approximate sampling frequency of the truncated signal
    fsred = 2.56*(freq[end]-freq[1])

    # Impulse function calculation
    nm, ne, nf = size(frf)
    FRF = reshape(frf, (nm*ne, nf))
    Hk = impulse_response(FRF, freq, fsred)

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

"""
    lscf(prob, order; stabdiag, weighting)

Perform Least Squares Complex Frequency (LSCF) method to extract complex poles from Frequency Response Function (FRF) data.

**Inputs**
- `prob`: EMAMdofProblem containing FRF data and frequency vector
- `order::Int`: Model order (number of poles to extract)
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)
- `weighting`: Boolean to indicate if the weighting based on the variance of each FRF is applied (default: true)

**Output**
- `poles`: Vector of extracted complex poles

**Reference**

[1] El-Kafafy M., Guillaume P., Peeters B., Marra F., Coppotelli G. (2012).Advanced Frequency-Domain Modal Analysis for Dealing with Measurement Noise and Parameter Uncertainty. In: Allemang R., De Clerck J., Niezrecki C., Blough J. (eds) Topics in Modal Analysis I, Volume 5. Conference Proceedings of the Society for Experimental Mechanics Series. Springer, New York, NY
"""
function lscf(prob::EMAMdofProblem, order::Int; stabdiag = false, weighting = true) :: Vector{eltype(prob.frf)}
    # Extract FRF and frequency from problem
    (; frf, freq) = prob

    # FRF post-processing - Reshape FRF matrix
    (nm, ne, nf) = size(frf)
    FRF = reshape(frf, (nm*ne, nf))  # each row = one FRF across frequencies

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
        M .+= (Tk .- Sk'*(Rk\Sk))/vk
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

    # Poles filtering (keep negative real part and negative imag to pick one of conjugates)
    valid_poles = p[@. imag(p) < 0. && real(p) ≤ 0. && !isreal(p)]
    poles = intersect(p, conj.(valid_poles))
    sort!(poles, by = abs)

    # Frequency offset correction (restore original frequency offset)
    if freq[1] != 0. && !isempty(poles)
        fn, ξn = poles2modal(poles)
        poles .= modal2poles(fn .+ freq[1], @. ξn * fn / (fn + freq[1]))
    end

    # Return format based on stabdiag flag
    return stabdiag ? [poles; fill(complex(NaN, NaN), order - length(poles))] : poles
end

"""
    plscf(prob order; stabdiag, weighting)

Perform Polyreference Least Squares Complex Frequency (pLSCF) method to extract complex poles from Frequency Response Function (FRF) data.

**Inputs**
- `prob`: EMAMdofProblem containing FRF data and frequency vector
- `order::Int`: Model order (number of poles to extract)
- `frange`: Frequency range for analysis (default: [freq[1], freq[end]])
- `stabdiag`: Boolean to indicate the function is used to build a stability diagram (default: false)
- `weighting`: Boolean to indicate if the weighting based on the variance of each FRF is applied (default: true)

**Output**
- `poles`: Vector of extracted complex poles

**Reference**

[1] El-Kafafy M., Guillaume P., Peeters B., Marra F., Coppotelli G. (2012).Advanced Frequency-Domain Modal Analysis for Dealing with Measurement Noise and Parameter Uncertainty. In: Allemang R., De Clerck J., Niezrecki C., Blough J. (eds) Topics in Modal Analysis I, Volume 5. Conference Proceedings of the Society for Experimental Mechanics Series. Springer, New York, NY
"""
function plscf(prob::EMAMdofProblem, order::Int; stabdiag = false, weighting = true) :: Vector{eltype(prob.frf)}
    # Extract FRF and frequency from problem
    (; frf, freq) = prob

    # FRF post-processing - Reshape FRF matrix
    (nm, ne, nf) = size(frf)
    if ne > nm
        frf = permutedims(frf, (2, 1, 3))
        (nm, ne, nf) = size(frf)
    end

    # pLSCF computation
    ω = 2(freq .- freq[1]) # Reduced angular frequency
    Δt = 1/(2(freq[end] - freq[1])) # Sampling interval
    modelOrder = 0:order
    nmodel = order + 1

    # Reduced normal equation computation
    Ω0 = cispi.(-ω*modelOrder'*Δt)
    # Precompute R0
    X0 = similar(Ω0, nf, ne*nmodel)
    for f in 1:nf
        X0[f, :] .= kron(Ω0[f, :], I(ne))
    end
    R0 = real(X0'X0)

    # Preallocation
    M = similar(freq, ne*nmodel, ne*nmodel)
    Y0 = similar(Ω0, nf, ne*nmodel)
    H0 = similar(Ω0, ne, nf)
    S0 = similar(M)
    T0 = similar(M)
    vk = ones(ne)

    for i in 1:nm
        H0 .= frf[i, :, :]
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
        M .+= (T0 .- S0'*(R0\S0))./vk
    end

    # Computation of the coefficients of the denominator
    A = -M[1:order*ne, 1:order*ne]
    b = M[1:order*ne, (order*ne + 1):nmodel*ne]
    α = A\b

    # Construct the companion matrix
    Id = I(ne*(order - 1))
    Ze = zeros(ne*(order - 1), ne)
    C = vcat(hcat(Ze, Id), -α')

    # Pole calculation
    p = similar(frf, ne*order)
    try
        p .= -log.(eigvals(C))/Δt
    catch e
        return fill(complex(NaN, NaN), order)
    end

    # Poles that do not appear as pairs of complex conjugate numbers with a positive real part or that are purely real are suppressed. For complex conjugate poles, only the pole with a positive imaginary part is retained.
    valid_poles = p[@. imag(p) < 0. && real(p) ≤ 0. && !isreal(p)]
    poles = intersect(p, conj.(valid_poles))

    # To avoid modifying the stability of the dynamic system in case of frequency offset, we set real(poles_new) = real(poles_old). This means that dn_new = dn_old * fn_old / fn_new.
    if freq[1] != 0.
        fn, ξn = poles2modal(poles)
        poles .= modal2poles(fn .+ freq[1], @. ξn*fn/(fn + freq[1]))
    end
    sort!(poles, by = abs)

    # If Stability diagram
    return stabdiag ? [poles; fill(complex(NaN, NaN), order - length(poles))] : poles
end

## Function for stabilization diagram analysis
"""
    stabilization(prob, max_order, method; weighting, stabcrit)

Perform stabilization diagram analysis using the specified modal extraction method (LSCE, LSCF, or PLSCF).

**Inputs**
- `prob::EMAMdofProblem`: EMA-MDOF problem containing FRF data and frequency vector
- `max_order::Int`: Maximum model order for the stabilization analysis
- `method`: Modal extraction method to use
    - `LSCE()`: Least Squares Complex Exponential method
    - `LSCF()`: Least Squares Complex Frequency method (default)
    - `PLSCF()`: Polyreference Least Squares Complex Frequency method
- `frange`: Frequency range for analysis (default: [freq[1], freq[end]])
- `weighting`: Boolean to indicate if the weighting based on the variance of each FRF is applied (default: true)
- `stabcrit`: Vector containing the stability criteria for natural frequencies and damping ratios (default: [0.01, 0.05])

**Output**
- `sol::EMAMdofStabilization`: Data structure containing the results of the stabilization analysis
"""
function stabilization(prob::EMAMdofProblem, max_order::Int, method::MdofModalExtraction = LSCF(); weighting = true, stabcrit = [0.01, 0.05])

    # Initialization
    max_order += 1 # For having stability information from order 1 to max_order
    poles = [similar(prob.frf, i) for i in 1:max_order]
    fn = [similar(prob.freq, i) for i in 1:max_order]
    dr = [similar(prob.freq, i) for i in 1:max_order]
    modefn = fill(NaN, max_order, max_order)
    mode_stabfn = falses(max_order, max_order)
    mode_stabdr = falses(max_order, max_order)

    for order in 1:max_order
        poles[order] .= poles_extraction(prob, order, method, stabdiag = true, weighting = weighting)

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
    return EMAMdofStabilization(prob, poles[1:max_order], modefn[1:max_order, 1:max_order], mode_stabfn[1:max_order, 1:max_order], mode_stabdr[1:max_order, 1:max_order])
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
    (; frf, freq) = prob

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
    modeshape_extraction(residues, poles, dpi; type = :complex)

Extract mode shapes using MDOF approximation

**Inputs**
- `residues`: Residues matrix of size (np, nm, ne)
- `poles`: Vector of complex poles
* `dpi`: Driving point indices - default = [1, 1]
    * `dpi[1]`: Driving point index on the measurement mesh
    * `dpi[2]`: Driving point index on the excitation mesh
- `type`: Type of mode shape
    * `:complex`: Complex mode shapes (default)
    * `:real`: Real mode shapes

**Outputs**
- `ms`: Mode shapes matrix
- `ci`: Scaling factors vector
"""
function modeshape_extraction(residues, poles::Vector{T}, dpi = [1, 1]; type = :complex) where {T <: Complex}

    # Data preparation
    np, nm, ne = size(residues)
    if nm < ne
        # Roving hammer test
        Res = permutedims(residues, (3, 1, 2))
        dpi = [dpi[2], dpi[1]]
    else
        # Roving accelerometer test
        Res = permutedims(residues, (2, 1, 3))
    end

    # Scaling factor calculation
    sqR = sqrt.(Res[dpi[1], :, dpi[2]])
    sqR[abs.(sqR) .≤ eps()] .= one(poles[1]) # Avoid division by zero

    # Mode shape extraction
    ms = Res[:, :, dpi[2]] ./ transpose(sqR)

    # Scaling mode shapes to unit modal mass if real mode shapes are requested
    ci = ones(eltype(ms), np)
    if type == :real
        fn, ξn = poles2modal(poles)
        ωn = 2π*fn
        @. ci /= 2im*ωn*√(1. - ξn^2)
        ms ./= sqrt.(transpose(ci))
    end

    return ms, ci
end

"""
    solve(prob::AutoEMAMdofProblem)

Solve automatically experimental modal analysis problem using Mdof methods

**Inputs**
* `prob`: Structure containing the input data for automatic experimental modal analysis using Mdof methods

**Outputs**
* `sol::EMAMdofSolution`: Structure containing the solution of the automatic experimental modal analysis using Mdof methods
"""
function solve(prob_ema::AutoEMAMdofProblem)
    # Unpack problem data
    (; prob, order, dpi, method) = prob_ema

    # Extraction of natural frequencies and damping ratios
    poles = poles_extraction(prob, order, method)

    # Extraction of mode shapes
    res = mode_residues(prob, poles)

    # Compute the residuals separately for accuracy
    lr, ur = compute_residuals(prob, res, poles)

    phi, ci = modeshape_extraction(res, poles, dpi, type = :real)

    return EMAMdofSolution(poles, c2r_modeshape(phi), ci, res, lr, ur)
end