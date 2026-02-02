"""
    poles_extraction(prob::EMAProblem, alg; width, min_prom, max_prom, pks_indices, scaling)
    poles_extraction(prob::MdofProblem, order, alg; stabdiag)

Extract poles using Sdof or Mdof experimental modal analysis methods

**Inputs**
* `prob::EMAProblem` or `prob::MdofProblem`: EMA problem containing FRF data and frequency vector
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
* `width::Int`: Half-width of the peaks (only for Sdof methods, default: 1)
* `min_prom::Real`: Minimum peak prominence (only for Sdof methods, default: 0. dB)
* `max_prom::Real`: Maximum peak prominence (only for Sdof methods, default: Inf)
* `pks_indices::Vector{Int}`: Indices of peaks to consider (only for Sdof methods, default: empty vector)

**Outputs**
* `poles`: Vector of extracted complex poles

**Note**
- For Sdof methods, the natural frequencies and damping ratios are extracted from each FRF (each row of the matrix) and then averaged. The number of FRF used for averaging are those having the maximum (and same) number of peaks detected.
"""
function poles_extraction(prob::EMAProblem, alg::SdofEMA; width::Int = 1, min_prom = 0., max_prom = Inf, pks_indices = Int[])
    # Extract FRF and frequency from problem
    (; frf, freq) = prob

    no, ni, nf = size(frf)
    Hr = reshape(frf, no*ni, nf)

    # Estimate the number of peaks in the FRF
    Hm = vec(20log10.(mean(abs, Hr, dims = 1)))
    if !isempty(pks_indices)
        pks = (indices = pks_indices, heights = Hm[pks_indices], data = Hm)
    else
        pks = findmaxima(Hm, width)
    end

    pks = peakproms!(pks, min = min_prom, max = max_prom) |> peakwidths!
    npeak = length(pks.indices)

    if npeak == 0
        throw(ErrorException("No peaks found in the FRF data. Adjust the peak detection parameters."))
    end

    # Initialization
    pk = similar(frf, no*ni, npeak)
    pn = similar(frf, npeak)
    # keeprow = trues(no*ni)
    poles = similar(frf, Complex{eltype(freq)}, npeak)
    for (k, Hv) in enumerate(eachrow(Hr))
        p = compute_poles(Hv, freq, alg, pks)

        nk = length(p)
        poles .= [p; fill(complex(NaN, NaN), npeak - nk)]
        # if any(isnan.(poles))
        #     keeprow[k] = false
        # end

        pk[k, :] .= poles
    end

    # Average the results from different FRFs
    for i in 1:npeak
        # pn[i] = mean(skipnan(pk[keeprow, i]))
        pn[i] = mean(skipnan(pk[:, i]))
    end

    return pn
end

"""
    compute_poles(H, freq, alg, width, min_prom, max_prom)


Extract poles from the peak picking method

**Inputs**
* `H::Array{Complex, 3}`: Frequency response function
* `freq::AbstractArray`: Frequency vector
* `width::Int`: Half-width of the peaks
* `alg`: Algorithm to extract the poles
    * `PeakPicking`: Peak picking method (default for Sdof methods)
    * `CircleFit`: Circle fitting method
    * `LSFit`: Least squares fitting method
* `pks`: NamedTuple containing indices, heights, prominences, widths, edges

**Outputs**
* `poles`: Extracted poles

**References**
[1] D.J. Ewins, "Modal Testing: Theory, Practice and Application", 2nd Edition, Research Studies Press, 2000.

[2] A. Brandt, "Noise and Vibration Analysis: Signal Analysis and Experimental Procedures", Wiley, 2011.
"""
function compute_poles(H, freq, alg::PeakPicking, pks)

    # Flat FRF - No peaks
    Habs = abs.(H)

    # Check if maximum and minimum are the same or if the coefficient of variation is too low
    if maximum(Habs) == minimum(Habs) || std(Habs[Habs .> 0]) < eltype(freq)(0.1)*mean(Habs[Habs .> 0])
        return Complex{eltype(freq)}[]
    end

    # Estimation of natural frequencies and damping ratios
    ftemp = freq[pks.indices]
    fn = similar(ftemp)
    ξn = similar(ftemp)

    # Damping ratios estimation from the half-bandwidth method
    nfreq_itp = 250
    freq_left = similar(fn, nfreq_itp)
    freq_right = similar(freq_left)
    Hleft = similar(freq_left)
    Hright = similar(freq_left)
    for (n, (f, idmax, edg)) in enumerate(zip(ftemp, pks.indices, pks.edges))
        Hmax = Habs[idmax]
        edge1 = floor(Int, edg[1])
        edge2 = ceil(Int, edg[2])

        # Left side of the peak
        itp_left = linear_interpolation(freq[edge1:idmax], Habs[edge1:idmax])
        freq_left .= LinRange(freq[edge1], f, nfreq_itp)
        @. Hleft = itp_left(freq_left)
        posleft = argmin(abs.(Hleft .- Hmax/√2))
        fmin = freq_left[posleft]

        # Right side of the peak
        itp_right = linear_interpolation(freq[idmax:edge2], Habs[idmax:edge2])
        freq_right .= LinRange(f, freq[edge2], nfreq_itp)
        @. Hright = itp_right(freq_right)
        posright = argmin(abs.(Hright .- Hmax/√2))
        fmax = freq_right[posright]

        # A structural damping is supposed η = 2ξ. So ξ = η/2
        ξn[n] = (fmax^2 - fmin^2)/4f^2
        if ξn[n] > 1. || ξn[n] < 0.
            ξn[n] = NaN
        end

        # Correction of the natural frequency due to damping
        fn[n] = ftemp[n]/√(1 - ξn[n]^2)
    end

    return modal2poles(fn, ξn)
end

function compute_poles(H, freq, alg::CircleFit, pks)

    # Flat FRF - No peaks
    Habs = abs.(H)

    # Check if maximum and minimum are the same or if the coefficient of variation is too low
    if maximum(Habs) == minimum(Habs) || std(Habs[Habs .> 0]) < eltype(freq)(0.1)*mean(Habs[Habs .> 0])
        return Complex{eltype(freq)}[]
    end

    # Estimation of natural frequencies and damping ratios
    fn = similar(freq[pks.indices])
    ξn = similar(fn)

    nfreq_itp = 500
    ReH_itp = similar(fn, nfreq_itp)
    ImH_itp = similar(ReH_itp)
    freq_itp = similar(ReH_itp)
    α = similar(ReH_itp)
    θ = similar(ReH_itp)
    for (n, edg) in enumerate(pks.edges)
        # Frequency range around the peak
        edge1 = floor(Int, edg[1])
        edge2 = ceil(Int, edg[2])

        freqs = freq[edge1:edge2]

        # Real and imaginary parts of the frequency response function
        ReH = real(H[edge1:edge2])
        ImH = imag(H[edge1:edge2])

        # Circle fitting
        itp_real = linear_interpolation(freqs, ReH)
        itp_imag = linear_interpolation(freqs, ImH)

        freq_itp .= LinRange(freqs[1], freqs[end], nfreq_itp)
        ReH_itp .= itp_real(freq_itp)
        ImH_itp .= itp_imag(freq_itp)

        # Angles
        @. α = atan(ImH_itp, ReH_itp)
        @. θ = π + 2α

        # Angular frequency such that θ₁ = π/2
        pos1 = argmin(@. abs(θ - π/2))
        f1 = freq_itp[pos1]
        θ1 = θ[pos1]

        # Angular frequency such that θ₂ = -π/2
        pos2 = argmin(@. abs(θ + π/2))
        f2 = freq_itp[pos2]
        θ2 = θ[pos2]

        # Calculation of the natural frequency and damping ratio
        fn[n] = √((f1^2*tan(θ2/2) - f2^2*tan(θ1/2))/(tan(θ2/2) - tan(θ1/2)))
        ξn[n] = (f2^2 - f1^2)/(2fn[n]^2*(tan(θ1/2) - tan(θ2/2)))

        if ξn[n] > 1. || ξn[n] < 0.
            ξn[n] = NaN
        end
    end

    return modal2poles(fn, ξn)
end

function compute_poles(H, freq, alg::LSFit, pks)

    # Flat FRF - No peaks
    Habs = abs.(H)

    # Check if maximum and minimum are the same or if the coefficient of variation is too low
    if maximum(Habs) == minimum(Habs) || std(Habs[Habs .> 0]) < eltype(freq)(0.1)*mean(Habs[Habs .> 0])
        return Complex{eltype(freq)}[]
    end

    # Estimation of natural frequencies and damping ratios
    fn = freq[pks.indices]
    ξn = similar(fn)

    nfreq_itp = 500
    freq_itp = similar(fn, nfreq_itp)
    ReH_itp = similar(freq_itp)
    ImH_itp = similar(freq_itp)
    Hitp = similar(H, nfreq_itp)

    A = similar(Hitp, nfreq_itp, 3)
    b = similar(Hitp, nfreq_itp)
    # Sa = similar(Hitp, 2nfreq_itp, 3)
    # Sb = similar(Hitp, 2nfreq_itp)
    res = similar(fn, 3)
    for (n, edg) in enumerate(pks.edges)
        # Frequency range around the peak
        edge1 = floor(Int, edg[1])
        edge2 = ceil(Int, edg[2])

        freqs = freq[edge1:edge2]

        # Real and imaginary parts of the frequency response function
        ReH = real(H[edge1:edge2])
        ImH = imag(H[edge1:edge2])

        itp_real = linear_interpolation(freqs, ReH)
        itp_imag = linear_interpolation(freqs, ImH)

        freq_itp .= LinRange(freqs[1], freqs[end], nfreq_itp)
        ReH_itp .= itp_real(freq_itp)
        ImH_itp .= itp_imag(freq_itp)
        @. Hitp = ReH_itp + 1im*ImH_itp

        # Construction of Least Squares problem
        A .= [Hitp 2im*freq_itp.*Hitp -one.(Hitp)]
        @. b = freq_itp^2*Hitp
        # Solve the system
        # Tips: Solving the complex system directly can lead to numerical issues
        # so we only use the real part of the system
        # Sa .= [real(A); imag(A)]
        # Sb .= [real(b); imag(b)]
        # res .= qr(Sa)\Sb
        res .= qr(real(A))\real(b)

        # Calculation of the natural frequency and damping ratio
        if real(res[1]) ≤ 0.
            fn[n] = NaN
        else
            fn[n] = √res[1]
        end

        ξn[n] = res[2]/fn[n]

        if ξn[n] > 1. || ξn[n] < 0.
            ξn[n] = NaN
        end
    end

    return modal2poles(fn, ξn)
end

"""
    modeshape_extraction(prob, poles, alg; dpi)
    modeshape_extraction(residues, poles, alg; dpi, modetype)

Extract mode shapes using Sdof approximation

**Inputs**
* `prob::EMAProblem`: EMA problem containing FRF data and frequency vector
* `poles`: Vector of complex poles
* `alg::SdofModalExtraction` or `alg::MdofEMA`: Modal extraction algorithm
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

**References**
[1] M. Géradin and D. J. Rixen. "Mechanical Vibrations: Theory and Application to Structural Dynamics". 3rd. Edition, Wiley, 2015.

[2] C. Ranieri and G. Fabbrocino. "Operational Modal Analysis of Civil Engineering Structures: An Introduction and Guide for Applications". Springer, 2014.
"""
function modeshape_extraction(prob::EMAProblem, poles::Vector{T}, alg::SdofEMA; dpi = [1, 1]) where {T <: Complex}
    # Extract FRF and frequency from problem
    (; frf, freq) = prob
    ω = 2π*freq

    # Conversion of poles to modal parameters
    fn, ξn = poles2modal(poles)
    ωn = 2π*fn

    # Data preparation
    nm, ne = size(frf)[1:2]
    if nm < ne
        # Roving hammer test
        H = frf[dpi[1], :, :]
        dpi = [dpi[2], dpi[1]]
    else
        # Roving accelerometer test
        H = frf[:, dpi[2], :]
    end

    #Initialization
    nx::Int = size(H, 1)
    nfreq = length(fn)
    ϕn = similar(fn, nx, nfreq)

    # Nodes other than driving point
    pos_ndi = findall(x -> x ∉ dpi[1], 1:nx)

    for (n, (ω0, ξ)) in enumerate(zip(ωn, ξn))
        # Trick for estimating the mode shape when ξ = 0
        η = 2ξ
        if η == 0.
            η = 1e-16
        end

        pos_max = argmin(@. abs(ω - ω0))

        # Mode shape estimation at driving point
        Himag = -imag(H[:, pos_max])
        ϕn[dpi[1], n] = √(η*ω0^2*Himag[dpi[1]])

        # Mode shape estimation at other nodes
        for pos in pos_ndi
            ϕn[pos, n] = η*ω0^2*Himag[pos]/ϕn[dpi[1], n]
        end
    end

    return ϕn
end

"""
    circfit(x, y)

Fit a circle to a set of points (x, y) using the least squares method

# Inputs
* `x`: x-coordinates of the points
* `y`: y-coordinates of the points

# Outputs
* `xc`: x-coordinate of the circle center
* `yc`: y-coordinate of the circle center
* `R`: Radius of the circle
"""
function circfit(x, y)
    # System matrix
    S = [x y ones(eltype(x), length(x))]

    # RHS
    b = @. -x^2 - y^2

    # Solve the system
    coeff = (S'S)\(S'b)

    # Extraction of circle parameters
    xc = -coeff[1]/2
    yc = -coeff[2]/2
    R = √(xc^2 + yc^2 - coeff[end])

    return xc, yc, R
end

"""
    solve(prob_ema::AutoEMASdofProblem)
    solve(prob_ema::AutoEMAMdofProblem, order::Int; stabdiag, weighting)

Solve automatically experimental modal analysis problem using Sdof or Mdof methods

**Inputs**
* `prob_ema`: Structure containing the input data for automatic experimental modal analysis using Sdof methods

**Outputs**
* `sol::EMASolution`: Structure containing the solution of the automatic experimental modal analysis using Sdof or Mdof methods
"""
function solve(prob_ema::AutoEMASdofProblem)
    # Unpack problem data
    (; prob, alg, dpi, idx_m, idx_e, width, min_prom, max_prom, pks_indices) = prob_ema

    # Extraction of natural frequencies and damping ratios
    poles = poles_extraction(prob, alg, width = width, min_prom = min_prom, max_prom = max_prom, pks_indices = pks_indices)

    # Extraction of mode shapes
    phi = modeshape_extraction(prob, poles, alg, dpi = dpi)

    res, ci = mode2residues(phi, poles, idx_m, idx_e)

    lr, ur = compute_residuals(prob, res, poles)

    return EMASolution(poles, phi, ci, res, lr, ur)
end