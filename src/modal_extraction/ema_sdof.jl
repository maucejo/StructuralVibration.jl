abstract type SdofModalExtraction end
struct PeakPicking <: SdofModalExtraction end
struct CircleFit <: SdofModalExtraction end
struct LSFit <: SdofModalExtraction end

"""
    EMASdofProblem(frf, freq; frange = [freq[1], freq[end]])

Data structure defining the inputs for EMA-MDOF modal extraction methods.

**Constructor parameters**
- `frf::Array{Complex,3}`: 3D array of Frequency Response Functions (FRF) (array nm x ne x nf)
- `freq::AbstractArray{Real}`: Vector of frequency values (Hz)
- `frange::Vector{Real}`: Frequency range for analysis (default: [freq[1], freq[end]])

**Fields**
- `frf::Array{Complex, 3}`: 3D array of Frequency Response Functions (FRF) (array nm x ne x nf)
- `freq::AbstractArray{Real}`: Vector of frequency values (Hz)
"""
@show_data struct EMASdofProblem{C <: Complex, R <: Real}
    frf::Array{C, 3}
    freq::AbstractArray{R}

    function EMASdofProblem(frf::Array{C,3}, freq::AbstractArray{R}; frange = [freq[1], freq[end]]) where {C <: Complex, R <: Real}

        # FRF post-processing - Frequency range reduction
        fidx = @. frange[1] ≤ freq ≤ frange[2]

        return new{C, R}(frf[:, :, fidx], freq[fidx])
    end
end

"""
    AutoEMASdofProblem(prob, dpi, method; type_frf)

Structure containing the input data for automatic experimental modal analysis using Sdof methods

**Fields**
* `prob::EMASdofProblem`: EMA-SDOF problem containing FRF data and frequency vector
* `dpi::Vector{Int}`: Driving point indices - default = [1, 1]
    * `dpi[1]`: Driving point index on the measurement mesh
    * `dpi[2]`: Driving point index on the excitation mesh
* `method::SdofModalExtraction`: Method to extract the poles
    * `PeakPicking`: Peak picking method (default)
    * `CircleFit`: Circle fitting method
    * `LSFit`: Least squares fitting method
* `type_frf::Symbol`: Type of FRF used for poles extraction
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance
"""
@show_data struct AutoEMASdofProblem
    prob::EMASdofProblem
    dpi:: Vector{Int}
    method::SdofModalExtraction
    type_frf::Symbol

    AutoEMASdofProblem(prob::EMASdofProblem, dpi::Vector{Int} = [1, 1], method::SdofModalExtraction = PeakPicking(); type_frf::Symbol = :dis) = new(prob, dpi, method, type_frf)
end

"""
    EMASdofSolution(poles, ϕn)

Structure containing the solution of the automatic experimental modal analysis using Sdof methods

**Fields**
* `poles`: Extracted poles
* `ms`: Mode shapes
"""
@show_data struct EMASdofSolution{Tp <: Complex, Tm <: Real}
    poles::Vector{Tp}
    ms::AbstractArray{Tm}
end

"""
    poles_extraction(prob, method::SdofModalExtraction; type = :dis)

Extract poles from the Bode diagram fitting method

**Inputs**
* `prob::EMASdofProblem`: EMA-SDOF problem containing FRF data and frequency vector
* `method::SdofModalExtraction`: Method to extract the poles
    * `PeakPicking`: Peak picking method (default)
    * `CircleFit`: Circle fitting method
    * `LSFit`: Least squares fitting method
* `frange`: Frequency range to consider for the extraction - default = [freq[1], freq[end]]
* `type`: Type of FRF used to extract the poles
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

**Outputs**
* `poles`: Vector of extracted complex poles

**Note**
The natural frequencies and damping ratios are extracted from each FRF (each row of the matrix) and then averaged. The number of FRF used for averaging are those having the maximum (and same) number of peaks detected.
"""
function poles_extraction(prob::EMASdofProblem, method::SdofModalExtraction;  type = :dis)

    # FRF post-processing - Frequency range reduction
    # Extract FRF and frequency from problem
    (; frf, freq) = prob

    nm, ne, nf = size(frf)
    Hr = reshape(frf, nm*ne, nf)

    # Estimate the number of peaks in the FRF
    npeak = 0
    np = zeros(Int, nm*ne)
    for (k, Hv) in enumerate(eachrow(Hr))
        np[k] = length(findmaxima(abs.(Hv)).indices)
        npeak = max(npeak, np[k])
    end

    est_method = let
        if method isa PeakPicking
            poles_ppm_extract
        elseif method isa CircleFit
            poles_cfm_extract
        elseif method isa LSFit
            poles_lsf_extract
        else
            throw(ArgumentError("Unknown modal extraction method"))
        end
    end

    # Initialization
    pk = similar(frf, nm*ne, npeak)
    pn = similar(frf, npeak)
    for (k, Hv) in enumerate(eachrow(Hr))
        p = est_method(Hv, freq, type = type)

        nk = length(p)
        pk[k, :] .= [p; fill(complex(NaN, NaN), npeak - nk)]
    end

    # Average the results from different FRFs
    keeprow = np .== npeak
    for i in 1:npeak
        pn[i] = mean(skipnan(pk[keeprow, i]))
    end

    return pn
end

"""
    poles_ppm_extract(H, freq; type = :dis)

Extract poles from the Bode diagram peak picking method

**Inputs**
* `H`: Frequency response function
* `freq`: Frequency vector
* `type`: Type of FRF used to extract the poles
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

**Outputs**
* `poles`: Extracted poles
"""
function poles_ppm_extract(H, freq; type = :dis)
    ω = 2π*freq

    Habs = abs.(H)
    pks = findmaxima(Habs)
    # Peaks not found
    if isempty(pks.indices)
        return Complex{eltype(freq)}[]
    end
    pks = peakproms!(pks) |> peakwidths!

    # Natural frequencies
    ftemp = freq[pks.indices]
    fn = similar(ftemp)
    ξn = similar(ftemp)

    # Convert to displacement
    if type == :vel
        H ./= 1im*ω
    elseif type == :acc
        H ./= -ω.^2
    end

    # Damping ratios estimation from the half-bandwidth method
    nfreq_itp = 250
    freq_left = similar(fn, nfreq_itp)
    freq_right = similar(freq_left)
    Hleft = similar(freq_left)
    Hright = similar(freq_left)
    for (n, (f, idmax, Hmax, edg)) in enumerate(zip(ftemp, pks.indices, pks.heights, pks.edges))
        edge1, edge2 = round.(Int, edg) .+ [-5, 5]

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

        # Correction of the natural frequency due to damping
        fn[n] = ftemp[n]/√(1 - ξn[n]^2)
    end

    return modal2poles(fn, ξn)
end

"""
    poles_cfm_extract(H, freq; type = :dis)

Extract poles from the Bode diagram circle fitting method

**Inputs**
* `H`: Frequency response function
* `freq`: Frequency vector
* `type`: Type of FRF used to extract the poles
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

**Outputs**
* `poles`: Extracted poles
"""
function poles_cfm_extract(H, freq; type = :dis)
    ω = 2π*freq

    pks = findmaxima(abs.(H))
    # Peaks not found
    if isempty(pks.indices)
        return Complex{eltype(freq)}[]
    end
    pks = peakproms!(pks) |> peakwidths!

    # Convert to displacement
    if type == :vel
        H ./= 1im*ω
    elseif type == :acc
        H ./= -ω.^2
    end

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
        edge1, edge2 = round.(Int, edg) .+ [-5, 5]
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
        f₁ = freq_itp[pos1]
        θ₁ = θ[pos1]

        # Angular frequency such that θ₂ = -π/2
        pos2 = argmin(@. abs(θ + π/2))
        f₂ = freq_itp[pos2]
        θ₂ = θ[pos2]

        # Calculation of the natural frequency and damping ratio
        fn[n] = √((f₁^2*tan(θ₂/2) - f₂^2*tan(θ₁/2))/(tan(θ₂/2) - tan(θ₁/2)))
        ξn[n] = (f₂^2 - f₁^2)/(2fn[n]^2*(tan(θ₁/2) - tan(θ₂/2)))

    end

    return modal2poles(fn, ξn)
end

"""
    poles_lsf_extract(H, freq; type = :dis)

Extract poles from the Bode diagram least squares fitting method

**Inputs**
* `H`: Frequency response function
* `freq`: Frequency vector
* `type`: Type of FRF used to extract the poles
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

**Outputs**
* `poles`: Extracted poles

**Source**
[1] Matlab documentation - `modalfit` function (https://www.mathworks.com/help/signal/ref/modalfit.html)
[2] A. Brandt, "Noise and Vibration Analysis: Signal Analysis and Experimental Procedures", Wiley, 2011.

"""
function poles_lsf_extract(H, freq; type = :dis)
    ω = 2π*freq

    pks = findmaxima(abs.(H))
    # Peaks not found
    if isempty(pks.indices)
        return Complex{eltype(freq)}[]
    end
    pks = peakproms!(pks) |> peakwidths!

    # Natural frequencies
    fn = freq[pks.indices]
    ξn = similar(fn)

    # Convert to displacement
    if type == :vel
        H ./= 1im*ω
    elseif type == :acc
        H ./= -ω.^2
    end

    nfreq_itp = 500
    freq_itp = similar(fn, nfreq_itp)
    ReH_itp = similar(freq_itp)
    ImH_itp = similar(freq_itp)
    Hitp = similar(H, nfreq_itp)

    A = similar(Hitp, nfreq_itp, 3)
    b = similar(Hitp, nfreq_itp)
    Sa = similar(Hitp, 2nfreq_itp, 3)
    Sb = similar(Hitp, 2nfreq_itp)
    res = similar(b, 3)
     for (n, edg) in enumerate(pks.edges)
        # Frequency range around the peak
        edge1, edge2 = round.(Int, edg) .+ [-5, 5]
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
        # so we separate the real and imaginary parts to solve a real system of double size using
        Sa .= [real(A); imag(A)]
        Sb .= [real(b); imag(b)]
        res .= (Sa'Sa)\(Sa'Sb)
        # Actually, only one of the two parts can be used to solve the system
        # res .= real(A)\real(b)

        # Calculation of the natural frequency and damping ratio
        fn[n] = √res[1]
        ξn[n] = res[2]/fn[n]
    end

    return modal2poles(fn, ξn)
end

"""
    modeshape_extraction(prob, poles, dpi; type = :dis)

Extract mode shapes using Sdof approximation

**Inputs**
* `prob::EMASdofProblem`: EMA-SDOF problem containing FRF data and frequency vector
* `poles`: Vector of complex poles
* `dpi`: Driving point indices - default = [1, 1]
    * `dpi[1]`: Driving point index on the measurement mesh
    * `dpi[2]`: Driving point index on the excitation mesh
* `type`: Type of FRF used to extract the natural frequencies and damping ratios
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

**Output**
* `ϕn`: Mode shapes

**Note**
If the number of measurement points is less than the number of excitation points, the mode shapes are estimated at the excitation points (roving hammer test). Otherwise, the mode shapes are estimated at the measurement points (roving accelerometer test).
"""
function modeshape_extraction(prob::EMASdofProblem, poles::Vector{T}, dpi = [1, 1]; type = :dis) where {T <: Complex}
    # Extract FRF and frequency from problem
    (; frf, freq) = prob

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

    # Convert to admittance
    ω = 2π*freq
    if type == :vel
        H ./= ω'
    elseif type == :acc
        H ./= ω'.^2
    end

    #Initialization
    nx::Int = size(H, 1)
    nfreq = length(fn)
    ϕn = similar(fn, nx, nfreq)

    # Nodes other than driving point
    pos_ndi = findall(x -> x ∉ dpi[1], 1:nx)

    for (n, (ω₀, ξ)) in enumerate(zip(ωn, ξn))
        # Trick for estimating the mode shape when ξ = 0
        η = 2ξ
        if η == 0.
            η = 1e-16
        end

        pos_max = argmin(@. abs(ω - ω₀))

        # Mode shape estimation at driving point
        Himag = -imag(H[:, pos_max])
        ϕn[dpi[1], n] = √(η*ω₀^2*Himag[dpi[1]])

        # Mode shape estimation at other nodes
        for pos in pos_ndi
            ϕn[pos, n] = η*ω₀^2*Himag[pos]/ϕn[dpi[1], n]
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
    solve(prob::AutoEMASdofProblem)

Solve automatically experimental modal analysis problem using Sdof methods

**Inputs**
* `prob_ema`: Structure containing the input data for automatic experimental modal analysis using Sdof methods

**Outputs**
* `sol`: Solution of problem:
    * `poles`: Extracted poles
    * `ms`: Mode shapes
"""
function solve(prob_ema::AutoEMASdofProblem)
    # Unpack problem data
    (; prob, dpi, method, type_frf) = prob_ema

    # Extraction of natural frequencies and damping ratios
    poles = poles_extraction(prob, method, type = type_frf)

    # Extraction of mode shapes
    phi = modeshape_extraction(prob, poles, dpi, type = type_frf)

    return EMASdofSolution(poles, phi)
end