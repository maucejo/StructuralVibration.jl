abstract type ModalExtraction end
struct PeakPicking <: ModalExtraction end
struct CircleFit <: ModalExtraction end

"""
    freq_extraction(freq, H, method::ModalExtraction; type = :dis)

Extract natural frequencies and damping ratios from the Bode diagram fitting method

# Inputs
* `freq`: Frequency vector
* `H`: Frequency response function
* `method`: Method to extract the natural frequencies and damping ratios
    * `PeakPickingExtract`: Peak picking method (default)
    * `CircleFitExtract`: Circle fitting method
* `type`: Type of FRF used to extract the natural frequencies and damping ratios
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

# Outputs
* `fn`: Natural frequencies
* `ξn`: Damping ratios
"""
function freq_extraction(freq, H, method::ModalExtraction = PeakPicking(); type = :dis)
    if method isa PeakPicking
        return freq_ppm_extract(freq, H, type = type)
    elseif method isa CircleFit
        return freq_cfm_extract(freq, H, type = type)
    end
end

function freq_ppm_extract(freq, H; type = :dis)
    ω = 2π*freq

    Habs = abs.(H)
    pks = findmaxima(Habs) |> peakproms! |> peakwidths!

    # Natural frequencies
    fn = freq[pks.indices]
    ξn = similar(fn)

    # Convert to displacement
    if type == :vel
        H ./= ω
    elseif type == :acc
        H ./= -ω.^2
    end

    # Damping ratios estimation from the half-bandwidth method
    nlocal = 500
    freq_left = similar(fn, nlocal)
    freq_right = similar(freq_left)
    Hleft = similar(freq_left)
    Hright = similar(freq_left)
    for (n, (f, idmax, Hmax, edg)) in enumerate(zip(fn, pks.indices, pks.heights, pks.edges))
        edge1, edge2 = round.(Int, edg) .+ [-5, 5]

        # Left side of the peak
        itp_left = linear_interpolation(freq[edge1:idmax], Habs[edge1:idmax])
        freq_left .= LinRange(freq[edge1], f, nlocal)
        @. Hleft = itp_left(freq_left)
        posleft = argmin(abs.(Hleft .- Hmax/√2))
        fmin = freq_left[posleft]

        # Right side of the peak
        itp_right = linear_interpolation(freq[idmax:edge2], Habs[idmax:edge2])
        freq_right .= LinRange(f, freq[edge2], nlocal)
        @. Hright = itp_right(freq_right)
        posright = argmin(abs.(Hright .- Hmax/√2))
        fmax = freq_right[posright]

        # A structural damping is supposed η = 2ξ. So ξ = η/2
        ξn[n] = (fmax^2 - fmin^2)/4f^2
    end

    return fn, ξn
end

function freq_cfm_extract(freq, H; type = :dis)
    ω = 2π*freq

    pks = findmaxima(abs.(H)) |> peakproms! |> peakwidths!

    # Convert to displacement
    if type == :vel
        H ./= ω
    elseif type == :acc
        H ./= ω.^2
    end

    fn = similar(pks.indices)
    ξn = similar(fn)

    nfreq_itp = 500
    ReH_itp = similar(fn, nfreq_itp)
    ImH_itp = similar(ReH_itp)
    freq_itp = similar(ReH_itp)
    α = similar(ReH_itp)
    θ = similar(ReH_itp)
    for (n, edg) in enumerate(pks.edges)
        # Frequency range around the peak
        edges = round.(Int, edg) .+ [-5, 5]
        freqs = freq[edges[1]:edges[2]]

        # Real and imaginary parts of the frequency response function
        ReH = real(H[edges[1]:edges[2]])
        ImH = imag(H[edges[1]:edges[2]])

        # Circle fitting
        xc = circfit(ReH, ImH)[1]

        itp_real = linear_interpolation(freqs, ReH)
        itp_imag = linear_interpolation(freqs, ImH)

        freq_itp .= LinRange(freqs[1], freqs[end], nfreq_itp)
        ReH_itp .= itp_real(freq_itp)
        ImH_itp .= itp_imag(freq_itp)

        # Angles
        @. α = atan(ImH_itp, ReH_itp - xc)
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

    return fn, ξn
end

"""
    modeshape_extraction(freq, H, fn, ξn, id_exc, method::ModalExtraction; type = :dis)

Extract mode shapes from the Bode's diagram of Nyquiste's circle fitting method

# Inputs
* `freq`: Frequency vector
* `H`: Frequency response function
* `fn`: Natural frequencies
* `ξn`: Damping ratios
* `id_exc`: Identifier of the excited nodes
* `method`: Method to extract the mode shapes
    * `BodeExtract`: Bode diagram (default)
    * `NyquistExtract`: Nyquist diagram
* `type`: Type of FRF used to extract the natural frequencies and damping ratios
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

# Output
* `ϕn`: Mode shapes
"""
function modeshape_extraction(freq, H, fn, ξn, id_exc, method::ModalExtraction = PeakPicking(); type = :dis)
    if method isa PeakPicking
        return modeshape_ppm_extract(freq, H, fn, ξn, id_exc; type = type)
    elseif method isa CircleFit
        return modeshape_nyquist_extract(freq, H, fn, ξn, id_exc; type = type)
    end
end

function modeshape_ppm_extract(freq, H, fn, ξn, id_exc; type = :dis)
    ω = 2π*freq
    ωn = 2π*fn

    ndim = ndims(H)
    if ndim == 1
        H = transpose(H)
    elseif ndim == 2
        H = H[:, id_exc, :]
    end

    # Convert to displacement
    if type == :vel
        H ./= ω'
    elseif type == :acc
        H ./= ω'.^2
    end

    #Initialization
    nx = size(H, 1)
    nfreq = length(fn)
    ϕn = similar(fn, nx, nfreq)

    # Unexcited nodes
    pos_nexc = findall(x -> x ∉ id_exc, 1:nx)

    for (n, (ω₀, ξ)) in enumerate(zip(ωn, ξn))
        # Trick for estimating the mode shape when ξ = 0
        η = 2ξ
        if η == 0.
            η = 1e-16
        end

        pos_max = argmin(abs.(ω - ω₀))

        # Mode shape estimation at driving point
        Hmax = H[id_exc, pos_max]
        signH = sign(angle(Hmax))
        ϕn[id_exc, n] = signH*√(η*ω₀^2*abs(Hmax))

        # Mode shape estimation at other nodes
        for pos in pos_nexc
            Hmax = H[pos, pos_max]
            signH = sign(angle(Hmax))
            ϕn[pos, n] = signH*η*ω₀^2*abs(Hmax)/abs(ϕn[id_exc, n])
        end
    end

    return ϕn
end

function modeshape_cfm_extract(freq, H, fn, ξn, id_exc; type = :dis)

    ω = 2π*freq
    ωn = 2π*fn

    ndim = ndims(H)
    if ndim == 1
        H = transpose(H)
    elseif ndim == 2
        H = H[:, id_exc, :]
    end

    pks = findmaxima(abs.(H[id_exc, :])) |> peakproms! |> peakwidths!

    # Convert to displacement
    if type == :vel
        H ./= ω'
    elseif type == :acc
        H ./= ω'.^2
    end

    #Initialization
    nx = size(H, 1)
    nfreq = length(fn)
    ϕn = similar(fn, nx, nfreq)

    # Unexcited nodes
    pos_nexc = findall(x -> x ∉ id_exc, 1:nx)

    for (n, (ω₀, ξ, edg)) in enumerate(zip(ωn, ξn, pks.edges))
        η = 2ξ
        edges = round.(Int, edg) .+ [-5, 5]
        freqs = freq[edges[1]:edges[2]]

        # Real and imaginary parts of the frequency response function
        Hd = H[id_exc, edges[1]:edges[2]]
        ReH = real(Hd)
        ImH = imag(Hd)

        # Circle fitting
        R = circfit(ReH, ImH)[end]
        An = 2R*η*ω₀^2

        # Mode shape estimation at driving point
        pos_max = argmin(abs.(ω - ω₀))
        Hmax = Hd[pos_max]
        signH = sign(angle(Hmax))
        ϕn[id_exc, n] = signH*√(An)

        # Mode shape estimation at other nodes
        for pos in pos_nexc
            Hd = H[pos, edges[1]:edges[2]]
            ReH = real(Hd)
            ImH = imag(Hd)

            R = circfit(ReH, ImH)[end]
            An = 2R*η*ω₀^2

            Hmax = Hd[pos_max]
            signH = sign(angle(Hmax))
            ϕn[pos, n] = signH*An/abs(ϕn[id_exc, n])
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