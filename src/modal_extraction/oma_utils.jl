"""
    xcorr(y, N)

Compute the cross-correlation matrix R of a signal y up to lag N-1.

**Inputs**
- `y::VecOrMat{Real}`: Input signal
- `N`: Number of lags

**Output**
- `R`: Cross-correlation matrix of size (ny, ny, N), where ny
         is the number of signals (ny = 1 if y is a vector)


"""
@views function xcorr(y, N)
    ndy = ndims(y)

    if ndy == 1
        ny = 1
        nt = length(y)
    else
        ny, nt = size(y)
    end

    R = similar(y, ny, ny, N)

    for k in 1:N
        if ndy == 1
            R[:, :, k] .= (y[1:nt - k + 1] ⋅ y[k:nt]) / (nt - k)
        else
            R[:, :, k] .= (y[:, 1:nt - k + 1] * y[:, k:nt]') / (nt - k)
        end
    end

    return R
end

"""
    psd_from_tf(H::Array{C, 3}, Gxx; id_input = 1:size(H, 1), id_output = id_input) where {C <: Complex}

Compute the output power spectral density matrix Gyy given the
frequency response function matrix H and the input power spectral density
matrix Gxx.

**Inputs**
- `H::Array{Complex, 3}`: Frequency response function matrix of size (no, ni, nf),
       where ni is the number of inputs, no is the number of outputs,
       and nf is the number of frequency points.
- `Gxx`: Input power spectral density matrix of size (ni, ni, nf)
- `id_input::AbstractVector`: Indices of input channels to consider (default: all inputs)
- `id_output::AbstractVector`: Indices of output channels to consider (default: `id_input`)

**Output**
- `Gyy::Array{Complex, 3}`: Output power spectral density matrix of size (ni, no, nf)

**References**

[1] B. Peeters and H. Van der Auweraer, "PolyMax: A revolution in operational modal analysis". In Proceedings of the 1st International Operational Modal Analysis Conference (IOMAC), Copenhagen, Denmark, 2005.
"""
@views function psd_from_tf(H::Array{C, 3}, Gxx; id_input = 1:size(H, 1), id_output = id_input) where {C <: Complex}
    ni = length(id_input)
    no = length(id_output)
    nf = size(H, 3)
    Gyy = similar(H, no, ni, nf)

    ndg = ndims(Gxx)
    if ndg == 3
        nx = size(Gxx, 1)
        Gxx_f = similar(Gxx, nx, nx)
    end

    for f in 1:nf
        Hi = H[id_input, :, f]
        Ho = H[id_output, :, f]

        if ndg == 3
            Gxx_f .= Gxx[:, :, f]
            Gyy[:, :, f] .= Hi * Gxx_f * Ho'
        else
            Gyy[:, :, f] .= Hi * Gxx * Ho'
        end
    end

    return Gyy
end

"""
    half_psd(Gyy, freq)
    half_psd(y, freq, fs, bs)

    Compute the half power spectral density matrix Gyy_half given the
    spectral density matrix Gyy.

**Inputs**
- `Gyy`: Power spectral density matrix of size (ni, no, nf),
       where ni is the number of inputs, no is the number of outputs,
       and nf is the number of frequency points.
- `freq::AbstractVector`: Frequency vector (Hz)

**Alternative inputs**
- `y::VecOrMat{Real}`: Output signal
- `freq::AbstractVector`: Frequency vector (Hz)
- `fs`: Sampling frequency (Hz)
- `bs`: Block size used to compute the cross-correlation matrix

**Output**
- `Gyy_half::Array{Complex, 3}`: Half power spectral density matrix

**References**
[1] S. Chauhan. "Parameter estimation algorithms in operational modal analysis". In Proceedings of the 6th International Operational Modal Analysis Conference (IOMAC), Gijon, Spain, 2015.

[2] B. Peeters and H. Van der Auweraer, "PolyMax: A revolution in operational modal analysis". In Proceedings of the 1st International Operational Modal Analysis Conference (IOMAC), Copenhagen, Denmark, 2005.

[3] R. Brincker and C. E. Ventura. "Introduction to operational modal analysis". Wiley, 2015.
"""
@views function half_psd(Gyy, freq::AbstractVector)
    # Initialization
    ndg = ndims(Gyy)
    if ndg == 3
        no, ni, nf = size(Gyy)
    elseif ndg == 1
        no = 1
        ni = 1
        nf = length(Gyy)
    else
        throw(ArgumentError("Gyy must be a 1D or 3D array"))
    end

    df = freq[2] - freq[1]                # Frequency resolution
    fs = 2.56*(freq[end] - freq[1])       # Effective sampling frequency
    # fs = 2.56*freq[end]
    freq_g = (0.:df:(fs - df)) .+ freq[1] # Frequency vector for Gyy_half

    Gyy_half = similar(Gyy, no, ni, length(freq_g))
    Gyy_ij = similar(Gyy, nf)

    # Number of points of IFFT
    N = Int(round(fs/df))
    Ryy = similar(freq, N)

    # Compute half power spectral density matrix
    n_half = iseven(N) ? N ÷ 2 + 1 : (N + 1) ÷ 2
    win = flattri(n_half, 0.5)
    for j in 1:ni
        for i in 1:no
            # Step 1: IFFT to get the correlation function
            if ndg == 3
                Gyy_ij .= Gyy[i, j, :]
            else
                Gyy_ij .= Gyy
            end

            Ryy .= impulse_response(Gyy_ij, freq, fs)

            # Step 2: FFT to get the half power spectral density
            # Step 2.1: Keep only positive time lags
            # Step 2.2: Apply windowing
            # Step 2.3: Zero-pad to N points and FFT
            Gyy_half[i, j, :] .= fft([win.*Ryy[1:n_half]; zeros(N - n_half)])
        end
    end

    fidx = freq[1] .<= freq_g .<= freq[end]

    return Gyy_half[:, :, fidx]
end

@views function half_psd(y, freq, fs, bs)
    # Compute cross-correlation matrices - Take only positive lags
    n_half = iseven(bs) ? bs ÷ 2 + 1 : (bs + 1) ÷ 2
    Ryy_pos = xcorr(y, n_half)

    # Compute half power spectral density matrix
    df = freq[2] - freq[1] # Frequency resolution
    freq_g = (0.:df:(fs - df)) # Frequency vector for Gyy_half
    if ndims(y) == 1
        ni = 1
        no = 1
    else
        no = size(y, 1)
        ni = no
    end

    Gyy_half = similar(complex.(Ryy_pos), no, ni, length(freq_g))
    Ryy_ij = similar(Ryy_pos, n_half)

    # Windowing
    win = flattri(n_half, 0.5)

    Npad = bs - n_half
    for j in 1:ni
        for i in 1:no
            # Step 1: FFT to get the half power spectral density
            Ryy_ij .= Ryy_pos[i, j, :]

            # Step 1.1: Apply windowing
            Ryy_ij .*= win

            # Step 1.2: Zero-pad to N points and FFT
            Gyy_half[i, j, :] .= fft([Ryy_ij; zeros(Npad)])
        end
    end

    fidx = freq[1] .<= freq_g .<= freq[end]

    return Gyy_half[:, :, fidx]/n_half
end

"""
    convert_Gyy(Gyy, freq; type = :dis)

Convert a power spectral density matrix Gyy to displacement type

**Inputs**
- `Gyy`: Power spectral density matrix of size (ni, no, nf),
         where ni is the number of inputs, no is the number of outputs,
         and nf is the number of frequency points.
- `freq::AbstractVector`: Frequency vector (Hz)
- `type::Symbol`: Type of Gyy
    - `:dis`: displacement (default) - no conversion
    - `:vel`: velocity
    - `:acc`: acceleration
**Output**
- `Gyy_conv`: Converted power spectral density matrix
"""
function convert_Gyy(Gyy, freq; type = :dis)
    ndg = ndims(Gyy)
    if ndg == 3
        no, ni, nf = size(Gyy)
    elseif ndg == 1
        no = 1
        ni = 1
        nf = length(Gyy)
        Gyy = reshape(Gyy,  (no, ni, length(Gyy)))
    else
        throw(ArgumentError("Gyy must be a 1D or 3D array"))
    end

    Gyy_conv = copy(Gyy)
    ω = 2π*freq
    for (f, ωf) in enumerate(ω)
        # Avoid division by zero
        ωf = ωf == 0. ? ω[2] : ωf
        if type == :vel
            Gyy_conv[:, :, f] /= -ωf^2
        elseif type == :acc
            Gyy_conv[:, :, f] /= ωf^4
        end
    end

    return ndg == 3 ? Gyy_conv : vec(Gyy_conv)
end

"""
    block_hankel(y, nbr)

Constructs the past and future block Hankel matrices for SSI-DATA.

**Inputs**
- `y::Matrix{Float64}`: Matrix of measured outputs (channels × samples)
- `nbr::Int`: Number of block rows

**Outputs**
- `Hp::Matrix{Float64}`: Past block Hankel matrix
- `Hf::Matrix{Float64}`: Future block Hankel matrix

**References**
[1] C. Rainieri and G. Fabbrocino. "Operational Modal Analysis of Civil Engineering Structures: An Introduction and Guide for Applications". Springer, 2014.
"""
function block_hankel(y::Matrix{Float64}, nbr::Int)
    nm, nt = size(y)
    nc = nt - 2nbr + 1 # Number of columns of Hankel matrices

    H = similar(y, nm*2nbr, nc) # Full Hankel matrix
    for i = 1:2nbr
        H[(i-1)*nm+1:i*nm, :] .= y[:, i:i+nc-1]/sqrt(nc)
    end

    return H
end

"""
    block_toeplitz(R, nbr)

Constructs a block Toeplitz matrix from correlation matrices.

**Inputs**
- `R::Array{Float64, 3}`: Correlation matrices of size (
- `nbr::Int`: Number of block rows

**Output**
- `T::Matrix{Float64}`: Block Toeplitz matrix
"""
function block_toeplitz(R::Array{Float64, 3}, nbr::Int)
    nm = size(R, 1)
    T = zeros(nm*nbr, nm*nbr)
    for j in 1:nbr
        for i in 1:nbr
            T[(i-1)*nm+1:i*nm, (j-1)*nm+1:j*nm] .= R[:, :, i + j - 1]
        end
    end

    return T
end