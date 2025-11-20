"""
    xcorr(y, yref, N)
    xcorr(y, N)

Compute the cross-correlation matrix R of a signal y up to lag N-1.

**Inputs**
- `y::VecOrMat{Real}`: Input signal
- `yref::VecOrMat{Real}`: Reference input signal
- `N`: Number of lags

**Output**
- `R`: Cross-correlation matrix of size (ny, nyref, N), where ny
         is the number of signals (ny = 1 if y is a vector)

**Note**
- xcorr(y, N) = xcorr(y, y, N)
"""
@views function xcorr(y, yref, N)
    # Work with proper matrices
    if y isa Vector
        y = reshape(y, 1, :)
    end

    if yref isa Vector
        yref = reshape(yref, 1, :)
    end

    # Check sizes
    ny, nt = size(y)
    nyref, ntref = size(yref)

    nt != ntref && throw(ArgumentError("Input signals must have the same number of samples"))

    # Compute cross-correlation matrices
    R = similar(y, ny, nyref, N)
    for k in 1:N
        R[:, :, k] .= (y[:, 1:nt - k + 1] * yref[:, k:nt]') / (nt - k)
    end

    return R
end

xcorr(y, N) = xcorr(y, y, N)

"""
    csd_from_tf(H, Gxx; id_input, id_output)

Compute the output cross-spectral density matrix Gyy given the
frequency response function matrix H and the input cross-spectral density
matrix Gxx.

**Inputs**
- `H::Array{Complex, 3}`: Frequency response function matrix of size (no, ni, nf),
- `Gxx`: Input power spectral density matrix of size (ni, ni, nf)
- `id_input::AbstractVector`: Indices of input channels to consider (default: all inputs)
- `id_output::AbstractVector`: Indices of output channels to consider (default: `id_input`)

**Output**
- `Gyy::Array{Complex, 3}`: Output cross-spectral density matrix of size (ni, no, nf)

**References**

[1] B. Peeters and H. Van der Auweraer, "PolyMax: A revolution in operational modal analysis". In Proceedings of the 1st International Operational Modal Analysis Conference (IOMAC), Copenhagen, Denmark, 2005.
"""
@views function csd_from_tf(H::Array{C, 3}, Gxx; id_input = 1:size(H, 1), id_output = id_input) where {C <: Complex}
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
            Gyy[:, :, f] .= Ho * Gxx_f * Hi'
        else
            Gyy[:, :, f] .= Ho * Gxx * Hi'
        end
    end

    return Gyy
end

"""
    half_psd(Gyy, freq)
    half_psd(y, yref, freq, fs, bs)
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
- `yref::VecOrMat{Real}`: Reference output signal
- `freq::AbstractVector`: Frequency vector (Hz)
- `fs`: Sampling frequency (Hz)
- `bs`: Block size used to compute the cross-correlation matrix

**Output**
- `Gyy_half::Array{Complex, 3}`: Half-spectrum matrix

**References**

[1] S. Chauhan. "Parameter estimation algorithms in operational modal analysis". In Proceedings of the 6th International Operational Modal Analysis Conference (IOMAC), Gijon, Spain, 2015.

[2] B. Peeters and H. Van der Auweraer, "PolyMax: A revolution in operational modal analysis". In Proceedings of the 1st International Operational Modal Analysis Conference (IOMAC), Copenhagen, Denmark, 2005.

[3] R. Brincker and C. E. Ventura. "Introduction to operational modal analysis". Wiley, 2015.
"""
@views function half_csd(Gyy, freq::AbstractVector)
    # Initialization
    if Gyy isa Vector
        Gyy = reshape(Gyy, 1, 1, length(Gyy))
    end

    no, ni, nf = size(Gyy)

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
            Gyy_ij .= Gyy[i, j, :]

            Ryy .= impulse_response(Gyy_ij, freq, fs)

            # Step 2: FFT to get the half power spectral density
            # Step 2.1: Keep only positive time lags
            # Step 2.2: Apply windowing
            # Step 2.3: Zero-pad to N points and FFT
            Gyy_half[i, j, :] .= fft([win.*Ryy[1:n_half]; zeros(N - n_half)])
        end
    end

    fidx = freq[1] .<= freq_g .<= freq[end]

    return Gyy_half[:, :, fidx]/rms(win)
end

@views function half_csd(y, yref, freq, fs, bs)
    if y isa Vector
        y = reshape(y, 1, :)
    end

    if yref isa Vector
        yref = reshape(yref, 1, :)
    end

    no = size(y, 1)
    ni = size(yref, 1)

    # Compute cross-correlation matrices - Take only positive lags
    nblocks = size(y, 2) ÷ bs
    n_half = iseven(bs) ? bs ÷ 2 + 1 : (bs + 1) ÷ 2

    # Compute half power spectral density matrix
    df = freq[2] - freq[1] # Frequency resolution
    freq_g = (0.:df:(fs - df)) # Frequency vector for Gyy_half

    # Windowing
    win = flattri(n_half, 0.5)

    Npad = bs - n_half
    Ryy_block = similar(y, no, ni, n_half)
    Gyy_half = zeros(Complex{eltype(y)}, no, ni, length(freq_g))
    Ryy_ij = similar(Ryy_block, n_half)
    for b in 1:nblocks
        Ryy_block .= xcorr(y[:, (b-1)*bs+1:b*bs], yref[:, (b-1)*bs+1:b*bs], n_half)
        for j in 1:ni
            for i in 1:no
                # Step 1: FFT to get the half power spectral density
                Ryy_ij .= Ryy_block[i, j, :]

                # Step 1.1: Apply windowing
                Ryy_ij .*= win

                # Step 1.2: Zero-pad to N points and FFT
                Gyy_half[i, j, :] .+= fft([Ryy_ij; zeros(Npad)])
            end
        end
    end

    fidx = freq[1] .<= freq_g .<= freq[end]

    return Gyy_half[:, :, fidx]/(n_half*nblocks*rms(win))
end

half_csd(y, freq, fs, bs) = half_csd(y, y, freq, fs, bs)

"""
    convert_csd(Gyy, freq; type = :dis)

Convert a cross-spectral density matrix to displacement type

**Inputs**
- `Gyy`: Cross-spectral density matrix of size (ni, no, nf),
         where ni is the number of inputs, no is the number of outputs,
         and nf is the number of frequency points.
- `freq::AbstractVector`: Frequency vector (Hz)
- `type::Symbol`: Type of Gyy
    - `:dis`: displacement (default) - no conversion
    - `:vel`: velocity
    - `:acc`: acceleration
**Output**
- `Gyy_conv`: Converted cross-spectral density matrix
"""
function convert_csd(Gyy, freq; type = :dis)
    if Gyy isa Vector
        Gyy = reshape(Gyy, 1, 1, length(Gyy))
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

    return Gyy_conv
end

"""
    block_hankel(y, yref, nbr)

Constructs the past and future block Hankel matrices for SSI-DATA.

**Inputs**
- `y::Matrix{Float64}`: Matrix of measured outputs (channels × samples)
- `yref::Matrix{Float64}`: Matrix of reference outputs (ref channels × samples)
- `nbr::Int`: Number of block rows

**Outputs**
- `Hp::Matrix{Float64}`: Past block Hankel matrix
- `Hf::Matrix{Float64}`: Future block Hankel matrix

**References**
[1] C. Rainieri and G. Fabbrocino. "Operational Modal Analysis of Civil Engineering Structures: An Introduction and Guide for Applications". Springer, 2014.
"""
function block_hankel(y::Matrix{Float64}, yref::Matrix{Float64}, nbr::Int)
    no, nt = size(y)
    nref, ntref = size(yref)

    nt ≠ ntref && throw(ArgumentError("Input signals must have the same number of samples."))

    nc = nt - 2nbr + 1 # Number of columns of Hankel matrices

    Yp = similar(y, nref*nbr, nc) # Past-Hankel matrix
    Yf = similar(y, no*nbr, nc)   # Future-Hankel matrix
    for i = 1:nbr
        Yp[(i-1)*nref+1:i*nref, :] .= yref[:, i:i+nc-1]/sqrt(nc)
        Yf[(i-1)*no+1:i*no, :] .= y[:, i+nbr:i+nbr+nc-1]/sqrt(nc)
    end

    return Yp, Yf
end

"""
    block_toeplitz(R, nbr)

Constructs a block Toeplitz matrix from correlation matrices.

**Inputs**
- `R::Array{Float64, 3}`: Correlation matrices
- `nbr::Int`: Number of block rows

**Output**
- `T::Matrix{Float64}`: Block Toeplitz matrix
"""
function block_toeplitz(R::Array{Float64, 3}, nbr::Int)
    no, ni = size(R)[1:2]
    T = zeros(no*nbr, ni*nbr)
    for j in 1:nbr
        for i in 1:nbr
            T[(i-1)*no+1:i*no, (j-1)*ni+1:j*ni] .= R[:, :, i + j - 1]
        end
    end

    return T
end

"""
    oma_modal_parameters(A, C, dt, freq, stabdiag)

Extract modal parameters from system matrices A and C.

**Inputs**
- `A::Matrix{Complex}`: System matrix
- `C::Matrix{Complex}`: Output matrix
- `dt::Real`: Time step
- `freq::AbstractVector{Float64}`: Frequency vector
- `stabdiag::Bool`: Whether to compute stabilization diagram

**Outputs**
- `poles::Vector{Complex}`: Extracted poles
- `ms::Array{Complex, 2}`: Extracted mode shapes
"""
function oma_modal_parameters(A, C, dt, order, freq, stabdiag)
    # Eigen decomposition of system matrix
    eigvals, eigvecs = eigen(A)

    # Extract poles and mode shapes
    p = log.(eigvals)/dt

    # Poles that do not appear as pairs of complex conjugate numbers with a positive real part or that are purely real are suppressed. For complex conjugate poles, only the pole with a positive imaginary part is retained.
    valid_poles = p[@. imag(p) < 0. && real(p) ≤ 0. && !isreal(p)]
    idx_valid_poles = findall(in(conj.(valid_poles)), p)
    p_valid = p[idx_valid_poles]
    idx_sort_poles = sortperm(p_valid, by = abs)
    p_sort = p_valid[idx_sort_poles]

    # Check if poles are in frange
    fn = poles2modal(p_sort)[1]
    fidx = @. freq[1] ≤ fn ≤ freq[end]
    poles = p_sort[fidx]

    # Extract mode shapes
    ms = C*eigvecs[:, idx_valid_poles[idx_sort_poles[fidx]]]

    return stabdiag ? [poles; fill(complex(NaN, NaN), order - length(poles))] : poles, ms
end

"""
    compute_ss_matrix(M, no, order)

Compute system matrices A and C from matrix M.

**Inputs**
- `M::Matrix{Complex}`: Input matrix
- `no::Int`: Number of outputs
- `order::Int`: Order of the system

**Outputs**
- `A::Matrix{Complex}`: System matrix
- `C::Matrix{Complex}`: Output matrix
"""
function compute_ss_matrix(M, no, order)
    F_svd = svd(M)
    U = F_svd.U[:, 1:order]
    S = Diagonal(F_svd.S[1:order])

    # Observability matrix
    Obs = U*sqrt.(S)

    # Output matrix
    C = Obs[1:no, :]

    # System matrix
    Obs_up = Obs[1:end-no, :]
    Obs_low = Obs[no+1:end, :]
    A = Obs_up\Obs_low

    return A, C
end