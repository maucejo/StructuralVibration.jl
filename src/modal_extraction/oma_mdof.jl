abstract type OMAModalExtraction end
struct CovSSI <: OMAModalExtraction end
struct DataSSI <: OMAModalExtraction end

"""
    OMAProblem(Gyy, freq; frange, type_data)

Data structure defining the inputs for Operational Modal Analysis (OMA) methods.

**Constructor parameters**
- `Gyy::Array{Complex, 3}`: Power spectral density matrix of size (no, no, nf),
       where no is the number of outputs and nf is the number of frequency points.
- `freq::AbstractArray{Real}`: Frequency vector corresponding to the PSD matrix
- `frange::AbstractVector{Real}`: Frequency range to consider for analysis (default: full range)
- `type_data::Symbol`: Type of measured data
    - `:dis`: Displacement (default)
    - `:vel`: Velocity
    - `:acc`: Acceleration

**Alternative constructor parameters**
- `y::AbstractMatrix{Real}`: Matrix of measured outputs (channels x samples)
- `t::AbstractArray{Real}`: Time vector corresponding to the measurements
- `fs::Real`: Sampling frequency (Hz)
- `bs::Int`: Block size for PSD estimation
- `frange::AbstractVector{Real}`: Frequency range to consider for analysis (default: [0., fs/2.56])
- `type_data::Symbol`: Type of measured data
    - `:dis`: Displacement (default)
    - `:vel`: Velocity
    - `:acc`: Acceleration

**Fields**
- `y::AbstractMatrix{Real}`: Matrix of measured outputs (channels x samples)
- `t::AbstractArray{Real}`: Time vector corresponding to the measurements
- `halfspec::Array{Complex, 3}`: Half power spectral density matrix
- `freq::AbstractArray{Real}`: Frequency vector corresponding to the half-spectrum
- `type_data::Symbol`: Type of measured data (:dis, :vel, :acc
"""
@show_data struct OMAProblem{C <: Complex, R <: Real} <: MdofProblem
    y::Matrix{R}
    t::AbstractArray{R}
    halfspec::Array{C, 3}
    freq::AbstractArray{R}
    type_data::Symbol

    function OMAProblem(Gyy::Array{C, 3}, freq::AbstractArray{R}; frange = [freq[1], freq[end]], type_data = :dis) where {C <: Complex, R <: Real}

        # Correct frange to avoid division by zero
        frange[1] == 0. ? frange[1] = 1. : nothing

        # FRF post-processing - Frequency range reduction
        fidx = @. frange[1] ≤ freq ≤ frange[2]
        freq_red = freq[fidx]

        ω = 2π*freq_red
        Gyy_red = Gyy[:, :, fidx]

        # Conversion to admittance
        for (f, ωf) in enumerate(ω)
            if type_data == :vel
                Gyy_red[:, :, f] ./= -ωf^2
            elseif type_data == :acc
                Gyy_red[:, :, f] ./= ωf^4
            end
        end

        # Compute half power spectral density
        halfspec = half_psd(Gyy_red, freq_red)

        return new{C, R}(Matrix{typeof(freq)(undef, 0, 0)}, similar(freq, 0), halfspec, freq_red, type_data)
    end

    function OMAProblem(y::Matrix{R}, t::AbstractArray{R}, fs, bs; frange = [0., fs/2.56], type_data = :dis) where {R <: Real}

        # Compute half power spectral density matrix
        df = fs/bs # Frequency resolution
        freq = (0.:df:(fs - df))

        # Correct frange to avoid division by zero
        frange[1] == 0. ? frange[1] = 1. : nothing

        # FRF post-processing - Frequency range reduction
        fidx = @. frange[1] ≤ freq ≤ frange[2]
        freq_red = freq[fidx]

        halfspec = half_psd(y, freq_red, fs, bs)

        return new{Complex{R}, R}(y, t, halfspec, freq_red, type_data)
    end
end


"""
    poles_extraction(prob, order, method; stabdiag)

Extract poles from half-spectrum data using Operational Modal Analysis (OMA) methods.

**Inputs**
- `prob::OMAProblem`: Structure containing half-spectrum data and frequency vector
- `order::Int`: Order of the system to identify
- `method::OMAModalExtraction`: OMA method to use for pole extraction
    - `CovSSI`: Covariance-based SSI (default)
    - `DataSSI`: Data-based SSI
- `stabdiag::Bool`: Whether to compute stabilization diagram (default: false)

**Output**
- `poles::Vector{Complex}`: Extracted poles

**References**
[1] C. Rainieri and G. Fabbrocino. "Operational Modal Analysis of Civil Engineering Structures: An Introduction and Guide for Applications". Springer, 2014.


"""
function poles_extraction(prob::OMAProblem, order::Int, method::OMAModalExtraction = CovSSI(); stabdiag = false)

    return modes_extraction(prob, order, method, stabdiag = stabdiag)[1]
end

"""
    modes_extraction(prob, order, alg = CovSSI(); stabdiag)

Extract modal parameters from Covariance-based SSI method.

**Inputs**
- `prob::OMAProblem`: Structure containing half-spectrum data and frequency vector
- `order::Int`: Order of the system to identify
- `alg::CovSSI`: OMA method to use for modal extraction
- `stabdiag::Bool`: Whether to compute stabilization diagram (default: false)

**Outputs**
- `poles::Vector{Complex}`: Extracted poles
- `ms::Array{Complex, 2}`: Extracted mode shapes
"""
function modes_extraction(prob::OMAProblem, order::Int, alg::CovSSI; stabdiag = false)

    # Unpack problem data
    (; y, t, freq) = prob
    dt = t[2] - t[1]

    # Initialization
    nm = size(y, 1)
    nbr = ceil(Int, order/nm)    # number of block rows
    nbr < 10 ? nbr = 10 : nothing

    # Compute correlation functions
    R = xcorr(y, 2nbr)

    # Block Toeplitz matrix
    T = block_toeplitz(R, nbr)

    # Singular Value Decomposition of the block Toeplitz matrix
    F_svd = svd(T)
    U = F_svd.U[:, 1:order]
    S = Diagonal(F_svd.S[1:order])

    # Observability matrix
    Obs = U*sqrt.(S)

    # Output matrix
    C = Obs[1:nm, :]

    # System matrix
    Obs_up = Obs[1:end-nm, :]
    Obs_low = Obs[nm+1:end, :]
    A = Obs_up\Obs_low

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
    solve_modes(prob, order, alg = DataSSI(); stabdiag)

Extract modal parameters from Data-based SSI method.

**Inputs**
- `prob::OMAProblem`: Structure containing half-spectrum data and frequency vector
- `order::Int`: Order of the system to identify
- `alg::DataSSI`: OMA method to use for modal extraction
- `stabdiag::Bool`: Whether to compute stabilization diagram (default: false)

**Outputs**
- `poles::Vector{Complex}`: Extracted poles
- `ms::Array{Complex, 2}`: Extracted mode shapes
"""
function modes_extraction(prob::OMAProblem, order::Int, alg::DataSSI; stabdiag = false)

    # Unpack problem data
    (; y, t, freq) = prob
    dt = t[2] - t[1]

    # Initialization
    nm, nt = size(y)
    nbr = 2order    # number of block rows
    nbr < 10 ? nbr = 10 : nothing
    nt < 2nbr - 1 && throw(ArgumentError("Order too high for the number of samples."))

    # Block Hankel matrices
    H = block_hankel(y, nbr)

    # Projection matrix from LQ decomposition
    F_lq = lq(H)
    L = F_lq.L[nm*nbr+1:end, 1:nm*nbr]
    Q = Matrix(F_lq.Q)[1:nm*nbr, :]
    Proj = L*Q

    # Singular Value Decomposition of the projection matrix
    F_svd = svd(Proj)
    U = F_svd.U[:, 1:order]
    S = Diagonal(F_svd.S[1:order])

    # Observability matrix
    Obs = U*sqrt.(S)

    # Output matrix
    C = Obs[1:nm, :]

    # System matrix
    Obs_up = Obs[1:end-nm, :]
    Obs_low = Obs[nm+1:end, :]
    A = Obs_up\Obs_low

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