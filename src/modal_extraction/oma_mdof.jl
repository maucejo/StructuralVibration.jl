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

    # # Eigen decomposition of system matrix
    # eigvals, eigvecs = eigen(A)

    # # Extract poles and mode shapes
    # p = log.(eigvals)/dt

    # # Poles that do not appear as pairs of complex conjugate numbers with a positive real part or that are purely real are suppressed. For complex conjugate poles, only the pole with a positive imaginary part is retained.
    # valid_poles = p[@. imag(p) < 0. && real(p) ≤ 0. && !isreal(p)]
    # idx_valid_poles = findall(in(conj.(valid_poles)), p)
    # p_valid = p[idx_valid_poles]
    # idx_sort_poles = sortperm(p_valid, by = abs)
    # p_sort = p_valid[idx_sort_poles]

    # # Check if poles are in frange
    # fn = poles2modal(p_sort)[1]
    # fidx = @. freq[1] ≤ fn ≤ freq[end]
    # poles = p_sort[fidx]

    # # Extract mode shapes
    # ms = C*eigvecs[:, idx_valid_poles[idx_sort_poles[fidx]]]

    # return stabdiag ? [poles; fill(complex(NaN, NaN), order - length(poles))] : poles, ms

    return oma_modal_parameters(A, C, dt, freq, stabdiag)
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

    # # Eigen decomposition of system matrix
    # eigvals, eigvecs = eigen(A)

    # # Extract poles and mode shapes
    # p = log.(eigvals)/dt

    # # Poles that do not appear as pairs of complex conjugate numbers with a positive real part or that are purely real are suppressed. For complex conjugate poles, only the pole with a positive imaginary part is retained.
    # valid_poles = p[@. imag(p) < 0. && real(p) ≤ 0. && !isreal(p)]
    # idx_valid_poles = findall(in(conj.(valid_poles)), p)
    # p_valid = p[idx_valid_poles]
    # idx_sort_poles = sortperm(p_valid, by = abs)
    # p_sort = p_valid[idx_sort_poles]

    # # Check if poles are in frange
    # fn = poles2modal(p_sort)[1]
    # fidx = @. freq[1] ≤ fn ≤ freq[end]
    # poles = p_sort[fidx]

    # # Extract mode shapes
    # ms = C*eigvecs[:, idx_valid_poles[idx_sort_poles[fidx]]]

    # return stabdiag ? [poles; fill(complex(NaN, NaN), order - length(poles))] : poles, ms

    return oma_modal_parameters(A, C, dt, freq, stabdiag)
end