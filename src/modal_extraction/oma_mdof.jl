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

[2] P. Peeters and G. De Roeck. "Reference-based stochastic subspace identification for output-only modal analysis". Mechanical Systems and Signal Processing, 13(6):855-878, 1999.

[3] L. Hermans and H. Van der Auweraer. "Modal testing and analysis of structures under operational conditions: Industrial applications". Mechanical Systems and Signal Processing, 13(2):193-216, 1999.
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
    (; y, yref, t, freq) = prob

    # Initialization
    dt = t[2] - t[1]
    no = size(y, 1)
    ni = size(yref, 1)
    nbr = ceil(Int, order/min(no, ni))    # number of block rows
    nbr < 10 ? nbr = 10 : nothing

    # Compute correlation functions
    R = xcorr(y, yref, 2nbr)

    # Block Toeplitz matrix
    T = block_toeplitz(R, nbr)

    # Compute system matrices A and C from block Toeplitz matrix
    A, C = compute_ss_matrix(T, no, order)

    return oma_modal_parameters(A, C, dt, order, freq, stabdiag)
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
@views function modes_extraction(prob::OMAProblem, order::Int, alg::DataSSI; stabdiag = false)

    # Unpack problem data
    (; y, yref, t, freq) = prob
    dt = t[2] - t[1]

    # Initialization
    no, nt = size(y)
    nbr = 2order    # number of block rows
    nbr < 10 ? nbr = 10 : nothing
    nt < 2nbr - 1 && throw(ArgumentError("Order too high for the number of samples."))

    # Block Hankel matrices
    Yp, Yf = block_hankel(y, yref, nbr)

    # Projection matrix from LQ decomposition)
    # H = [Yp; Yf]
    # F_lq = lq(H)
    # L = F_lq.L[no*nbr+1:end, 1:no*nbr]
    # Q = Matrix(F_lq.Q)[1:no*nbr, :]
    # Proj = L*Q

    # SVD of Yp for stable pseudo-inverse
    svd_tol = 1e-12
    Up, Sp, Vp = svd(Yp*Yp')
    tol = svd_tol * Sp[1]
    r = count(s -> s > tol, Sp)
    Sp_inv = Diagonal(1.0 ./ Sp[1:r])
    Yp_pinv = Vp[:, 1:r] * Sp_inv * Up[:, 1:r]'

    # Projection of future onto past
    Proj = Yf * Yp'Yp_pinv * Yp

    # Compute system matrices A and C from projection matrix
    A, C = compute_ss_matrix(Proj, no, order)

    return oma_modal_parameters(A, C, dt, order, freq, stabdiag)
end