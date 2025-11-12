abstract type OMAModalExtraction end
struct CovSSI <: OMAModalExtraction end
struct DataSSI <: OMAModalExtraction end

"""
    OMAMdofProblem(Gyy, freq; frange, type_frf)

Data structure defining the inputs for EMA-MDOF modal extraction methods.

**Constructor parameters**
- `Gyy::Array{Complex, 3}`: 3D power spectral density matrix (array nm x ne x nf)
- `freq::AbstractArray{Real}`: Vector of frequency values (Hz)
- `frange::Vector{Real}`: Frequency range for analysis (default: [freq[1], freq[end]])
- `type_spec::Symbol`: Type of spectral density used in the analysis
    * `:dis`: Displacement (default)
    * `:vel`: Velocity
    * `:acc`: Acceleration

**Fields**
- `halfspec::Array{Complex, 3}`: 3D admittance matrix (array nm x ne x nf)
- `freq::AbstractArray{Real}`: Vector of frequency values (Hz)
- `type_spec::Symbol`: Type of spectral density used in the analysis

**Note**
"""
@show_data struct OMAMdofProblem{C <: Complex, R <: Real} <: MdofProblem
    halfspec::Array{C, 3}
    freq::AbstractArray{R}
    type_spec::Symbol

    function OMAMdofProblem(Gyy::Array{C, 3}, freq::AbstractArray{R}; frange = [freq[1], freq[end]], type_spec = :dis) where {C <: Complex, R <: Real}

        # Correct frange to avoid division by zero
        frange[1] == 0. ? frange[1] = 1. : nothing

        # FRF post-processing - Frequency range reduction
        fidx = @. frange[1] ≤ freq ≤ frange[2]
        ω = 2π*freq[fidx]
        Gyy_red = Gyy[:, :, fidx]

        # Conversion to admittance
        for (f, ωf) in enumerate(ω)
            if type_spec == :vel
                Gyy_red[:, :, f] ./= -ωf^2
            elseif type_spec == :acc
                Gyy_red[:, :, f] ./= ωf^4
            end
        end

        return new{C, R}(Gyy_red, freq[fidx], type_spec)
    end
end

"""
    OMAProblem(y, yref = y, t)

Data structure defining the inputs for Operational Modal Analysis (OMA) methods.

**Fields**
- `y::AbstractMatrix`: Measured output response - dimensions nm x nt
- `yref::AbstractMatrix`: Reference output responses - dimensions nr x nt (default: y)
- `t::AbstractArray`: Time vector (nt)
"""
@show_data struct OMAProblem{C <: Complex, R <: Real}
    y::AbstratctMatrix{R}
    yref::AbstractMatrix{R}
    t::AbstractArray{R}
    halfspec::Array{C, 3}
    freq::AbstractArray{R}

    function OMAProblem(y::AbstractMatrix{R}, yref::AbstractMatrix{R} = y, t::AbstractArray{R}, freq::AbstractArray{R}, fs, bs) where {R <: Real}

        # Compute half power spectral density matrix
        halspec = half_psd(y, yref, freq, fs, bs)

        return new{C, R}(y, yref, t, halspec, freq)
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
"""
function poles_extraction(prob::OMAProblem, order::Int, method::OMAModalExtraction = CovSSI(); stabdiag = false)

    return solve_modes(prob, order, method, stabdiag = stabdiag)[1]
end

"""
    solve_modes(prob, order, alg = CovSSI(); stabdiag)

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
function solve_modes(prob::OMAProblem, order::Int, alg::CovSSI; stabdiag = false)

    # Unpack problem data
    (; y, yref) = prob

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
function solve_modes(prob::OMAProblem, order::Int, alg::DataSSI; stabdiag = false)

end

"""
    build_hankel(y, nbr)

Constructs the past and future block Hankel matrices for SSI-DATA.

**Inputs**
- `y::Matrix{Float64}`: Matrix of measured outputs (channels × samples)
- `nbr::Int`: number of block rows

**Outputs**
- `Hp::Matrix{Float64}`: Past block Hankel matrix
- `Hf::Matrix{Float64}`: Future block Hankel matrix

**References**
[1] C. Ranieri and G. Fabbrocino. "Operational Modal Analysis of Civil Engineering Structures: An Introduction and Guide for Applications". Springer, 2014.
"""
function build_hankel(y::Matrix{Float64}, nbr::Int)
    nm, nt = size(y)
    nc = nt - 2nbr + 1 # number of columns of Hankel matrices

    Yp = similar(y, nm*n, nc) # Past Hankel matrix
    Yf = similar(Hp)          # Future Hankel matrix
    for row = 1:nbr
        Yp[(row-1)*nm+1:row*nm, :] .= y[:, row:row+nc-1]
        Yf[(row-1)*nm+1:row*nm, :] .= y[:, row+nbr:row+nbr+nc-1]
    end

    return Yp, Yf
end