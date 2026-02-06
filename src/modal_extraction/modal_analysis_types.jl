abstract type MdofProblem end

# Experimental Modal Analysis (EMA) - Sdof
abstract type SdofEMA end

"""
    PeakPicking

Peak picking method for Sdof modal extraction

**Fields**
* `nothing`

**References**

[1] D.J. Ewins, "Modal Testing: Theory, Practice and Application", 2nd Edition, Research Studies Press, 2000.

[2] D. J. Inman, "Engineering Vibration", 4th Edition, Pearson, 2013.
"""
struct PeakPicking <: SdofEMA end

"""
    PeakPicking

Circle fitting method for Sdof modal extraction

**Fields**
* `nothing`

**References**

[1] D.J. Ewins, "Modal Testing: Theory, Practice and Application", 2nd Edition, Research Studies Press, 2000.

[2] D. J. Inman, "Engineering Vibration", 4th Edition, Pearson, 2013.
"""
struct CircleFit <: SdofEMA end

"""
    LSFit

Least squares fitting method for Sdof modal extraction

**Fields**
* `nothing`

**Reference**

[1] A. Brandt, "Noise and Vibration Analysis: Signal Analysis and Experimental Procedures", Wiley, 2011.
"""
struct LSFit <: SdofEMA end

# Experimental Modal Analysis (EMA) - Mdof
abstract type MdofEMA end

"""
    LSCE

Least Squares Complex Exponential method for Mdof modal extraction

**Fields**
* `nothing`

**Reference**

[1] D. L. Brown, R. J. Allemang, Ray Zimmerman and M. Mergeay. "Parameter Estimation Techniques for Modal Analysis". SAE Transactions, vol. 88, pp. 828-846, 1979.
"""
struct LSCE <: MdofEMA end

"""
    LSCF

Least Squares Complex Frequency method for Mdof modal extraction

**Fields**
* `nothing`

**References**

[1] P. Guillaume, P. Verboven, S. Vanlanduit, H. Van der Auweraer, and B. Peeters, "A poly-reference implementation of the least-squares complex frequency-domain estimator". In Proceedings of IMAC XXI. Kissimmee, FL, 2003.

[2] El-Kafafy M., Guillaume P., Peeters B., Marra F., Coppotelli G. (2012).Advanced Frequency-Domain Modal Analysis for Dealing with Measurement Noise and Parameter Uncertainty. In: Allemang R., De Clerck J., Niezrecki C., Blough J. (eds) Topics in Modal Analysis I, Volume 5. Conference Proceedings of the Society for Experimental Mechanics Series. Springer, New York, NY
"""
struct LSCF <: MdofEMA end

"""
    pLSCF

Polyreference Least Squares Complex Frequency method for Mdof modal extraction

**Fields**
* `nothing`

**References**

[1] P. Guillaume, P. Verboven, S. Vanlanduit, H. Van der Auweraer, and B. Peeters, "A poly-reference implementation of the least-squares complex frequency-domain estimator". In Proceedings of IMAC XXI. Kissimmee, FL, 2003.

[2] El-Kafafy M., Guillaume P., Peeters B., Marra F., Coppotelli G. (2012).Advanced Frequency-Domain Modal Analysis for Dealing with Measurement Noise and Parameter Uncertainty. In: Allemang R., De Clerck J., Niezrecki C., Blough J. (eds) Topics in Modal Analysis I, Volume 5. Conference Proceedings of the Society for Experimental Mechanics Series. Springer, New York, NY
"""
struct pLSCF <: MdofEMA end

# Operational Modal Analysis (OMA) - Sdof
abstract type SdofOMA end

"""
    FSDD

Frequency-Spatial Domain Decomposition method for OMA modal extraction

**Fields**
* `nothing`

**Reference**

[1] L. Zhang, T. Wang and Y. Tamura. "A frequency–spatial domain decomposition (FSDD) method for operational modal analysis". Mechanical Systems and Signal Processing, 24: 1227-1239, 2010.
"""
struct FSDD <: SdofOMA end

# Operational Modal Analysis (OMA) - Mdof
abstract type MdofOMA end

"""
    CovSSI

Covariance-driven Stochastic Subspace Identification method for OMA modal extraction

**Fields**
* `nothing`

**References**

[1] C. Rainieri and G. Fabbrocino. "Operational Modal Analysis of Civil Engineering Structures: An Introduction and Guide for Applications". Springer, 2014.

[2] P. Peeters and G. De Roeck. "Reference-based stochastic subspace identification for output-only modal analysis". Mechanical Systems and Signal Processing, 13(6):855-878, 1999.

[3] L. Hermans and H. Van der Auweraer. "Modal testing and analysis of structures under operational conditions: Industrial applications". Mechanical Systems and Signal Processing, 13(2):193-216, 1999.
"""
struct CovSSI <: MdofOMA end

"""
    DataSSI

Data-driven Stochastic Subspace Identification method for OMA modal extraction

**Fields**
* `nothing`

**References**

[1] C. Rainieri and G. Fabbrocino. "Operational Modal Analysis of Civil Engineering Structures: An Introduction and Guide for Applications". Springer, 2014.

[2] P. Peeters and G. De Roeck. "Reference-based stochastic subspace identification for output-only modal analysis". Mechanical Systems and Signal Processing, 13(6):855-878, 1999.

[3] L. Hermans and H. Van der Auweraer. "Modal testing and analysis of structures under operational conditions: Industrial applications". Mechanical Systems and Signal Processing, 13(2):193-216, 1999.
"""
struct DataSSI <: MdofOMA end

## Common structures for EMA
"""
    EMAProblem(frf, freq; frange, type_frf)

Data structure defining the inputs for EMA modal extraction methods.

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
@show_data struct EMAProblem{C <: Complex, R <: Real} <: MdofProblem
    frf::Array{C, 3}
    freq::Vector{R}
    type_frf::Symbol

    function EMAProblem(frf::Array{C, 3}, freq::AbstractArray{R}; frange = [freq[1], freq[end]], type_frf = :dis) where {C <: Complex, R <: Real}

        # Correct frange to avoid division by zero
        frange[1] == 0. ? frange[1] = 1. : nothing

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

        return new{C, R}(frf_red, freq[fidx], type_frf)
    end
end

"""
    EMASolution(poles, ms, ci, res, lr, ur)

Structure containing the solution of the automatic experimental modal analysis using Mdof methods

**Fields**
* `poles::Vector{Complex}`: Extracted poles
* `ms::AbstractArray{Real}`: Mode shapes
* `ci::Vector{Complex}`: Scaling constants for each mode
* `res::AbstractArray{Complex}`: Residues for each mode
* `lr::Matrix{Complex}`: Lower residual
* `ur::Matrix{Complex}`: Upper residual
"""
@show_data struct EMASolution{Tp <: Complex, Tm <: Real}
    poles::Vector{Tp}
    ms::AbstractArray{Tm}
    ci::Vector{Tp}
    res::AbstractArray{Tp}
    lr::Matrix{Tp}
    ur::Matrix{Tp}
end

## Structures for EMA-SDOF modal extraction
# """
#     AutoEMASdofProblem(prob, alg; dpi, idx_m, idx_e, width, min_prom, max_prom, pks_indices)

# Structure containing the input data for automatic experimental modal analysis using Sdof methods

# **Fields**
# * `prob::EMASdofProblem`: EMA-SDOF problem containing FRF data and frequency vector
# * `alg::SdofModalExtraction`: Method to extract the poles
#     * `PeakPicking`: Peak picking method (default)
#     * `CircleFit`: Circle fitting method
#     * `LSFit`: Least squares fitting method
# * `dpi::Vector{Int}`: Driving point indices - default = [1, 1]
#     * `dpi[1]`: Driving point index on the measurement mesh
#     * `dpi[2]`: Driving point index on the excitation mesh
# * `idx_m::AbstractArray{Int}`: Indices of measurement DOFs used for residues computation (default: all measurement DOFs)
# * `idx_e::AbstractArray{Int}`: Indices of excitation DOFs used for residues computation (default: all excitation DOFs)
# * `width::Int`: Half-width of the peaks (default: 1)
# * `min_prom::Real`: Minimum peak prominence (default: 0.)
# * `max_prom::Real`: Maximum peak prominence (default: Inf)
# * `pks_indices::Vector{Int}`: Predefined peak indices (default: empty vector)
# """
# @show_data struct AutoEMASdofProblem{R <: Real}
#     prob::EMAProblem
#     alg::SdofEMA
#     dpi:: Vector{Int}
#     idx_m::AbstractArray{Int}
#     idx_e::AbstractArray{Int}
#     width::Int
#     min_prom::R
#     max_prom::R
#     pks_indices::Vector{Int}

#     AutoEMASdofProblem(prob::EMAProblem, alg::SdofEMA = PeakPicking(); dpi::Vector{Int} = [1, 1], idx_m::AbstractArray{Int} = 1:size(prob.frf, 1), idx_e::AbstractArray{Int} = 1:size(prob.frf, 2), width::Int = 1, min_prom::R = 0., max_prom::R = Inf, pks_indices::Vector{Int} = Int[]) where {R <: Real} = new{R}(prob, alg, dpi, idx_m, idx_e, width, min_prom, max_prom, pks_indices)
# end

## Structures for EMA-MDOF modal extraction
# """
#     AutoEMAMdofProblem(prob, dpi, method; modetype)

# Structure containing the input data for automatic experimental modal analysis using Mdof methods

# **Fields**
# * `prob::EMAProblem`: EMA problem containing FRF data and frequency vector
# * `order::Int`: Model order (number of poles to extract)
# * `dpi::Vector{Int}`: Driving point indices - default = [1, 1]
#     * `dpi[1]`: Driving point index on the measurement mesh
#     * `dpi[2]`: Driving point index on the excitation mesh
# * `alg::MdofEMA`: Method to extract the poles
#     * `LSCE`: Least Squares Complex Exponential method
#     * `LSCF``: Least Squares Complex Frequency method (default)
#     * `PLSCF`: Polyreference Least Squares Complex Frequency method
# * `modetype::Symbol`: Type of mode shapes to extract
#     * `:real`: Real mode shapes (default)
#     * `:complex`: Complex mode shapes
# """
# @show_data struct AutoEMAMdofProblem
#     prob::EMAProblem
#     order::Int
#     alg::MdofEMA
#     dpi:: Vector{Int}

#     AutoEMAMdofProblem(prob::EMAProblem, order::Int, alg::MdofEMA = LSCF(); dpi::Vector{Int} = [1, 1]) = new(prob, order, alg, dpi)
# end

"""
    StabilizationAnalysis(prob, poles, modefn, mode_stabfn, mode_stabdr)

Data structure summarizing the results of the stabilization analysis.

**Fields**
- `prob::MdofProblem`: EMA-MDOF problem containing FRF data and frequency vector
- `frange::Vector{Real}`: Frequency range used for the stabilization analysis
- `poles::Vector{Vector{Complex}}`: Vector of vectors containing extracted poles at each model order
- `modefn::Matrix{Real}`: Matrix containing the natural frequencies (useful for plotting)
- `mode_stabfn::Matrix{Bool}`: Matrix indicating the stability of natural frequencies
- `mode_stabdr::Matrix{Bool}`: Matrix indicating the stability of damping ratios

**Note**

This structure is returned by the `stabilization` function after performing a stabilization diagram analysis and used by `stabilization_plot` for visualization.
"""
@show_data struct StabilizationAnalysis{Tc <: Complex, Tf <: Real}
    prob::MdofProblem           # EMA-MDOF problem containing FRF data and frequency vector
    poles::Vector{Vector{Tc}}   # Extracted poles at each model order
    modefn::Matrix{Tf}          # Natural frequencies (used for plotting)
    mode_stabfn::BitMatrix      # Stability of natural frequencies
    mode_stabdr::BitMatrix      # Stability of damping ratios
end

## Structures for OMA
"""
    OMAProblem(Gyy, freq; frange, type_data)
    OMAProblem(y, yref, t, fs, bs; frange, type_data, win, overlap)
    OMAProblem(y, t, fs, bs; frange, type_data, win, overlap)

Data structure defining the inputs for Operational Modal Analysis (OMA) methods.

**Constructor parameters**
- `Gyy::Array{Complex, 3}`: Cross-spectral density matrix of size (no, nref, nf),
- `freq::AbstractArray{Real}`: Frequency vector corresponding to the CSD matrix
- `frange::AbstractVector{Real}`: Frequency range to consider for analysis (default: full range)

**Alternative constructor parameters**
- `y::AbstractMatrix{Real}`: Matrix of measured outputs (no x nt)
- `yref::AbstractMatrix{Real}`: Matrix of reference outputs (nref x nt)
- `t::AbstractArray{Real}`: Time vector corresponding to the measurements
- `fs::Int`: Sampling frequency (Hz)
- `bs::Int`: Block size for CSD estimation
- `frange::AbstractVector{Real}`: Frequency range to consider for analysis (default: [0., fs/2.56])
- `win::Function`: Window function for CSD estimation (default: hanning)
- `overlap::Real`: Overlap ratio for CSD estimation (default: 0.5)

**Fields**
- `y::AbstractMatrix{Real}`: Matrix of measured outputs (no x nt)
- `yref::AbstractMatrix{Real}`: Matrix of reference outputs (nref x nt)
- `t::AbstractArray{Real}`: Time vector corresponding to the measurements
- `fullspec::Array{Complex, 3}`: Full cross-spectral density matrix
- `halfspec::Array{Complex, 3}`: Half spectral density matrix
- `freq::AbstractArray{Real}`: Frequency vector corresponding to the half-spectrum

**Note**
- `OMAProblem(y, t, ...) = OMAProblem(y, y, t, ...)`
"""
@show_data struct OMAProblem{C <: Complex, R <: Real} <: MdofProblem
    y::Matrix{R}
    yref::Matrix{R}
    t::Vector{R}
    fullspec::Array{C, 3}
    halfspec::Array{C, 3}
    freq::Vector{R}

    function OMAProblem(Gyy::Array{C, 3}, freq::AbstractArray{R}; frange = [freq[1], freq[end]]) where {C <: Complex, R <: Real}

        # Correct frange to avoid division by zero
        frange[1] == 0. ? frange[1] = 1. : nothing

        # FRF post-processing - Frequency range reduction
        fidx = @. frange[1] ≤ freq ≤ frange[2]
        freq_red = freq[fidx]

        ω = 2π*freq_red
        Gyy_red = Gyy[:, :, fidx]

        # Compute half spectrum
        halfspec = half_csd(Gyy_red, freq_red)

        return new{C, R}(Matrix{eltype(freq)}(undef, 0, 0), Matrix{eltype(freq)}(undef, 0, 0), similar(freq_red, 0), Gyy_red, halfspec, freq_red)
    end

    function OMAProblem(y::Matrix{R}, yref::Matrix{R}, t::AbstractArray{R}, fs, bs; frange = [0., fs/2.56], win = hanning, overlap = 0.5) where {R <: Real}

        # Compute half power spectral density matrix
        Gyy, freq = csd(yref, y, bs, win, fs = fs, overlap = overlap)

        # Correct frange to avoid division by zero
        frange[1] == 0. ? frange[1] = 1. : nothing

        # FRF post-processing - Frequency range reduction
        fidx = @. frange[1] ≤ freq ≤ frange[2]
        freq_red = freq[fidx]
        Gyy_red = Gyy[:, :, fidx]

        # Compute half spectrum
        halfspec = half_csd(y, yref, freq_red, fs, bs)

        return new{Complex{R}, R}(y, yref, collect(t), Gyy_red, halfspec, freq_red)
    end

    OMAProblem(y::Matrix{R}, t::AbstractArray{R}, fs, bs; frange = [0., fs/2.56], win = hanning, overlap = 0.5) where {R <: Real} = OMAProblem(y, y, t, fs, bs; frange = frange, win = win, overlap = overlap)
end