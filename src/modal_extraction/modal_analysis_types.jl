abstract type MdofProblem end

# Experimental Modal Analysis (EMA) - Sdof
abstract type SdofModalExtraction end
struct PeakPicking <: SdofModalExtraction end
struct CircleFit <: SdofModalExtraction end
struct LSFit <: SdofModalExtraction end

# Experimental Modal Analysis (EMA) - Mdof
abstract type MdofModalExtraction end
struct LSCE <: MdofModalExtraction end
struct LSCF <: MdofModalExtraction end
struct PLSCF <: MdofModalExtraction end

# Operational Modal Analysis (OMA)
abstract type OMAModalExtraction end
struct CovSSI <: OMAModalExtraction end
struct DataSSI <: OMAModalExtraction end

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
    freq::AbstractArray{R}
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
"""
    AutoEMASdofProblem(prob, alg; dpi, idx_m, idx_e)

Structure containing the input data for automatic experimental modal analysis using Sdof methods

**Fields**
* `prob::EMASdofProblem`: EMA-SDOF problem containing FRF data and frequency vector
* `alg::SdofModalExtraction`: Method to extract the poles
    * `PeakPicking`: Peak picking method (default)
    * `CircleFit`: Circle fitting method
    * `LSFit`: Least squares fitting method
* `dpi::Vector{Int}`: Driving point indices - default = [1, 1]
    * `dpi[1]`: Driving point index on the measurement mesh
    * `dpi[2]`: Driving point index on the excitation mesh
* `idx_m::AbstractArray{Int}`: Indices of measurement DOFs used for residues computation (default: all measurement DOFs)
* `idx_e::AbstractArray{Int}`: Indices of excitation DOFs used for residues computation (default: all excitation DOFs)
"""
@show_data struct AutoEMASdofProblem
    prob::EMAProblem
    alg::SdofModalExtraction
    dpi:: Vector{Int}
    idx_m::AbstractArray{Int}
    idx_e::AbstractArray{Int}

    AutoEMASdofProblem(prob::EMAProblem, alg::SdofModalExtraction = PeakPicking(); dpi::Vector{Int} = [1, 1], idx_m::AbstractArray{Int} = 1:size(prob.frf, 1), idx_e::AbstractArray{Int} = 1:size(prob.frf, 2)) = new(prob, alg, dpi, idx_m, idx_e)
end

## Structures for EMA-MDOF modal extraction
"""
    AutoEMAMdofProblem(prob, dpi, method; modetype)

Structure containing the input data for automatic experimental modal analysis using Mdof methods

**Fields**
* `prob::EMAProblem`: EMA problem containing FRF data and frequency vector
* `order::Int`: Model order (number of poles to extract)
* `dpi::Vector{Int}`: Driving point indices - default = [1, 1]
    * `dpi[1]`: Driving point index on the measurement mesh
    * `dpi[2]`: Driving point index on the excitation mesh
* `alg::MdofModalExtraction`: Method to extract the poles
    * `LSCE`: Least Squares Complex Exponential method
    * `LSCF``: Least Squares Complex Frequency method (default)
    * `PLSCF`: Polyreference Least Squares Complex Frequency method
* `modetype::Symbol`: Type of mode shapes to extract
    * `:real`: Real mode shapes (default)
    * `:complex`: Complex mode shapes
"""
@show_data struct AutoEMAMdofProblem
    prob::EMAProblem
    order::Int
    dpi:: Vector{Int}
    alg::MdofModalExtraction

    AutoEMAMdofProblem(prob::EMAProblem, order::Int, dpi::Vector{Int} = [1, 1], alg::MdofModalExtraction = LSCF()) = new(prob, order, dpi, alg)
end

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
@show_data struct StabilizationAnalysis{Tp <: Complex, Tf <: Real}
    prob::MdofProblem           # EMA-MDOF problem containing FRF data and frequency vector
    poles::Vector{Vector{Tp}}   # Extracted poles at each model order
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
- `win::Function`: Window function for PSD estimation (default: hanning)
- `overlap::Real`: Overlap ratio for PSD estimation (default: 0.5)

**Fields**
- `y::AbstractMatrix{Real}`: Matrix of measured outputs (channels x samples)
- `t::AbstractArray{Real}`: Time vector corresponding to the measurements
- `halfspec::Array{Complex, 3}`: Half power spectral density matrix
- `freq::AbstractArray{Real}`: Frequency vector corresponding to the half-spectrum
- `type_data::Symbol`: Type of measured data (:dis, :vel, :acc
"""
@show_data struct OMAProblem{C <: Complex, R <: Real} <: MdofProblem
    y::Matrix{R}
    yref::Matrix{R}
    t::AbstractArray{R}
    fullspec::Array{C, 3}
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

        return new{C, R}(Matrix{typeof(freq)(undef, 0, 0)}, Matrix{typeof(freq)(undef, 0, 0)}, similar(freq, 0), Gyy_red, halfspec, freq_red, type_data)
    end

    function OMAProblem(y::Matrix{R}, yref::Matrix{R}, t::AbstractArray{R}, fs, bs; frange = [0., fs/2.56], type_data = :dis, win = hanning, overlap = 0.5) where {R <: Real}

        # Compute half power spectral density matrix
        Gyy, freq = csd(y, yref, bs, win, fs = fs, overlap = overlap)

        # Correct frange to avoid division by zero
        frange[1] == 0. ? frange[1] = 1. : nothing

        # FRF post-processing - Frequency range reduction
        fidx = @. frange[1] ≤ freq ≤ frange[2]
        freq_red = freq[fidx]

        halfspec = half_psd(y, yref, freq_red, fs, bs)

        return new{Complex{R}, R}(y, yref, t, Gyy[:, :, fidx], halfspec, freq_red, type_data)
    end

    OMAProblem(y::Matrix{R}, t::AbstractArray{R}, fs, bs; frange = [0., fs/2.56], type_data = :dis, win = hanning, overlap = 0.5) where {R <: Real} = OMAProblem(y, y, t, fs, bs; frange = frange, type_data = type_data, win = win, overlap = overlap)
end