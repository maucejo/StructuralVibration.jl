"""
    frf_reconstruction(res, poles, freq; lr, ur, type)

Reconstruct a frequency response function (FRF) from its residues and poles.

**Inputs**
- `res::Array{Complex, 3}`: Residues corresponding to each pole
- `poles::Vector{Complex}`: Poles extracted from the FRF
- `freq::Vector{Float64}`: Frequency vector
- `lr::Matrix{Complex, 2}`: Lower residuals (default: zeros)
- `ur::Matrix{Complex, 2}`: Upper residuals (default: zeros)
- `type::Symbol`: Type of FRF to reconstruct
    - `:dis`: displacement (default)
    - `:vel`: velocity
    - `:acc`: acceleration

**Output**
- `H_rec::Array{Complex, 3}`: Reconstructed FRF
"""
function frf_reconstruction(res::Array{T, 3}, poles::Vector{T}, freq; lr = zeros(eltype(res), size(res)[2:end]), ur = zeros(eltype(res), size(res)[2:end]), type = :dis) where {T <: Complex}

    # Initialization
    ω = 2π*freq
    Res = [res; conj.(res)]
    p = [poles; conj.(poles)]

    nf = length(freq)
    nm, ne = size(res)[2:end]
    H_rec = similar(res, nm, ne, nf)
    for (f, ωf) in enumerate(ω)
        for i in 1:nm
            for j in 1:ne
                H_rec[i, j, f] = sum(Res[:, i, j] ./ (im*ωf .- p)) - lr[i, j]/ωf^2 + ur[i, j]

                if type == :vel
                    H_rec[i, j, f] *= 1im*ωf
                elseif type == :acc
                    H_rec[i, j, f] *= -ωf^2
                end
            end
        end
    end

    return H_rec
end

"""
    compute_residuals(prob, res, poles)

Compute lower and upper residuals of a frequency response function (FRF) given its residues and poles.

**Inputs**
- `prob::Union{EMASdofProblem, EMAMdofProblem}`: Structure containing FRF data and frequency vector
- `res::Array{Complex, 3}`: Residues corresponding to each pole
- `poles::Vector{T}`: Poles extracted from the FRF

**Outputs**
- `lr::Matrix{Complex, 2}`: Lower residuals
- `ur::Matrix{Complex, 2}`: Upper residuals
"""
function compute_residuals(prob:: Union{EMASdofProblem, EMAMdofProblem}, res::Array{T, 3}, poles::Vector{T}) where {T <: Complex}

    # Initialization
    (; frf, freq) = prob

    # Correct frange to avoid division by zero
    if freq[1] < 1.
        freq_red = freq[2:end]
        frf_red = frf[:, :, 2:end]
    else
        freq_red = freq
        frf_red = frf
    end

    nm, ne, nf = size(frf_red)
    ω = 2π*freq_red

    # FRF reconstruction
    H_rec = frf_reconstruction(res, poles, freq_red)

    # Residual computation
    Y = permutedims(frf_red .- H_rec, (3, 1, 2))

    R = [(-1. ./ ω.^2) ones(nf)]
    Comp = similar(res, 2, nm, ne)
    for i in 1:ne
        Comp[:, :, i] = R\Y[:, :, i]
    end

    return Comp[1, :, :], Comp[2, :, :]
end

"""
    mode2residues(ms, poles, idx_m, idx_e)

Compute the mode residues from mode shapes and poles.

**Inputs**
- `ms::Matrix{Complex, 2}`: Mode shapes matrix
- `poles::Vector{Complex}`: Poles extracted from the FRF
- `idx_m::Vector{Int}`: Indices of measurement locations
- `idx_e::Vector{Int}`: Indices of excitation locations

**Output**
- `res::Array{Complex, 3}`: Residues corresponding to each pole

**Note**
The mode shapes are supposed to be real and mass-normalized.
"""
function mode2residues(ms, poles, idx_m, idx_e)
    nm = length(idx_m)
    ne = length(idx_e)
    n_modes = length(poles)

    fn, ξn = poles2modal(poles)
    Ωn = @. 2π*fn*√(1 - ξn^2)

    res = zeros(eltype(poles), n_modes, nm, ne)
    for (i, mes) in enumerate(idx_m)
        for (j, exc) in enumerate(idx_e)
            res[:, i, j] .= ms[mes, :].*ms[exc, :]./(2im*Ωn)
        end
    end

    return res
end