"""
    half_csd_reconstruction(res, poles, freq; lr, ur, type)

Reconstruct a half-spectrum from its residues and poles.

**Inputs**
- `res::Array{Complex, 3}`: Residues corresponding to each pole
- `poles::Vector{Complex}`: Poles extracted from the half-spectrum
- `freq::Vector{Float64}`: Frequency vector
- `lr::Matrix{Complex, 2}`: Lower residuals (default: zeros)
- `ur::Matrix{Complex, 2}`: Upper residuals (default: zeros)

**Output**
- `H_rec::Array{Complex, 3}`: Reconstructed FRF
"""
function half_csd_reconstruction(res::Array{T, 3}, poles::Vector{T}, freq; lr = zeros(eltype(res), size(res)[2:end]), ur = zeros(eltype(res), size(res)[2:end])) where {T <: Complex}

    # Initialization
    ω = 2π*freq
    Res = [res; conj.(res)]
    p = [poles; conj.(poles)]

    nf = length(freq)
    no, ni = size(res)[2:end]
    S_rec = similar(res, no, ni, nf)
    for (f, ωf) in enumerate(ω)
        for j in 1:ni
            for i in 1:no
                S_rec[i, j, f] = sum(Res[:, i, j] ./ (im*ωf .- p)) + lr[i, j]/(im*ωf) + 1im*ωf*ur[i, j]
            end
        end
    end

    return S_rec
end

"""
    compute_residuals(prob, res, poles)

Compute lower and upper residuals of a half-spectrum given its residues and poles.

**Inputs**
- `prob::OMAProblem`: Structure containing half-spectrum data and frequency vector
- `res::Array{Complex, 3}`: Residues corresponding to each pole
- `poles::Vector{T}`: Poles extracted from the half-spectrum

**Outputs**
- `lr::Matrix{Complex, 2}`: Lower residuals
- `ur::Matrix{Complex, 2}`: Upper residuals
"""
function compute_residuals(prob:: OMAProblem, res::Array{T, 3}, poles::Vector{T}) where {T <: Complex}

    # Initialization
    (; halfspec, freq) = prob

    # Correct frange to avoid division by zero
    if freq[1] < 1.
        freq_red = freq[2:end]
        Gyy_red = halfspec[:, :, 2:end]
    else
        freq_red = freq
        Gyy_red = halfspec
    end

    no, ni = size(Gyy_red)[1:2]
    ω = 2π*freq_red

    # FRF reconstruction
    Gyy_rec = half_csd_reconstruction(res, poles, freq_red)

    # Residual computation
    Y = permutedims(Gyy_red .- Gyy_rec, (3, 1, 2))

    R = [(1. ./ (1im*ω)) (1im*ω)]
    Comp = similar(res, 2, no, ni)
    for i in 1:ni
        Comp[:, :, i] = R\Y[:, :, i]
    end

    return Comp[1, :, :], Comp[2, :, :]
end