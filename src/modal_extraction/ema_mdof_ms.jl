"""
    residues(frf, freq, poles; frange = [freq[1], freq[end]], type = :dis)

Compute residues and lower and upper residuals of a frequency response function (FRF) given its poles.

**Inputs**
- `frf`: Frequency response function (vector)
- `freq`: Frequency vector (Hz)
- `poles`: Vector of complex poles
**Keyword Arguments**
- `frange`: Frequency range for residue calculation (default: full range)
- `type`: Type of FRF
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance

**Outputs**
- `res`: Residues corresponding to each pole
- `lr`: Lower residual
- `ur`: Upper residual
"""
function residues(frf, freq, poles; frange = [freq[1], freq[end]], type = :dis)

    # Correct frange to avoid division by zero
    frange[1] = frange[1] < 1. ? 1. : frange[1]

    # FRF post-processing - Frequency range reduction
    fidx = @. frange[1] ≤ freq ≤ frange[2]
    frf_red = permutedims(frf[:, :, fidx], (3, 1, 2)) # Put the last dimension first for subsequent calculations
    freq_red = freq[fidx]
    nm, ne = size(frf_red)[2:3]

    # Convert frf to admittance if needed
    om = 2π*freq_red
    if type != :dis
        for (f, omega) in enumerate(om)
            if type == :vel
                frf_red[f, :, :] ./= 1im*omega
            elseif type == :acc
                frf_red[f, :, :] ./= -omega^2
            end
        end
    end

    # Residue calculation
    valid_poles = !isnan.(poles)
    p = [poles[valid_poles]; conj.(poles[valid_poles])]
    np = length(p)

    res = similar(frf, np + 2, nm, ne)
    P = [(1 ./(1im*om .- transpose(p))) (-1 ./om.^2) ones(length(om))]
    for i in 1:ne
        res[:, :, i] .= P\frf_red[:, :, i]
    end

    # Return residues, lower and upper residuals
    return res[1:Int(np/2), : , :], res[end - 1, : , :], res[end, : , :]
end

"""
    modeshape_extraction(residues, poles, di = [1, 1]; type = :complex)

Extract mode shapes using MDOF approximation

**Inputs**
- `residues`: Residues matrix of size (np, nm, ne)
- `poles`: Vector of complex poles
* `dpi`: Driving point indices - default = [1, 1]
    * `dpi[1]`: Driving point index on the measurement mesh
    * `dpi[2]`: Driving point index on the excitation mesh
- `type`: Type of mode shape
    * `:complex`: Complex mode shapes (default)
    * `:real`: Real mode shapes

**Outputs**
- `ϕ`: Mode shapes matrix
- `Q`: Scaling factors vector
"""
function modeshape_extraction(residues, poles::Vector{T}, dpi = [1, 1]; type = :complex) where {T <: Complex}


end