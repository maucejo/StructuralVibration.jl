"""
    mof(fn, dr)

Compute the modal overlap factor.

**Inputs**
- `fn`: Natural frequencies (Hz)
- `dr`: Damping ratios

**Output**
- `mof`: Modal overlap factor

**References**

[1] Srikantha Phani A, Woodhouse J (2007) Viscous damping identification in linear vibration. Journal of Sound and Vibration. 303: 475 – 500
"""
function mof(fn, dr)
    mof = similar(fn)

    for i in eachindex(fn)
        if i == 1
            mof[i] = fn[i]*dr[i]/(fn[i + 1] - fn[i])
        else
            mof[i] = fn[i]*dr[i]/(fn[i] - fn[i - 1])
        end
    end

    return length(mof) == 1 ? mof[1] : mof
end

## Mode complexity
"""
    mov(poles, ms, ci)

Compute the mode overcomplexity value

This indicator is a weighted percentage of the degrees of freedom of the response for which adding mass leads to a negative frequency shift. A value close to 1 indicates physical modes, while a low value indicates non-physical modes (numerical or noise-related)

**Inputs**
- `poles`: Poles of the system
- `ms`: Mode shapes (each column corresponds to a mode)
- `ci`: Scaling factors associated to the mode shapes

**Output**
- `mov`: Mode overcomplexity values

**Reference**

[1] M+P Analyzer manual. Rev. 5.1. 2017
"""
function mov(poles, ms, ci)
    if ms isa Vector
        ms = reshape(ms, :, 1)
    end

    n, m = size(ms)
    mov = zeros(m)
    s = zeros(n)
    Msen = similar(s)
    sig = similar(s)

    for (i, (p, phi, Qi)) in enumerate(zip(poles, eachcol(ms), ci))
        @. Msen = imag(-p^2*phi^2*Qi) # Poles sensitivity w.r.t. mass
        @. sig = sign(Msen)
        s[sig .> 0.] .= 0.
        s[sig .<= 0.] .= 1.

        mov[i] = sum(s.*abs2.(phi))/sum(abs2, phi)
    end

    return length(mov) == 1 ? mov[1] : mov
end

"""
    mpc(ms)

Compute the mode phase collinearity.

This indicator aims to measure the complexity of a mode. Its value ranges from 0 (no collinearity) to 1 (perfect collinearity). For real modes, mpc tends towards 1.

**Input**
- `ms`: Mode shapes (each column corresponds to a mode)

**Output**
- `mpc`: Mode phase collinearity values

**References**

[1] J.-N. Juang and R. Pappa: "An eigensystem realization algorithm for modal parameter identification and model reduction", Journal of Guidance, Control, and Dynamics, Vol. 8, No. 5, Sept.-Oct. 1985, pp. 620-627.

[2] M+P Analyzer manual. Rev. 5.1. 2017
"""
function mpc(ms)
    if ms isa Vector
        ms = reshape(ms, :, 1)
    end

    mpc = zeros(size(ms, 2))

    phi = ms .- mean(ms, dims = 1)

    phi_r = real(phi)
    phi_i = imag(phi)

    for (i, (pr, pi)) in enumerate(zip(eachcol(phi_r), eachcol(phi_i)))
        Sxx = pr'pr
        Syy = pi'pi
        Sxy = pr'pi

        λ1 = (Sxx + Syy)/2 + sqrt(((Sxx - Syy)/2)^2 + Sxy^2)
        λ2 = (Sxx + Syy)/2 - sqrt(((Sxx - Syy)/2)^2 + Sxy^2)

        mpc[i] = ((λ1 - λ2)/(λ1 + λ2))^2
    end

    return length(mpc) == 1 ? mpc[1] : mpc
end

"""
    mcf(ms)

Compute the mode complexity factor.

This indicator aims to measure the complexity of a mode. Its value ranges from 1 (purely real mode) to 0 (purely complex mode).

**Input**
- `ms`: Mode shapes (each column corresponds to a mode)

**Output**
- `mcf`: Mode complexity factor values

**References**

[1] J.-N. Juang and R. Pappa: "An eigensystem realization algorithm for modal parameter identification and model reduction", Journal of Guidance, Control, and Dynamics, Vol. 8, No. 5, Sept.-Oct. 1985, pp. 620-627.

[2] M+P Analyzer manual. Rev. 5.1. 2017
"""
function mcf(ms)
    if ms isa Vector
        ms = reshape(ms, :, 1)
    end

    mcf = zeros(size(ms, 2))

    phi_r = real(ms)
    phi_i = imag(ms)

    for (i, (pr, pi)) in enumerate(zip(eachcol(phi_r), eachcol(phi_i)))
        Sxx = pr'pr
        Syy = pi'pi
        Sxy = pr'pi

        mcf[i] = 1 - ((Sxx - Syy)^2 + 4Sxy^2)/(Sxx + Syy)^2
    end

    return length(mcf) == 1 ? mcf[1] : mcf
end

"""
    mpd(ms)

Compute the mode phase deviation.

This indicator aims to measure the complexity of a mode. Its value is close to 0 for real modes, and increases as the mode becomes more complex.

**Input**
- `ms`: Mode shapes (each column corresponds to a mode)

**Output**
- `mpd`: Mode phase deviation values

**Reference**

[1] E. Reynders, J. Houbrechts and G. De Roeck. Fully automated (operational) modal analysis. Mechanical Systems and Signal Processing. 29: 228-250. 2012

[2] A. C. Dederichs and O. Oiseth. Experimental comparison of automatic operational modal analysis algorithms for application to long-span road bridges. Mechanical Systems and Signal Processing. 199: 110485. 2023
"""
function mpd(ms)
    if ms isa Vector
        ms = reshape(ms, :, 1)
    end

    n, m = size(ms)
    mpd = zeros(m)

    S = similar(ms, n, 2)
    num = similar(ms, n)
    for (i, phi) in enumerate(eachcol(ms))
        S[:, 1] .= real(phi)
        S[:, 2] .= imag(phi)

        V = svd(S).V
        V12 = V[1, 2]
        V22 = V[2, 2]
        @. num = abs(phi)*acos((S[:, 1]*V22 - S[:, 2]*V12)/(sqrt(V12^2 + V22^2)*abs(phi)))
        den = sum(abs, phi)

        mpd[i] = real(sum(skipnan(num))/den)
    end

    return length(mpd) == 1 ? only(mpd) : mpd
end

## Correlation functions
"""
    msf(ms_exp, ms_ref)

Compute the modal scale factor between experimental and reference mode shapes.

**Inputs**
- `ms_exp`: Experimental mode shapes (nmes x nmodes array)
- `ms_ref`: Reference mode shapes (nmes x nmodes array)
**Output**
- `msf`: Modal scale factors (nmodes array)

**Reference**

[1] R. J. Allemang. The modal assurance criterion twenty years of use and abuse. Sound & Vibration. 37 (8): 14-23. 2003

[2] M. Rades. Comparison of vibration properties: Comparison of modal properties. Encyclopedia of Vibration. pp. 265-273. 2007.

**Note**

The modal scale factor is computed as the least squares error estimate of the proportionality constant between the corresponding elements of each modal vector
"""
function msf(ms_exp, ms_ref)
    if ms_exp isa Vector
        ms_exp = reshape(ms_exp, :, 1)
    end

    if ms_ref isa Vector
        ms_ref = reshape(ms_ref, :, 1)
    end

    ne, me = size(ms_exp)
    nr, mr = size(ms_ref)
    if (ne != nr) || (me != mr)
        throw(DimensionMismatchError("Experimental and reference mode shapes must have the same dimensions"))
    end

    msf = similar(ms_exp, me)
    for (i, (pe, pr)) in enumerate(zip(eachcol(ms_exp), eachcol(ms_ref)))
        num = pe'pr
        den = pe'pe
        msf[i] = num/den
    end

    return length(msf) == 1 ? only(msf) : msf
end

"""
    comac(ms_exp, ms_ref)

Compute the coordinate modal assurance criterion (COMAC) between experimental and reference mode shapes.

**Inputs**
- `ms_exp`: Experimental mode shapes (nmes x nmodes array)
- `ms_ref`: Reference mode shapes (nmes x nmodes array)
**Output**
- `comac`: Coordinate modal assurance criterion values (nmes array)

**Reference**

[1] R. J. Allemang. The modal assurance criterion twenty years of use and abuse. Sound & Vibration. 37 (8): 14-23. 2003
"""
function comac(ms_exp, ms_ref)
    if ms_exp isa Vector
        ms_exp = reshape(ms_exp, :, 1)
    end

    if ms_ref isa Vector
        ms_ref = reshape(ms_ref, :, 1)
    end

    ne, me = size(ms_exp)
    nr, mr = size(ms_ref)

    if (ne != nr) || (me != mr)
        throw(DimensionMismatch("Experimental and reference mode shapes must have the same dimensions"))
    end

    # Scale experimental mode shapes w.r.t. reference ones
    ms_exp_scaled = ms_exp .* msf(ms_exp, ms_ref)'

    comac = zeros(ne)
    for (p, (phi, psi)) in enumerate(zip(eachrow(ms_exp_scaled), eachrow(ms_ref)))
        num =  abs2(dot(phi, conj.(psi)))
        den = abs(dot(psi, conj.(psi))) * abs(dot(phi, conj.(phi)))
        comac[p] = real(num/den)
    end

    return length(comac) == 1 ? only(comac) : vec(comac)
end

"""
    ecomac(ms_exp, ms_ref)

Compute the enhanced coordinate modal assurance criterion (eCOMAC) between experimental and reference mode shapes.

**Inputs**
- `ms_exp`: Experimental mode shapes (nmes x nmodes array)
- `ms_ref`: Reference mode shapes (nmes x nmodes array)
**Output**
- `ecomac`: Enhanced coordinate modal assurance criterion values (nmes array)

**Reference**
[1] D. L. Hunt. Application of an Enhanced Coordinate Modal Assurance Criterion (ECOMAC). Proceedings of International Modal Analysis Conference, pp. 66-71, 1992.

[2] G. Martin, E. Balmes and T. Chancelier. Improved Modal Assurance Criterion using a quantification of identification errors per mode/sensor. Proceedings of ISMA 2014, pp. 2509-2519. 2014.
"""
function ecomac(ms_exp, ms_ref)
    if ms_exp isa Vector
        ms_exp = reshape(ms_exp, :, 1)
    end

    if ms_ref isa Vector
        ms_ref = reshape(ms_ref, :, 1)
    end

    ne, me = size(ms_exp)
    nr, mr = size(ms_ref)

    if (ne != nr) || (me != mr)
        throw(DimensionMismatch("Experimental and reference mode shapes must have the same dimensions"))
    end

    # Scale experimental mode shapes w.r.t. reference ones
    ms_exp_scaled = ms_exp .* msf(ms_exp, ms_ref)'
    ecomac = mean(abs, ms_ref .- ms_exp_scaled, dims = 2)/2

    return length(ecomac) == 1 ? only(ecomac) : vec(ecomac)
end

"""
    mac(ms_exp, ms_ref)

Compute the modal assurance criterion (MAC) between experimental and reference mode shapes.

**Inputs**
- `ms_exp`: Experimental mode shapes (nmes x nmodes array) if `ms_exp` is a matrix
- `ms_ref`: Reference mode shapes (nmes x nmodes array) if `ms_ref` is a matrix

**Output**
- `mac`: Modal assurance criterion values (nmodes x nmodes array)

**Reference**

[1] R. J. Allemang. The modal assurance criterion twenty years of use and abuse. Sound & Vibration. 37 (8): 14-23. 2003
"""
function mac(ms_exp, ms_ref)
    if ms_exp isa Vector
        ms_exp = reshape(ms_exp, :, 1)
    end

    if ms_ref isa Vector
        ms_ref = reshape(ms_ref, :, 1)
    end

    ne, me = size(ms_exp)
    nr, mr = size(ms_ref)

    if (ne != nr)
        throw(DimensionMismatch("The number of degrees of freedom of the experimental and reference mode shapes must be the same"))
    end

    mac = zeros(mr, me)
    for i in 1:mr, j in 1:me
        num = abs2(ms_exp[:, j]'ms_ref[:, i])
        den = (ms_exp[:, j]'ms_exp[:, j]) * (ms_ref[:, i]'ms_ref[:, i])
        mac[i, j] = num/real(den)
    end

    return sum(size(mac)) == 2 ? only(mac) : mac
end

"""
    frac(frf_exp, frf_ref)

Compute the frequency response assurance criterion (FRAC) between experimental and reference frequency response functions.

**Inputs**
- `frf_exp`: Experimental frequency response functions (nmes x nexc x nf array
- `frf_ref`: Reference frequency response functions (nmes x nexc x nf array)

**Output**
- `frac`: Frequency response assurance criterion values (nmes x nexc array)

**Reference**

[1] R. J. Allemang. The modal assurance criterion twenty years of use and abuse. Sound & Vibration. 37 (8): 14-23. 2003
"""
function frac(frf_exp, frf_ref)
    ne, me, nfe = size(frf_exp)
    nr, mr, nfr = size(frf_ref)

    if (ne != nr) || (me != mr) || (nfe != nfr)
        throw(DimensionMismatch("Experimental and reference FRFs must have the same dimensions"))
    end

    frac = zeros(ne, me)
    for i in 1:ne, j in 1:me
        num = abs2(frf_exp[i, j, :]'frf_ref[i, j, :])
        den = (frf_exp[i, j, :]'frf_exp[i, j, :]) * (frf_ref[i, j, :]'frf_ref[i, j, :])
        frac[i, j] = num/real(den)
    end

    return ne == 1 || me == 1 ? vec(frac) : frac
end

## Indicator functions
"""
    cmif(frf; type = :dis)

Compute the complex mode indicator function (CMIF)

**Inputs**
- `frf`: Frequency response functions (nmes x nexc x nf array)
- `type`: Type of FRF (:dis for displacement, :vel for velocity, :acc for acceleration)

**Output**
- `cmif`: Complex mode indicator function (min(nmes, nexc) x nf array)

**Reference**

[1] R. J. Allemang and D. L. Brown. A Complete Review of the Complex Mode  Indicator Function (CMIF) with Applications. ISMA 2006. 2006
"""
function cmif(frf; type = :dis)
    nmes, nexc, nf = size(frf)
    cmif = zeros(min(nmes, nexc), nf)

    for f in 1:nf
        if type == :dis || type == :acc
            cmif[:, f] = svd(imag(frf[:, :, f])).S
        elseif type == :vel
            cmif[:, f] = svd(real(frf[:, :, f])).S
        end
    end

    return size(cmif, 1) == 1 ? vec(cmif) : cmif
end

"""
    psif(frf)

Compute the power spectrum indicator function (PSIF).

**Input**
- `frf`: Frequency response functions (nmes x nexc x nf array or nmes x nf array)

**Output**
- `psif`: Power spectrum indicator function (nf array)

**Reference**

[1] M+P Analyzer manual. Rev. 5.1. 2017
"""
function psif(frf)
    n = ndims(frf)

    return vec(sum(abs2, frf, dims = 1:n-1))
end