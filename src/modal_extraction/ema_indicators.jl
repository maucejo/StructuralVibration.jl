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

    return mof
end

## Mode complexity
"""
    mov(poles, ms, Q)

Compute the mode overcomplexity value

This indicator is a weighted percentage of the degrees of freedom of the response for which adding mass leads to a negative frequency shift. A value close to 1 indicates physical modes, while a low value indicates non-physical modes (numerical or noise-related)

**Inputs**
- `poles`: Poles of the system
- `ms`: Mode shapes (each column corresponds to a mode)
- `Q`: Scaling factors associated to the mode shapes

**Output**
- `mov`: Mode overcomplexity values

**Reference**

[1] M+P Analyzer manual. Rev. 5.1. 2017
"""
function mov(poles, ms, Q)
    n, m = size(ms)
    mov = zeros(m)
    s = zeros(n)
    Msen = similar(s)
    sig = similar(s)

    for (i, (p, phi, Qi)) in enumerate(zip(poles, eachcol(ms), Q))
        @. Msen = imag(-p^2*phi^2*Qi) # Poles sensitivity w.r.t. mass
        @. sig = sign(Msen)
        s[sig .> 0.] .= 0.
        s[sig .<= 0.] .= 1.

        mov[i] = sum(s.*abs2.(phi))/sum(abs2, phi)
    end

    return mov
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
    mpc = zeros(size(ms, 2))

    phi = ms .- mean(ms, dims = 1)

    phi_r = real(phi)
    phi_i = imag(phi)

    for (i, (pr, pi)) in enumerate(zip(eachcol(phi_r), eachcol(phi_i)))
        Sxx = pr'pr
        Syy = pi'pi
        Sxy = pr'pi

        λ1 = (Sxx + Syy)/2 + sqrt(((Sxx - Syy)/2)^4 + Sxy^2)
        λ2 = (Sxx + Syy)/2 - sqrt(((Sxx - Syy)/2)^4 + Sxy^2)

        mpc[i] = (2(abs(λ1)/(abs(λ1) + abs(λ2))) - 0.5)^2
    end

    return mpc
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
    mcf = zeros(size(ms, 2))

    phi_r = real(ms)
    phi_i = imag(ms)

    for (i, (pr, pi)) in enumerate(zip(eachcol(phi_r), eachcol(phi_i)))
        Sxx = pr'pr
        Syy = pi'pi
        Sxy = pr'pi

        mcf[i] = 1 - ((Sxx - Syy)^2 + 4Sxy^2)/(Sxx + Syy)^2
    end

    return mcf
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

        mpd[i] = sum(num)/den
    end

    return mpd
end

## Correlation functions
"""
    msf(ms_exp, ms_th)

Compute the modal scale factor between experimental and theoretical mode shapes.

**Inputs**
- `ms_exp`: Experimental mode shapes (nmes x nmodes array)
- `ms_th`: Theoretical mode shapes (nmes x nmodes array)

**Output**
- `msf`: Modal scale factors (nmodes array)

**Reference**

[1] R. J. Allemang. The modal assurance criterion twenty years of use and abuse. Sound & Vibration. 37 (8): 14-23. 2003
"""
function msf(ms_exp, ms_th)
    ne, me = size(ms_exp)
    nt, mt = size(ms_th)

    if (ne != nt) || (me != mt)
        throw(DimensionMismatchError("Experimental and theoretical mode shapes must have the same dimensions"))
    end

    msf = zeros(me)
    for (i, (pe, pt)) in enumerate(zip(eachcol(ms_exp), eachcol(ms_th)))

        num = pt'conj(pe)
        den = pt'conj(pt)
        msf[i] = num/den
    end

    return msf
end

"""
    comac(ms_exp, ms_th)

Compute the coordinate modal assurance criterion (COMAC) between experimental and theoretical mode shapes.

**Inputs**
- `ms_exp`: Experimental mode shapes (nmes x nmodes array)
- `ms_th`: Theoretical mode shapes (nmes x nmodes array)

**Output**
- `comac`: Coordinate modal assurance criterion values (nmes array)

**Reference**

[1] R. J. Allemang. The modal assurance criterion twenty years of use and abuse. Sound & Vibration. 37 (8): 14-23. 2003
"""
function comac(ms_exp, ms_th)
    ne, me = size(ms_exp)
    nt, mt = size(ms_th)

    if (ne != nt) || (me != mt)
        throw(DimensionMismatch("Experimental and theoretical mode shapes must have the same dimensions"))
    end

    # Scale experimental mode shapes w.r.t. theoretical ones
    ms_exp .*= transpose(msf(ms_exp, ms_th))

    num = sum(abs2, ms_exp.*ms_th, dims = 2)
    den = sum(abs2, ms_th, dims = 2) .* sum(abs2, ms_exp, dims = 2)

    return num./den
end


"""
    ecomac(ms_exp, ms_th)

Compute the enhanced coordinate modal assurance criterion (eCOMAC) between experimental and theoretical mode shapes.

**Inputs**
- `ms_exp`: Experimental mode shapes (nmes x nmodes array)
- `ms_th`: Theoretical mode shapes (nmes x nmodes array)

**Output**
- `ecomac`: Enhanced coordinate modal assurance criterion values (nmes array)

**Reference**
[1] D. L. Hunt. Application of an Enhanced Coordinate Modal Assurance Criterion (ECOMAC). Proceedings of International Modal Analysis Conference, pp. 66-71, 1992.

[2] G. Martin, E. Balmes and T. Chancelier. Improved Modal Assurance Criterion using a quantification of identification errors per mode/sensor. Proceedings of ISMA 2014, pp. 2509-2519. 2014.
"""
function ecomac(ms_exp, ms_th)
    ne, me = size(ms_exp)
    nt, mt = size(ms_th)

    if (ne != nt) || (me != mt)
        throw(DimensionMismatch("Experimental and theoretical mode shapes must have the same dimensions"))
    end

    # Scale experimental mode shapes w.r.t. theoretical ones
    ms_exp .*= transpose(msf(ms_exp, ms_th))

    return mean(abs, ms_th .- ms_exp, dims = 2)/2.
end

"""
    mac(ms_exp, ms_th)

Compute the modal assurance criterion (MAC) between experimental and theoretical mode shapes.

**Inputs**
- `ms_exp`: Experimental mode shapes (nmes x nmodes array)
- `ms_th`: Theoretical mode shapes (nmes x nmodes array)

**Output**
- `mac`: Modal assurance criterion values (nmodes x nmodes array)

**Reference**

[1] R. J. Allemang. The modal assurance criterion twenty years of use and abuse. Sound & Vibration. 37 (8): 14-23. 2003
"""
function mac(ms_exp, ms_th)
    ne, me = size(ms_exp)
    nt, mt = size(ms_th)

    if (ne != nt)
        throw(DimensionMismatch("The number of degrees of freedom of the experimental and theoretical mode shapes must be the same"))
    end

    mac = zeros(nt, ne)
    for i in 1:nt, j in 1:ne
        num = abs2(ms_exp[:, j]'ms_th[:, i])
        den = (ms_exp[:, j]'ms_exp[:, j]) * (ms_th[:, i]'ms_th[:, i])
        mac[i, j] = num/real(den)
    end


    return mac
end


"""
    frac(frf_exp, frf_th)

Compute the frequency response assurance criterion (FRAC) between experimental and theoretical frequency response functions.

**Inputs**
- `frf_exp`: Experimental frequency response functions (nmes x nexc x nf array
- `frf_th`: Theoretical frequency response functions (nmes x nexc x nf array)

**Output**
- `frac`: Frequency response assurance criterion values (nmes x nexc array)

**Reference**

[1] R. J. Allemang. The modal assurance criterion twenty years of use and abuse. Sound & Vibration. 37 (8): 14-23. 2003
"""
function frac(frf_exp, frf_th)
    ne, me, nfe = size(frf_exp)
    nt, mt, nft = size(frf_th)

    if (ne != nt) || (me != mt) || (nfe != nft)
        throw(DimensionMismatch("Experimental and theoretical FRFs must have the same dimensions"))
    end

    frac = zeros(ne, me)
    for i in 1:ne, j in 1:me
        num = abs2(frf_exp[i, j, :]'frf_th[i, j, :])
        den = (frf_exp[i, j, :]'frf_exp[i, j, :]) * (frf_th[i, j, :]'frf_th[i, j, :])
        frac[i, j] = num/real(den)
    end

    return frac
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

    return cmif
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

    return sum(abs2, frf, dims = 1:n-1)[:]
end