"""
    peak_value(x; dims = :)

Compute the peak value of a signal `x`.

The peak value is defined as:

x_peak = max(|x(t)|)

**Inputs**
* `x`: Signal (displacement, velocity, acceleration, ...)
* `dims`: Dimension along which the peak value is computed (default: all dimensions)

**Output**
* `x_peak`: Peak value

**Note**

The peak value is often used to characterize short-duration phenomena (e.g. shocks), but it does not account for the time evolution of the vibration.
"""
peak_value(x; dims = :) = maximum(abs, x; dims)

"""
    peak2peak(x; dims = :)

Compute the peak-to-peak value of a signal `x`.

The peak-to-peak value is defined as:

x_pp = max(x(t)) - min(x(t))

**Inputs**
* `x`: Real signal (displacement, velocity, acceleration, ...)
* `dims`: Dimension along which the peak-to-peak value is computed (default: all dimensions)

**Output**
* `x_pp`: Peak-to-peak value

**Note**

The peak-to-peak value is often used for displacement measurements, to be compared with a critical displacement.
"""
peak2peak(x; dims = :) = maximum(x; dims) .- minimum(x; dims)

"""
    rms_value(x; dims = :)

Compute the root mean square (RMS) value of a signal `x`.

The RMS value is defined as:

x_RMS = sqrt((1/T)∫₀ᵀ x²(t) dt)

For a sampled signal, it is computed as:

x_RMS = sqrt((1/N)∑ₙ |xₙ|²)

**Inputs**
* `x`: Signal (displacement, velocity, acceleration, ...)
* `dims`: Dimension along which the RMS value is computed (default: all dimensions)

**Output**
* `x_rms`: RMS value

**Note**

The RMS value is related to the vibratory energy and accounts for the time evolution of the signal.
"""
rms_value(x; dims = :) = sqrt.(mean(abs2, x; dims))
