"""
    to_dB(x, ref = 1.; type = :amplitude)

Convert a signal `x` to decibels (dB) with respect to the reference value `ref`.

**Inputs**
* `x`: Signal to convert (scalar or array, real or complex)
* `ref`: Reference value, expressed in the same unit as `x` (default: 1.)
* `type`: Type of quantity contained in `x`
    * `:amplitude` (default): Root-power (field) quantity - displacement, velocity, acceleration, pressure, FRF, RMS values, ...

        L = 20 log₁₀(|x|/ref)

    * `:power`: Power (or energy) quantity - squared amplitude, power, energy, PSD, autospectrum, ...

        L = 10 log₁₀(|x|/ref)

**Output**
* `L`: Level of `x` in dB

**Notes**

Both definitions give the same level for a given physical quantity, provided that `ref` is consistent with `x`. For a power quantity, `ref` must therefore be a power reference, i.e. the square of the amplitude reference (e.g. `ref = (2e-5)^2` Pa² for a pressure autospectrum).

Usual reference values are given by `ref_dB`.

Null values of `x` give `-Inf` dB.
"""
function to_dB(x, ref = 1.; type = :amplitude)
    ref > 0 || throw(ArgumentError("The reference value must be strictly positive."))

    if type === :amplitude
        factor = 20
    elseif type === :power
        factor = 10
    else
        throw(ArgumentError("type must be :amplitude or :power."))
    end

    return @. factor*log10(abs(x)/ref)
end

"""
    ref_dB(quantity::Symbol; type = :amplitude)

Return the usual reference value used to express a `quantity` in decibels (SI units).

**Inputs**
* `quantity`: Type of physical quantity
    * `:mobility`: 1 m/s/N
    * `:accelerance`: 1 m/s²/N
    * `:displacement`: 1 pm = 1e-12 m
    * `:velocity`: 1 nm/s = 1e-9 m/s
    * `:acceleration`: 1 µm/s² = 1e-6 m/s²
    * `:pressure`: 20 µPa = 2e-5 Pa
* `type`: Type of quantity to be converted (see `to_dB`)
    * `:amplitude` (default): Amplitude reference value
    * `:power`: Power reference value, i.e. the square of the amplitude reference value

**Output**
* `ref`: Reference value
```
"""
function ref_dB(quantity::Symbol; type = :amplitude)
    if quantity === :admittance
        ref = 1.
    elseif quantity === :mobility
        ref = 1.
    elseif quantity === :accelerance
        ref = 1.
    elseif quantity === :displacement
        ref = 1e-12
    elseif quantity === :velocity
        ref = 1e-9
    elseif quantity === :acceleration
        ref = 1e-6
    elseif quantity === :pressure
        ref = 2e-5
    else
        throw(ArgumentError("Unknown quantity $quantity. Available quantities are :mobility, :accelerance, :displacement, :velocity, :acceleration and :pressure."))
    end

    if type === :amplitude
        return ref
    elseif type === :power
        return ref^2
    else
        throw(ArgumentError("type must be :amplitude or :power."))
    end
end
