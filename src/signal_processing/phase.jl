"""
    unwrap_phase(x; dims = 1, deg = false)

Compute the unwrapped phase of a signal `x`.

The phase returned by `angle` is wrapped in the interval (-π, π]. This function removes the artificial 2π jumps to obtain a continuous phase.

**Inputs**
* `x`: Signal
    * Complex signal (e.g. FRF, spectrum): the phase is first computed with `angle`
    * Real signal: `x` is assumed to be a wrapped phase expressed in **radians**
* `dims`: Dimension along which the phase is unwrapped (default: 1, i.e. each column is considered as a channel)
* `deg`: If `true`, the unwrapped phase is returned in degrees (default: `false`)

**Output**
* `ϕ`: Unwrapped phase [rad or deg]

**Notes**

Common pitfalls avoided by this function:
* The unwrapping is always performed in radians. Converting the phase to degrees **before** unwrapping it with a 2π range leads to wrong results.
* For multichannel signals, the unwrapping is performed independently along `dims`, so that the channels are not mixed.

The unwrapping relies on the assumption that the phase difference between two successive samples is lower than π. If the frequency resolution is too coarse, the unwrapped phase may be wrong.
"""
function unwrap_phase(x::AbstractArray; dims = 1, deg = false)
    ϕ = eltype(x) <: Complex ? angle.(x) : float.(x)
    unwrap!(ϕ; dims)

    return deg ? rad2deg.(ϕ) : ϕ
end
