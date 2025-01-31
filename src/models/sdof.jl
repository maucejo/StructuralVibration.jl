"""
    Sdof(m, ω₀, ξ)

Structure containing the data of a sdof system

# Constructor
* `m`: Mass [kg]
* `f₀`: Natural frequency [rad/s]
* `ξ`: Damping ratio

# Fields
* `m`: Mass [kg]
* `ω₀`: Natural frequency [rad/s]
* `ξ`: Damping ratio
"""
@with_kw struct Sdof
    m :: Float64
    f₀ ::Float64
    ξ :: Float64

    Sdof(m, f₀, ξ) = new(m, 2π*f₀, ξ)
end