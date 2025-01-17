"""
    Sdof(m, ω₀, ξ)

Structure containing the data of a sdof system

# Fields
* `m`: Mass [kg]
* `ω₀`: natural angular frequency [rad/s]
* `ξ`: Damping ratio
"""
@with_kw struct Sdof
    m :: Float64
    ω₀ ::Float64
    ξ :: Float64
end