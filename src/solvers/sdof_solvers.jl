"""
    SdofFreeTimeProblem(sdof, u0, t)

Structure containing the data of a time problem for a sdof system

**Fields**
* `sdof::Sdof`: Sdof structure
* `u0::AbstractVector`: Initial conditions
    * `x0`: Initial displacement [m]
    * `v0`: Initial velocity [m/s]
* `t::AbstractRange`: Time points at which to evaluate the response
"""
@show_data struct SdofFreeTimeProblem{Tu <: AbstractVector, Tt <: AbstractRange}
    sdof::Sdof
    u0::Tu
    t::Tt

    SdofFreeTimeProblem(sdof, u0::Tu, t::Tt) where {Tu, Tt} = new{Tu, Tt}(sdof, u0, t)
end

"""
    SdofHarmonicTimeProblem(sdof, u0, t, F, ω, type_exc)

Structure containing the data of a time problem for a sdof system subject to a harmonic excitation

**Constructor parameters**
* `sdof::Sdof`: Sdof structure
* `F`: Amplitude of the force excitation [N] or base motion [m]
* `f`: Frequency of the excitation [Hz]
* `u0::AbstractVector`: Initial conditions
    * `x0`: Initial displacement [m]
    * `v0`: Initial velocity [m/s]
* `t::AbstractRange`: Time points at which to evaluate the response
* `type_exc::Symbol`: Type of excitation
    * `:force`: External force (default)
    * `:base`: Base motion

**Fields**
* `sdof::Sdof`: Sdof structure
* `F`: Amplitude of the force excitation [N] or base motion [m]
* `ω`: Frequency of the excitation [rad/s]
* `u0::AbstractVector`: Initial conditions
    * `x0`: Initial displacement [m]
    * `v0`: Initial velocity [m/s]
* `t::AbstractRange`: Time points at which to evaluate the response
* `type_exc::Symbol`: Type of excitation
    * `:force`: External force (default)
    * `:base`: Base motion
"""
@show_data struct SdofHarmonicTimeProblem{Tu <: AbstractVector, Tt <: AbstractRange, T <: Real}
    sdof::Sdof
    u0::Tu
    t::Tt
    F::T
    ω::T
    type_exc::Symbol

    SdofHarmonicTimeProblem(sdof, u0::Tu, t::Tt, F::T, f::T, type_exc = :force) where {Tu, Tt, T} = new{Tu, Tt, T}(sdof, u0, t, F, 2π*f, type_exc)
end

"""
    SdofForcedTimeProblem(sdof, u0, t, F, ω, type_exc)

Structure containing the data of a time problem for a sdof system subject to an arbitrary excitation

**Fields**
* `sdof::Sdof`: Sdof structure
* `F::AbstractVector`: Amplitude of the force excitation [N] or base motion [m]
* `u0::AbstractVector`: Initial conditions
    * `x0`: Initial displacement [m]
    * `v0`: Initial velocity [m/s]
* `t::AbstractRange`: Time points at which to evaluate the response
* `type_exc::Symbol`: Type of excitation
    * `:force`: External force (default)
    * `:base`: Base motion
"""
@show_data struct SdofForcedTimeProblem{Tu <: AbstractVector, Tt <: AbstractRange, Tf <: AbstractVector}
    sdof::Sdof
    u0::Tu
    t::Tt
    F::Tf
    type_exc::Symbol

    SdofForcedTimeProblem(sdof, u0::Tu, t::Tt, F::Tf, type_exc = :force) where {Tu, Tt, Tf} = new{Tu, Tt, Tf}(sdof, u0, t, F, type_exc)
end

"""
    SdofTimeSolution(u, du, ddu)

Structure containing the data of the solution of the forced response of a sdof system

**Fields**
* `u::Vector{Real}`: Displacement solution
* `du::Vector{Real}`: Velocity solution
* `ddu::Vector{Real}`: Acceleration solution
"""
@show_data struct SdofTimeSolution{T <: Real}
    u::Vector{T}
    du::Vector{T}
    ddu::Vector{T}
end

"""
    SdofFRFProblem(sdof, freq, type_exc, type_resp)

Structure containing the data for computing the FRF a sdof system

**Fields**
* `sdof::Sdof`: Sdof structure
* `freq::AbstractRange`: Vector of frequencies [Hz]
* `type_exc::Symbol`: Type of excitation
    * `:force`: External force (default)
    * `:base`: Base motion
* `type_resp::Symbol`: Type of response
    * `:dis`: Admittance (default)
    * `:vel`: Mobility
    * `:acc`: Accelerance
"""
@show_data struct SdofFRFProblem{Tf <: AbstractRange}
    sdof::Sdof
    freq::Tf
    type_exc::Symbol
    type_resp::Symbol

    SdofFRFProblem(sdof, freq::Tf; type_exc = :force, type_resp = :dis) where Tf = new{Tf}(sdof, freq, type_exc, type_resp)
end

"""
    SdofFrequencyProblem(sdof, F, type_exc, type_resp)

Structure containing the data for computing the frequency response of a sdof system

**Fields**
* `sdof::Sdof`: Sdof structure
* `freq::AbstractRange`: Vector of frequencies [Hz]
* `F::AbstractVector`: Vector of the force excitation [N] or base motion [m]
* `type_exc::Symbol`: Type of excitation
    * `:force`: External force (default)
    * `:base`: Base motion
* `type_resp::Symbol`: Type of response
    * `:dis`: Displacement spectrum(default)
    * `:vel`: Velocity spectrum
    * `:acc`: Acceleration spectrum
"""
@show_data struct SdofFrequencyProblem{Tf <: AbstractRange, TF <: AbstractVector}
    sdof::Sdof
    freq::Tf
    F::TF
    type_exc::Symbol
    type_resp::Symbol

    SdofFrequencyProblem(sdof, freq::Tf, F::TF; type_exc = :force, type_resp = :dis) where {Tf, TF} = new{Tf, TF}(sdof, freq, F, type_exc, type_resp)
end

"""
    SdofFrequencySolution(u)

Structure containing the data of the solution of a frequency problem for a sdof system

**Field**
* `u::Vector{Complex}`: Solution of the frequency problem
   * Response spectrum (displacement, velocity, acceleration) [m, m/s, m/s²]
   * Or Frequency response function (FRF) (Admittance, Mobility, Accelerance) [m/N, m.s/N, m.s²/N]
"""
@show_data struct SdofFrequencySolution{T <: Complex}
    u::Vector{T}
end

"""
    solve(prob::SdofFreeTimeProblem)

Compute the free response of a single degree of freedom (Sdof) system.

**Input**
* `prob`: Structure containing the parameters of the Sdof problem

**Output**
* `sol`: The response of the system at the given time points
    * `u`: Displacement solution
    * `du`: Velocity solution
    * `ddu`: Acceleration solution
"""
function solve(prob::SdofFreeTimeProblem)
    (; sdof, u0, t) = prob
    (; ω0, ξ) = sdof
    x0, v0 = u0

    cache = (
        xh = similar(t),
        vh = similar(t),
        ah = similar(t)
    )

    free_response_sdof!(cache, ω0, ξ, x0, v0, t)

    return SdofTimeSolution(cache.xh, cache.vh, cache.ah)
end

"""
    free_response_sdof!(cache, ω0, ξ, x0, v0, t)

Compute the free response of a single degree of freedom (Sdof) system.

**Inputs**
* `cache`: Cache for the solution
    * `xh`: Displacement
    * `vh`: Velocity
    * `ah`: Acceleration
* `ω0`: Natural angular frequency [rad/s]
* `ξ`: Damping ratio
* `x0`: Initial displacement [m]
* `v0`: Initial velocity [m/s]
* `t`: Time vector
"""
function free_response_sdof!(cache, ω0, ξ, x0, v0, t)
    if ω0 == 0.
        @. cache.xh = x0 + v0*t
        cache.vh .= v0*one(cache.x)
        cache.ah .= zero(cache.x)
    elseif ξ < 1.
        Ω0 = ω0*√(1 - ξ^2)
        A = x0
        B = (v0 + ξ*ω0*x0)/Ω0

        @. cache.xh = (A*cos(Ω0*t) + B*sin(Ω0*t))*exp(-ξ*ω0*t)
        @. cache.vh = Ω0*(-A*sin(Ω0*t) + B*cos(Ω0*t))*exp(-ξ*ω0*t) - ξ*ω0*cache.xh
        @. cache.ah = -2ξ*ω0*cache.vh - ω0^2*cache.xh

    elseif ξ == 1.
        A = x0
        B = v0 + ω0*x0

        @. cache.xh = (A + B*t)*exp(-ω0*t)
        @. cache.vh = B*exp(-ω0*t) - ω0*cache.xh
        @. cache.ah = -2ω0*cache.vh - ω0^2*cache.xh
    else
        β = ω0*√(ξ^2 - 1)
        A = x0
        B = (v0 + ξ*ω0*x0)/β

        @. cache.xh = exp(-ξ*ω0*t)*(A*cosh(β*t) + B*sinh(β*t))
        @. cache.vh = β*(A*sinh(β*t) + B*cosh(β*t))*exp(-ξ*ω0*t) - ξ*ω0*cache.xh
        @. cache.ah = -2ξ*ω0*cache.vh - ω0^2*cache.xh
    end
end

"""
    solve(prob::SdofHarmonicTimeProblem)

Computes the forced response of a single degree of freedom (Sdof) system due to an harmonic external force or base motion

**Inputs**
* `prob`: Structure containing the parameters of the Sdof harmonic problem

**Output**
* `sol`: The response of the system at the given time points
    * `u`: Displacement solution
    * `du`: Velocity solution
    * `ddu`: Acceleration solution
"""
function solve(prob::SdofHarmonicTimeProblem)
    (; sdof, u0, t, F, ω, type_exc) = prob
    (; m, ω0, ξ) = sdof
    x0, v0 = u0

    cache = (
        x = similar(t),
        v = similar(t),
        a = similar(t),
        xh = similar(t),
        vh = similar(t),
        ah = similar(t)
    )

    harmonic_response_sdof!(cache, F, ω, x0, v0, m, ω0, ξ, t, type_exc)

    return SdofTimeSolution(cache.x, cache.v, cache.a)
end

"""
    harmonic_response_sdof!(cache, xh, vh, ah, F, ω, x0, v0, m, ω0, ξ, t, type_exc)

Compute the forced response of a single degree of freedom (Sdof) system due to an harmonic external force or base motion

**Inputs**
* `cache`: Cache for the solution
    * `x`: Displacement
    * `v`: Velocity
    * `a`: Acceleration
    * `xh`: Displacement - Homogeneous solution (cache)
    * `vh`: Velocity - Homogeneous solution (cache)
    * `ah`: Acceleration - Homogeneous solution (cache)
* `F`: Amplitude of the force excitation [N] or base motion [m]
* `ω`: Frequency of the excitation [rad/s]
* `x0`: Initial displacement [m]
* `v0`: Initial velocity [m/s]
* `m`: Mass [kg]
* `ω0`: Natural angular frequency [rad/s]
* `ξ`: Damping ratio
* `t::AbtractRange`: Time vector
* `type_exc::Symbol`: Type of excitation
    * `:force`: External force (default)
    * `:base`: Base motion
"""
function harmonic_response_sdof!(cache, F, ω, x0, v0, m, ω0, ξ, t, type_exc)
    if ξ == 0. && ω0 == ω
        # Variation parameters
        if type_exc == :force
            α = 1/m
        else
            α = ω^2
        end

        ρ0 = α*F/2ω
        A = x0
        B = v0/ω
        @. cache.x = A*cos(ω*t) + B*sin(ω*t) + ρ0*t*sin(ω*t)
        @. cache.v = -A*ω*sin(ω*t) + B*ω*cos(ω*t) + ρ0*(sin(ω*t) + ω*t*cos(ω*t))
        @. cache.a = -A*ω^2*cos(ω*t) - B*ω^2*sin(ω*t) + ρ0*(2ω*cos(ω*t) - ω^2*t*sin(ω*t))
    else
        if type_exc == :force
            X = F/m/(ω0^2 - ω^2 + 2im*ξ*ω*ω0)
        else
            X = F*(ω0^2 + 2im*ξ*ω*ω0)/(ω0^2 - ω^2 + 2im*ξ*ω*ω0)
        end

        ρ0 = abs.(X)
        ϕ = angle.(X)

        A = x0 - ρ0*cos(ϕ)
        if ω0 == 0.
            B = v0 + ρ0*ω*sin(ϕ)
            @. cache.xh = A + B*t
            @. cache.vh =  B*one(xh)
            @. cache.ah = zero(xh)
        elseif ξ < 1.
            Ω0 = ω0*√(1 - ξ^2)
            B = (v0 + ξ*ω0*A + ρ0*ω*sin(ϕ))/Ω0
            @. cache.xh = (A*cos(Ω0*t) + B*sin(Ω0*t))*exp(-ξ*ω0*t)
            @. cache.vh = Ω0*(-A*sin(Ω0*t) + B*cos(Ω0*t))*exp(-ξ*ω0*t) - ξ*ω0*cache.xh
            @. cache.ah = -2ξ*ω0*cache.vh - ω0^2*cache.xh
        elseif ξ == 1.
            B = v0 + ω0*A + ρ0*ω*sin(ϕ)
            @. cache.xh = (A + B*t)*exp(-ω0*t)
            @. cache.vh = B*exp(-ω0*t) - ω0*cache.xh
            @. cache.ah = -2ω0*cache.vh - ω0^2*cache.xh
        else
            β = ω0*√(ξ^2 - 1)
            B = (v0 + ξ*ω0*A + ρ0*ω*sin(ϕ))/β
            @. cache.xh = (A*cosh(β*t) + B*sinh(β*t))*exp(-ξ*ω0*t)
            @. cache.vh = β*(A*sinh(β*t) + B*cosh(β*t))*exp(-ξ*ω0*t) - ξ*ω0*cache.xh
            @. cache.ah = -2ξ*ω0*cache.vh - ω0^2*cache.xh
        end

        @. cache.x = cache.xh + ρ0*cos(ω*t + ϕ)
        @. cache.v = cache.vh - ρ0*ω*sin(ω*t + ϕ)
        @. cache.a = cache.ah - ρ0*ω^2*cos(ω*t + ϕ)
    end
end

"""
    solve(prob::SdofForcedTimeProblem)

Computes the forced response of a single degree of freedom (Sdof) system due to an arbitrary external force or base motion

**Inputs**
* `prob`: Structure containing the parameters of the Sdof forced problem
* `method`: Method to compute the Duhamel's integral
    * `:filt`: Filtering using the Z-transform of the impulse response (default)
    * `:interp`: Interpolation + Gaussian quadrature
    * `:conv`: Convolution product

**Output**
* `sol`: The response of the system at the given time points
    * `u`: Displacement solution
    * `du`: Velocity solution
    * `ddu`: Acceleration solution
"""
function solve(prob::SdofForcedTimeProblem; method = :filt)

    (; sdof, u0, t, F, type_exc) = prob
    (; m, ω0, ξ) = sdof
    x0, v0 = u0

    # Cache for the free response
    cache_h = (
        xh = similar(t),
        vh = similar(t),
        ah = similar(t),
    )

    # Cache for the particular solution
    cache_p = (
        xp = similar(t),
        vp = similar(t),
        ap = similar(t),
        h = similar(t),
        num = similar(t, 3),
        denom = similar(t, 3)
    )

    free_response_sdof!(cache_h, ω0, ξ, x0, v0, t)

    forced_response_sdof!(cache_p, F, m, ω0, ξ, t, type_exc, method)

    x = cache_h.xh .+ cache_p.xp
    v = cache_h.vh .+ cache_p.vp
    a = cache_h.ah .+ cache_p.ap

    return SdofTimeSolution(x, v, a)
end

"""
    forced_response_sdof!(cache, F, m, ω0, ξ, t, type_exc, method)

Compute the forced response of a single degree of freedom (Sdof) system due to an arbitrary external force or base motion

**Inputs**
* `cache`: Cache for the solution
    * `xp`: Displacement
    * `vp`: Velocity
    * `ap`: Acceleration
    * `h`: Impulse response
    * `num`: Numerator of the transfer function in z-domain
    * `denom`: Denominator of the transfer function in z-domain
* `F`: Amplitude of the force excitation [N] or base motion [m]
* `m`: Mass [kg]
* `ω0`: Natural angular frequency [rad/s]
* `ξ`: Damping ratio
* `t`: Time vector
* `type_exc`: Type of excitation
    * `:force`: External force (default)
    * `:base`: Base motion
* `method`: Method to compute the Duhamel's integral
    * `:filt`: Filtering using the Z-transform of the impulse response (default)
    * `:interp`: Interpolation + Gaussian quadrature
    * `:conv`: Convolution product
"""
function forced_response_sdof!(cache, F, m, ω0, ξ, t, type_exc, method)
    # Impulse response
    Δt = t[2] - t[1]
    if ω0 == 0.
        if method == :interp || method == :conv
            @. cache.h = t/m
        else
            cache.num .= [0., Δt/m, 0.]
            cache.denom .= [1., -2., 1.]
        end
    elseif ξ < 1.
        Ω0 = ω0*√(1 - ξ^2)
        if method == :interp || method == :conv
            @. cache.h = exp(-ξ*ω0*t)*sin(Ω0*t)/m/Ω0
        else
            α = exp(-ξ*ω0*Δt)
            β = Ω0*Δt
            # Transfer function in the z-domain
            cache.num .= [0., α*sin(β)/m/Ω0, 0.]
            cache.denom .= [1., -2*α*cos(β), α^2]
        end
    elseif ξ == 1.
        if method == :interp || method == :conv
            @. cache.h = t*exp(-ω0*t)/m
        else
            α = exp(-ω0*Δt)
            cache.num .= [0., α*Δt/m, 0.]
            cache.denom .= [1., -2*α, α^2]
        end
    else
        β = ω0*√(ξ^2 - 1)
        if method == :interp || method == :conv
            @. cache.h = exp(-ξ*ω0*t)*sinh(β*t)/m/β
        else
            α = exp(-ξ*ω0*Δt)
            γ = β*Δt
            cache.num .= [0., α*sinh(γ)/m/β, 0.]
            cache.denom .= [1., -2*α*cosh(γ), α^2]
        end
    end

    if type_exc == :base
        k, c = ω0^2*m, 2ξ*ω0*m
        xb = F
        vb = gradient(xb, t)

        F .= k*xb .+ c*vb
    end

    if method == :interp || method == :conv
        if method == :interp
            cache.xp .= duhamel_integral(F, cache.h, t)
        else
            cache.xp .= Δt*conv(F, cache.h)[1:length(F)]
        end
    else
        cache.xp .= Δt*filt(cache.num, cache.denom, F)
    end

    cache.vp .= gradient(cache.xp, t)
    cache.ap .= gradient(cache.vp, t)
end

"""
    duhamel_integral(F, h, t)

Compute the Duhamel's integral

**Inputs**
* `F`: Force excitation
* `h`: Impulse response
* `t`: Time vector

**Output**
* `x`: Solution of the Duhamel's integral
"""
function duhamel_integral(F, h, t)
    # Interpolation
    fc = cubic_spline_interpolation(t, F)
    hc = cubic_spline_interpolation(t, h)

    # Gauss-Legendre quadrature
    nodes, weights = gausslegendre(500)

    # Initialization
    x = similar(t)
    scaled_nodes = similar(nodes)
    scaled_weights = similar(nodes)
    for (i, ti) in enumerate(t)
        @. scaled_nodes = (nodes + 1)*ti/2
        @. scaled_weights = weights*ti/2
        x[i] = sum(fc(τ) * hc(ti - τ) * w for (τ, w) in zip(scaled_nodes, scaled_weights))
    end

    return x
end

"""
    solve(prob::SdofFRFProblem)

Compute the FRF of a single degree of freedom (Sdof) system

**Inputs**
* `prob`: Structure containing the parameters of the Sdof FRF problem

**Output**
* `sol`: Solution of the FRF problem
    * `u`: FRF of the system
"""
function solve(prob::SdofFRFProblem)
    (; sdof, freq, type_exc, type_resp) = prob
    (; m, ω0, ξ) = sdof
    ω = 2π*freq

    x = similar(ω, Complex{eltype(ω)})
    if type_exc == :force
        @. x =  1/m/(ω0^2 - ω^2 + 2im*ξ*ω0*ω)
    else
        @. x = (ω0^2 + 2im*ξ*ω0*ω)/(ω0^2 - ω^2 + 2im*ξ*ω0*ω)
    end

    if type_resp == :vel
        @. x *= 1im*ω
    elseif type_resp == :acc
        @. x *= -ω^2
    end

    return SdofFrequencySolution(x)
end

"""
    solve(prob::SdofFrequencyProblem)

Compute the frequency response function of a single degree of freedom (Sdof) system

**Inputs**
* `prob`: Structure containing the parameters of the Sdof frequency problem

**Output**
* `sol`: Solution of the problem
    * `u`: Response spectrum of the system
"""
function solve(prob::SdofFrequencyProblem)
    (; sdof, freq, F, type_exc, type_resp) = prob
    (; m, ω0, ξ) = sdof
    ω = 2π*freq

    x = similar(freq, Complex{Real})
    if type_exc == :force
        @. x = F/m/(ω0^2 - ω^2 + 2im*ξ*ω0*ω)
    else
        @. x = F*(ω0^2 + 2im*ξ*ω0*ω)/(ω0^2 - ω^2 + 2im*ξ*ω0*ω)
    end

    if type_resp == :vel
        @. x *= 1im*ω
    elseif type_resp == :acc
        @. x *= -ω^2
    end

    return SdofFrequencySolution(x)
end

"""
    impulse_response(sdof, t)

Compute the impulse response of a single degree of freedom (Sdof) system

**Inputs**
* `sdof`: Sdof structure
* `t`: Time points at which to evaluate the response

**Output**
* `h`: Impulse response
"""
function impulse_response(sdof, t)
    (; m, ω0, ξ) = sdof

    cache = (
        xh = similar(t),
        vh = similar(t),
        ah = similar(t)
    )

    free_response_sdof!(cache, ω0, ξ, zero(ξ), 1/m, t)

    return cache.xh
end

abstract type SRSalg end
struct BasicSRS <: SRSalg end
struct RecursiveInt <: SRSalg end

"""
    srs(base_acc::ArbitraryExc, freq, t, ξ = 0.05,
        type = (instance = :primary, amplitude = :abs), alg = :Smallwood)

Compute the Shock Response Spectrum (SRS)

**Inputs**
* `base_acc::ArbitraryExc`: Base acceleration type - see `excitation` function for detatils
* `freq`: Vector of frequencies [Hz]
* `t`: Time points at which to evaluate the response
* `ξ`: Damping ratio (default = 0.05)
* `type`: Type of SRS
    * `instance`: Instance of the SRS
        * `:primary`: Primary instance (default)
        * `:secondary`: Secondary instance
    * `amplitude`: Amplitude used of the computing SRS
        * `:abs`: Maximum absolute amplitude (default)
        * `:pos`: Maximum positive amplitude
        * `:neg`: Maximum negative amplitude
* `alg`: Algorithm to compute the SRS
    * `:Basic`: Basic algorithm
    * `:RecursiveInt`: Recursive integration algorithm
    * `:RecursiveFilt`: Recursive filtering algorithm
    * `:Smallwood`: Smallwood algorithm (default)

**Output**
* `srs`: Vector of the SRS values

*Note*
* Primary instance - response of the system during the application of the base acceleration
* Secondary instance - response of the system after the application of the base acceleration
"""
function srs(base_acc::ArbitraryExc, freq, t, ξ = 0.05; type_srs = (instance = :primary, amplitude = :abs), alg = :Smallwood)

    # Some checks
    t_srs = typeof(type_srs)
    if !hasfield(t_srs, :instance)
        type_instance = :primary
    else
        type_instance = type_srs.instance
    end

    if !hasfield(t_srs, :amplitude)
        type_amplitude = :abs
    else
        type_amplitude = type_srs.amplitude
    end

    # Inner function to compute the required amplitude
    amplitude = let
        if type_amplitude == :abs
            x ->  maximum(abs, x)
        elseif type_amplitude == :pos
            x -> maximum(x)
        else
            x -> abs(minimum(x))
        end
    end

    (; tstart, duration) = base_acc
    t_primary = tstart + duration

    pos_primary = argmin(@. abs(t - t_primary)^2)

    acc = excitation(base_acc, t)

    srs = similar(freq)
    x = similar(t)
    for (f, f0) in enumerate(freq)
        if alg == :Basic
            x .= srs_basic(f0, ξ, acc, t)
        elseif alg == :RecursiveInt
            x .= srs_recursive_int(f0, ξ, acc, t)
        elseif alg == :RecursiveFilt
            x .= srs_recursive_filt(f0, ξ, acc, t)
        elseif alg == :Smallwood
            x .= srs_smallwood(f0, ξ, acc, t)
        end

        if type_instance == :primary
            srs[f] = amplitude(x[1:pos_primary])
        else
            srs[f] = amplitude(x[pos_primary + 1:end])
        end
    end

    return srs
end

function srs_basic(f0, ξ, acc, t)
    ω0 = 2π*f0
    sdof = Sdof(1., f0, ξ)

    prob = SdofForcedTimeProblem(sdof, [0., 0.], t, -acc, :conv)
    (; u, du) = solve(prob)

    return -2ξ*ω0*du - ω0^2*u
end

function srs_recursive_int(f0, ξ, acc, t)
    # Time step
    h = t[2] - t[1]
    nt = length(t)

    if ξ ≥ 1.
        error("Damping ratio must be less than 1")
    end

    # Constants
    ω0 = 2π*f0
    Ω0 = ω0*√(1 - ξ^2)

    # Pre-compute integration constants
    eh = exp(-ξ*ω0*h)
    eCos = eh*cos(Ω0*h)
    eSin = eh*sin(Ω0*h)

    # Integration coefficients of Duhamel's integral
    c1 = eCos + ξ*eSin/(1 - ξ^2)
    c2 = eSin/Ω0
    c3 = (c1 - 1.)/ω0^2
    c4 = (2ξ*(1. - eCos)/ω0/h + (1. - 2ξ^2)*eSin/Ω0/h - 1.)/ω0^2
    c5 = (4ξ/ω0/h + ((2 - 8ξ^2)/ω0^2/h^2 - 2ξ/ω0/h)*(1. - eCos)
        - ((1. - 2ξ^2)/ω0/h + (6ξ- 8ξ^3)/ω0^2/h^2)*eSin/√(1 - ξ^2))/2ω0^2
    c6 = -c2*Ω0
    c7 = (eCos - ξ*eSin/√(1 - ξ^2))/ω0
    c8 = -c2/ω0
    c9 = (c1 - 1.)/ω0^3/h
    c10 = ((1/ω0/h + 4ξ/ω0^2/h^2)*(1. - eCos) + ((2. - 4ξ^2)/ω0^2/h^2 - ξ/ω0/h)*eSin/√(1 - ξ^2) - 2/ω0/h)/2ω0^2

    # Acceleration response
    x = similar(t)

    # Initial conditions - Relative displacement and velocity
    z = 0.
    dz = 0.
    jerk = acc[1]
    snap = acc[1]
    for k in 1:nt - 1
        z = c1*z + c2*dz + c3*acc[k] + c4*jerk + c5*snap
        dz = ω0*(c6*z + c7*dz + c8*acc[k] + c9*jerk + c10*snap)
        x[k] = -2ξ*ω0*dz - ω0^2*z

        jerk = acc[k + 1] - acc[k]
        if k > 1
            snap = acc[k+1] - 2*acc[k] + acc[k-1]
        else
            snap = acc[k+1] - 2acc[k]
        end
    end

    z = c1*z + c2*dz + c3*acc[end] + c4*jerk + c5*snap
    dz = ω0*(c6*z + c7*dz + c8*acc[end] + c9*jerk + c10*snap)
    x[end] = -2ξ*ω0*dz - ω0^2*z

    return x
end

function srs_recursive_filt(f0, ξ, acc, t)
    # Time step
    h = t[2] - t[1]

    if ξ ≥ 1.
        error("Damping ratio must be less than 1")
    end

    # Constants
    ω0 = 2π*f0
    Ω0 = ω0*√(1 - ξ^2)

    # Pre-compute integration constants
    eh = exp(-ξ*ω0*h)
    eCos = eh*cos(Ω0*h)
    eSin = eh*sin(Ω0*h)

    # filter coefficients
    b = similar(t, 3)
    a = similar(b)

    b[1] = 2ξ*ω0*h
    b[2] = ω0*h*(ω0*(1. -2ξ^2)*eSin/Ω0 - 2ξ*eCos)
    b[3] = 0.

    a[1] = 1.
    a[2] = -2eCos
    a[3] = exp(-2ξ*ω0*h)

    return filt(b, a, acc)
end

function srs_smallwood(f0, ξ, acc, t)
    # Time step
    h = t[2] - t[1]

    if ξ ≥ 1.
        error("Damping ratio must be less than 1")
    end

    # Constants
    ω0 = 2π*f0
    Ω0 = ω0*√(1 - ξ^2)

    # Pre-compute integration constants
    eh = exp(-ξ*ω0*h)
    eCos = eh*cos(Ω0*h)
    eSin = eh*sin(Ω0*h)

    # filter coefficients
    b = similar(t, 3)
    a = similar(b)

    b[1] = 1 - eSin/Ω0/h
    b[2] = 2(eSin/Ω0/h - eCos)
    b[3] = exp(-2ξ*ω0*h) - eSin/Ω0/h

    a[1] = 1.
    a[2] = -2eCos
    a[3] = exp(-2ξ*ω0*h)

    return filt(b, a, acc)
end