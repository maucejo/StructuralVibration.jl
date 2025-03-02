"""
    SdofFreeTimeProblem(sdof, u0, t, type_exc)

Structure containing the data of a time problem for a sdof system

# Fields
* `sdof`: Sdof structure
* `u0`: Initial conditions
    * `x₀`: Initial displacement [m]
    * `v₀`: Initial velocity [m/s]
* `t`: Time points at which to evaluate the response
* `type_exc`: Type of excitation
    * :force: External force (default)
    * :base: Base motion
"""
@with_kw struct SdofFreeTimeProblem
    sdof :: Sdof
    u0 :: Vector{Float64}
    t
    type_exc :: Symbol

    SdofFreeTimeProblem(sdof, u0, t, type_exc = :force) = new(sdof, u0, t, type_exc)
end

"""
    SdofHarmonicTimeProblem(sdof, u0, t, F, ω, type_exc)

Structure containing the data of a time problem for a sdof system subject to a harmonic excitation

# Constructor
* `sdof`: Sdof structure
* `u0`: Initial conditions
    * `x₀`: Initial displacement [m]
    * `v₀`: Initial velocity [m/s]
* `t`: Time points at which to evaluate the response
* `F`: Amplitude of the force excitation [N] or base motion [m]
* `f`: Frequency of the excitation [Hz]
* `type_exc`: Type of excitation
    * :force: External force (default)
    * :base: Base motion

# Fields
* `sdof`: Sdof structure
* `u0`: Initial conditions
    * `x₀`: Initial displacement [m]
    * `v₀`: Initial velocity [m/s]
* `t`: Time points at which to evaluate the response
* `F`: Amplitude of the force excitation [N] or base motion [m]
* `ω`: Frequency of the excitation [rad/s]
* `type_exc`: Type of excitation
    * :force: External force (default)
    * :base: Base motion
"""
@with_kw struct SdofHarmonicTimeProblem
    sdof :: Sdof
    u0 :: Vector{Float64}
    t
    F :: Float64
    ω :: Float64
    type_exc :: Symbol

    SdofHarmonicTimeProblem(sdof, u0, t, F, f, type_exc = :force) = new(sdof, u0, t, F, 2π*f, type_exc)
end

"""
    SdofForcedTimeProblem(sdof, u0, t, F, ω, type_exc)

Structure containing the data of a time problem for a sdof system subject to an arbitrary excitation

# Fields
* `sdof`: Sdof structure
* `u0`: Initial conditions
    * `x₀`: Initial displacement [m]
    * `v₀`: Initial velocity [m/s]
* `t`: Time points at which to evaluate the response
* `F`: Amplitude of the force excitation [N] or base motion [m]
* `type_exc`: Type of excitation
    * :force: External force (default)
    * :base: Base motion
"""
@with_kw struct SdofForcedTimeProblem
    sdof :: Sdof
    u0 :: Vector{Float64}
    t
    F :: Vector{Float64}
    type_exc :: Symbol

    SdofForcedTimeProblem(sdof, u0, t, F, type_exc = :force) = new(sdof, u0, t, F, type_exc)
end

"""
    SdofTimeSolution(u, du, ddu)

Structure containing the data of the solution of the forced response of a sdof system

# Fields
* `u`: Displacement solution
* `du`: Velocity solution
* `ddu`: Acceleration solution
"""
@with_kw struct SdofTimeSolution
    u :: Vector{Float64}
    du :: Vector{Float64}
    ddu :: Vector{Float64}
end

"""
    SdofImpulseSolution(u)

Structure containing the data of the solution of the impulse response of a sdof system

# Field
* `u`: Solution of the impulse response
"""
@with_kw struct SdofImpulseSolution
    u :: Vector{Float64}
end

"""
    SdofFRFProblem(sdof, u0, t, F, type_exc, type_resp)

Structure containing the data for computing the FRF a sdof system

# Fields
* `sdof`: Sdof structure
* `freq``: Vector of frequencies [Hz]
* `type_exc`: Type of excitation
    * :force: External force (default)
    * :base: Base motion
* `type_resp`: Type of response
    * :dis: Displacement spectrum or Admittance (default)
    * :vel: Velocity spectrum or Mobility
    * :acc: Acceleration spectrum or Accelerance
"""
@with_kw struct SdofFRFProblem
    sdof :: Sdof
    freq
    type_exc :: Symbol
    type_resp :: Symbol

    SdofFRFProblem(sdof, freq; type_exc = :force, type_resp = :dis) = new(sdof, freq, type_exc, type_resp)
end

"""
    SdofFrequencyProblem(sdof, F, type_exc, type_resp)

Structure containing the data for computing the frequency response of a sdof system

# Fields
* `sdof`: Sdof structure
* `freq``: Vector of frequencies [Hz]
* `F`: Vector of the force excitation [N] or base motion [m]
* `type_exc`: Type of excitation
    - :force: External force (default)
    - :base: Base motion
* `type_resp`: Type of response
    - :dis: Displacement spectrum or Admittance (default)
    - :vel: Velocity spectrum or Mobility
    - :acc: Acceleration spectrum or Accelerance
"""
@with_kw struct SdofFrequencyProblem
    sdof :: Sdof
    freq
    F :: Vector{Float64}
    type_exc :: Symbol
    type_resp :: Symbol

    SdofFrequencyProblem(sdof, freq, F; type_exc = :force, type_resp = :dis) = new(sdof, freq, F, type_exc, type_resp)
end

"""
    SdofFrequencySolution(u)

Structure containing the data of the solution of a frequency problem for a sdof system

# Field
* `u`: Solution of the frequency problem
   * Response spectrum (displacement, velocity, acceleration) [m, m/s, m/s²]
   * Or Frequency response function (FRF) (Admittance, Mobility, Accelerance) [m/N, m.s/N, m.s²/N]
"""
@with_kw struct SdofFrequencySolution
    u :: Vector{Complex{Float64}}
end

"""
    solve(prob::SdofFreeTimeProblem)

Compute the free response of a single degree of freedom (Sdof) system.

# Input
* `prob`: Structure containing the parameters of the Sdof problem

# Output
* `sol`: The response of the system at the given time points
"""
function solve(prob::SdofFreeTimeProblem)
    (; sdof, u0, t) = prob
    (; ω₀, ξ) = sdof
    x₀, v₀ = u0

    if ω₀ == 0.
        x = @. x₀ + v₀*t
        v = @. v₀*ones(length(t))
        a = @. zeros(length(t))
    elseif ξ < 1.
        Ω₀ = ω₀*√(1 - ξ^2)
        A = x₀
        B = (v₀ + ξ*ω₀*x₀)/Ω₀

        x = @. (A*cos(Ω₀*t) + B*sin(Ω₀*t))*exp(-ξ*ω₀*t)
        v = @. Ω₀*(-A*sin(Ω₀*t) + B*cos(Ω₀*t))*exp(-ξ*ω₀*t) - ξ*ω₀*x
        a = @. -2ξ*ω₀*v - ω₀^2*x

    elseif ξ == 1.
        A = x₀
        B = v₀ + ω₀*x₀

        x = @. (A + B*t)*exp(-ω₀*t)
        v = @. B*exp(-ω₀*t) - ω₀*x
        a = @. -2ω₀*v - ω₀^2*x
    else
        β = ω₀*√(ξ^2 - 1)
        A = x₀
        B = (v₀ + ξ*ω₀*x₀)/β

        x = @. exp(-ξ*ω₀*t)*(A*cosh(β*t) + B*sinh(β*t))
        v = @. β*(A*sinh(β*t) + B*cosh(β*t))*exp(-ξ*ω₀*t) - ξ*ω₀*x
        a = @. -2ξ*ω₀*v - ω₀^2*x
    end

    return SdofTimeSolution(x, v, a)
end

"""
    solve(prob::SdofHarmonicTimeProblem)

Computes the forced response of a single degree of freedom (Sdof) system due to an harmonic external force or base motion

# Inputs
* `prob`: Structure containing the parameters of the Sdof harmonic problem

# Output
* `sol`: The response of the system at the given time points
"""
function solve(prob::SdofHarmonicTimeProblem)

    (; sdof, u0, t, F, ω, type_exc) = prob
    (; m, ω₀, ξ) = sdof
    x₀, v₀ = u0

    if ξ == 0. && ω₀ == ω
        # Variation parameters
        if type_exc == :force
            α = 1/m
        else
            α = ω^2
        end

        ρ₀ = α*F/2ω
        A = x₀
        B = v₀/ω
        x = @. A*cos(ω*t) + B*sin(ω*t) + ρ₀*t*sin(ω*t)
        v = @. -A*ω*sin(ω*t) + B*ω*cos(ω*t) + ρ₀*(sin(ω*t) + ω*t*cos(ω*t))
        a = @. -A*ω^2*cos(ω*t) - B*ω^2*sin(ω*t) + ρ₀*(2ω*cos(ω*t) - ω^2*t*sin(ω*t))
    else
        if type_exc == :force
            X = F/m/(ω₀^2 - ω^2 + 2im*ξ*ω*ω₀)
        else
            X = F*(ω₀^2 + 2im*ξ*ω*ω₀)/(ω₀^2 - ω^2 + 2im*ξ*ω*ω₀)
        end

        ρ₀ = abs.(X)
        ϕ = angle.(X)

        A = x₀ - ρ₀*cos(ϕ)
        if ω₀ == 0.
            B = v₀ + ρ₀*ω*sin(ϕ)
            xh = @. A + B*t
            vh = @. B*ones(length(t))
            ah = @. zeros(length(t))
        elseif ξ < 1.
            Ω₀ = ω₀*√(1 - ξ^2)
            B = (v₀ + ξ*ω₀*A + ρ₀*ω*sin(ϕ))/Ω₀
            xh = @. (A*cos(Ω₀*t) + B*sin(Ω₀*t))*exp(-ξ*ω₀*t)
            vh = @. Ω₀*(-A*sin(Ω₀*t) + B*cos(Ω₀*t))*exp(-ξ*ω₀*t) - ξ*ω₀*xh
            ah = @. -2ξ*ω₀*vh - ω₀^2*xh
        elseif ξ == 1.
            B = v₀ + ω₀*A + ρ₀*ω*sin(ϕ)
            xh = @. (A + B*t)*exp(-ω₀*t)
            vh = @. B*exp(-ω₀*t) - ω₀*xh
            ah = @. -2ω₀*vh - ω₀^2*xh
        else
            β = ω₀*√(ξ^2 - 1)
            B = (v₀ + ξ*ω₀*A + ρ₀*ω*sin(ϕ))/β
            xh = @. (A*cosh(β*t) + B*sinh(β*t))*exp(-ξ*ω₀*t)
            vh = @. β*(A*sinh(β*t) + B*cosh(β*t))*exp(-ξ*ω₀*t) - ξ*ω₀*xh
            ah = @. -2ξ*ω₀*vh - ω₀^2*xh
        end

        x = @. xh + ρ₀*cos(ω*t + ϕ)
        v = @. vh - ρ₀*ω*sin(ω*t + ϕ)
        a = @. ah - ρ₀*ω^2*cos(ω*t + ϕ)
    end

    return SdofTimeSolution(x, v, a)
end

"""
    solve(prob::SdofForcedTimeProblem)

Computes the forced response of a single degree of freedom (Sdof) system due to an arbitrary external force or base motion

# Inputs
* `prob`: Structure containing the parameters of the Sdof forced problem
* `method`: Method to compute the Duhamel's integral
    * :interp: Interpolation + Gaussian quadrature (default)
    * :conv: Convolution

# Output
* `sol`: The response of the system at the given time points
"""
function solve(prob::SdofForcedTimeProblem; method = :interp)

    (; sdof, u0, t, F, type_exc) = prob
    (; m, ω₀, ξ) = sdof
    x₀, v₀ = u0

    # Time step
    Δt = t[2] - t[1]

    # Impulse response
    A = x₀
    if ω₀ == 0.
        B = v₀
        xh = @. A + B*t
        h = @. t/m
    elseif ξ < 1.
        Ω₀ = ω₀*√(1 - ξ^2)
        B = (v₀ + ξ*ω₀*x₀)/Ω₀
        xh = @. (A*cos(Ω₀*t) + B*sin(Ω₀*t))*exp(-ξ*ω₀*t)
        h = @. exp(-ξ*ω₀*t)*sin(Ω₀*t)/m/Ω₀
    elseif ξ == 1.
        B = v₀ + ω₀*x₀
        xh = @. (A + B*t)*exp(-ω₀*t)
        h = @. t*exp(-ω₀*t)/m
    else
        β = ω₀*√(ξ^2 - 1)
        B = (v₀ + ξ*ω₀*x₀)/β
        xh = @. (A*cosh(β*t) + B*sinh(β*t))*exp(-ξ*ω₀*t)
        h = @. exp(-ξ*ω₀*t)*sinh(β*t)/m/β
    end

    # Duhamel's integral
    if type_exc == :base
        k, c = ω₀^2*m, 2ξ*ω₀*m
        xb = F
        vb = gradient(xb, t)

        F = k*xb .+ c*vb
    end

    if method == :interp
        xp = duhamel_integral(F, h, t)
    else
        xp = Δt*conv(F, h)[1:length(F)]
    end

    x = xh .+ xp
    v = gradient(x, t)
    a = gradient(v, t)

    return SdofTimeSolution(x, v, a)
end

function duhamel_integral(F, h, t)
    # Interpolation
    fc = cubic_spline_interpolation(t, F)
    hc = cubic_spline_interpolation(t, h)

    # Gauss-Legendre quadrature
    nodes, weights = gausslegendre(500)
    n = length(nodes)
    nt = length(t)

    # Initialization
    x = undefs(nt)
    scaled_nodes = undefs(n)
    scaled_weights = undefs(n)
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

# Inputs
* `prob`: Structure containing the parameters of the Sdof FRF problem

# Output
* sol: Solution of the FRF problem
"""
function solve(prob::SdofFRFProblem)
    (; sdof, freq, type_exc, type_resp) = prob
    (; m, ω₀, ξ) = sdof
    ω = 2π*freq

    if type_exc == :force
        x = @. 1/m/(ω₀^2 - ω^2 + 2im*ξ*ω₀*ω)
    else
        x = @. (ω₀^2 + 2im*ξ*ω₀*ω)/(ω₀^2 - ω^2 + 2im*ξ*ω₀*ω)
    end

    if type_resp == :vel
        x .*= 1im*ω
    elseif type_resp == :acc
        @. x *= -ω^2
    end

    return SdofFrequencySolution(x)
end

"""
    solve(prob::SdofFrequencyProblem)

Compute the frequency response function of a single degree of freedom (Sdof) system

# Inputs
* `prob`: Structure containing the parameters of the Sdof frequency problem

# Output
* sol: Solution of the frequency problem
"""
function solve(prob::SdofFrequencyProblem)
    (; sdof, freq, F, type_exc, type_resp) = prob
    (; m, ω₀, ξ) = sdof
    ω = 2π*freq

    if type_exc == :force
        x = @. F/m/(ω₀^2 - ω^2 + 2im*ξ*ω₀*ω)
    else
        x = @. F*(ω₀^2 + 2im*ξ*ω₀*ω)/(ω₀^2 - ω^2 + 2im*ξ*ω₀*ω)
    end

    if type_resp == :vel
        x .*= 1im*ω
    elseif type_resp == :acc
        @. x *= -ω^2
    end

    return SdofFrequencySolution(x)
end

"""
    impulse_response(sdof, t)

Compute the impulse response of a single degree of freedom (Sdof) system

# Inputs
* `sdof`: Sdof structure
* `t`: Time points at which to evaluate the response

# Output
* `sol`: SdofImpulseSolution
"""
function impulse_response(sdof, t)
    prob = SdofFreeTimeProblem(sdof, [0., 1/sdof.m], t)

    return SdofImpulseSolution(solve(prob).u)
end

function srs(base_exc, freq, t, ξ = 0.05)
    nf = length(freq)

    srs = undefs(nf)
    x = undefs(nt)
    for (f, f₀) in freq
        ω₀ = 2π*f₀
        sdof = Sdof(1., f₀, ξ)
        prob = SdofForcedTimeProblem(sdof, [0., 0.], t, -base_exc)
        (; u, du) = solve(prob)

        x .= -2ξ*ω₀*du - ω₀^2*u
        srs[f] = maximum(abs, x)
    end

    return srs
end