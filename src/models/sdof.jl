"""
Structure containing the data of a sdof system

# Parameters
* m: Mass [kg]
* ω₀: natural angular frequency [rad/s]
* ξ: Damping ratio
"""
@with_kw struct SDOF
    m :: Float64
    ω₀ ::Float64
    ξ :: Float64
end

"""
    free_response(s::SDOF, t, x₀ = 0., v₀ = 1.)

Compute the free response of a single degree of freedom (SDOF) system.

# Inputs
- s: Structure containing the parameters of the SDOF system
- t: A vector of time points at which to evaluate the response
- x₀: Initial displacement (default is 0.0)
- v₀: Initial velocity (default is 1.0)

# Outputs
- rep: The response of the system at the given time points
- env: The envelope of the response
"""
function free_response(s::SDOF, t, x₀ = 0., v₀ = 1.)
    (; ω₀, ξ) = s
    nt = length(t)

    if ξ < 1.
        Ω₀ = ω₀*√(1 - ξ^2)
        A = x₀
        B = (v₀ + ξ*ω₀*x₀)/Ω₀

        rep = @. exp(-ξ*ω₀*t)*(A*cos(Ω₀*t) + B*sin(Ω₀*t))
        env = @. exp(-ξ*ω₀*t)*√(A^2 + B^2)

    elseif ξ == 1.
        A = x₀
        B = v₀ + ω₀*x₀

        rep = @. (A + B*t)*exp(-ω₀*t)
        env = zeros(nt)
    else
        β = ω₀*√(ξ^2 - 1)
        A = x₀
        B = (v₀ + ξ*ω₀*x₀)/β

        rep = @. exp(-ξ*ω₀*t)*(A*cosh(β*t) + B*sinh(β*t))
        env = zeros(nt)
    end

    return rep, env
end

"""
    forced_response(s::SDOF, Amp::Float64, ω::Float64, t, x₀ = 0., v₀ = 0., type= :force)

    forced_response(s::SDOF, F::Vector{Float64}, t, x₀ = 0., v₀ = 0.; type = :force)

Computes the forced response of a single degree of freedom (SDOF) system due to a harmonic external force or base motion

# Inputs
- s: Structure containing the parameters of the SDOF system
- t: A vector of time points at which to evaluate the response
- x₀: Initial displacement (default is 0.)
- v₀: Initial velocity (default is 0.)
- type: Type of excitation
    - :force: External force (default)
    - :base: Base motion

For a harmonic force:
- Amp: Amplitude of the force excitation [N]or base motion [m]
- ω: Frequency of the excitation [rad/s]

For any time dependence:
- F: Vector of the force excitation [N] or base motion [m]

# Outputs
- x: Response of the system at the given time points
- xh: Homogeneous solution
- xp: Particular solution
"""
function forced_response(s::SDOF, Amp::Float64, ω::Float64, t, x₀ = 0., v₀ = 0.; type = :force)
    (; m, ω₀, ξ) = s

    if type == :force
        A₀ = Amp/m/√((ω₀^2 - ω^2)^2 + (2ξ*ω*ω₀)^2)
        ϕ = atan(2ξ*ω*ω₀, ω₀^2 - ω^2)
    else
        A₀ = X₀*√(ω₀^4 + (2ξ*ω*ω₀)^2)/√((ω₀^2 - ω^2)^2 + (2ξ*ω*ω₀)^2)
        ϕ = atan(2ξ*ω*ω₀, ω₀^2 - ω^2) - atan(2ξ*ω, ω₀)
    end

    A = x₀ - A₀*cos(ϕ)
    if ξ < 1.
        Ω₀ = ω₀*√(1 - ξ^2)
        B = (v₀ + ξ*ω₀*A - A₀*ω*sin(ϕ))/Ω₀
        xh = @. (A*cos(Ω₀*t) + B*sin(Ω₀*t))*exp(-ξ*ω₀*t)
    elseif B == 1.
        B = v₀ + ω₀*A - A₀*ω*sin(ϕ)
        xh = @. (A + B*t)*exp(-ω₀*t)
    else
        β = ω₀*√(ξ^2 - 1)
        B = (v₀ + ξ*ω₀*A - A₀*ω*sin(ϕ))/β
        xh = @. (A*cosh(β*t) + B*sinh(β*t))*exp(-ξ*ω₀*t)
    end

    xp = A₀*cos.(ω*t .- ϕ)
    x = xh .+ xp

    return x, xh, xp
end

function forced_response(s::SDOF, F::Vector{Float64}, t, x₀ = 0., v₀ = 0.; type = :force)
    (; m, ω₀, ξ) = s
    # Time step
    Δt = t[2] - t[1]

    # Impulse response
    if ξ < 1.
        Ω₀ = ω₀*√(1 - ξ^2)
        h = @. exp(-ξ*ω₀*t)*sin(Ω₀*t)/m/Ω₀
    elseif ξ == 1.
        h = @. t*exp(-ω₀*t)/m
    else
        β = ω₀*√(ξ^2 - 1)
        h = @. exp(-ξ*ω₀*t)*sinh(β*t)/m/β
    end

    # Free response
    xh = free_response(s, t, x₀, v₀)[1]

    # Duhamel's integral
    if type == :base
        k, c = ω₀^2*m, 2ξ*ω₀*m
        xb = F
        vb = gradient(xb, t)

        F = k*xb .+ c*vb
    end

    xp = Δt*conv(F, h)[1:length(F)]

    x = xh .+ xp

    return x, xh, xp
end

"""
    frf(s::SDOF, freq; type = :force)

Compute the frequency response function of a single degree of freedom (SDOF) system

# Inputs
- s: Structure containing the parameters of the SDOF system
- freq: Vector of frequencies [Hz]
- type: Type of excitation
    - :force: Transfert function (default)
    - :base: Transmissibility function

# Output
- H: FRF of the system at the given frequencies
"""
function frf(s::SDOF, freq; type = :force)
    (; m, ω₀, ξ) = s
    ω = 2π*freq

    if type == :force
        return @. 1/m/(ω₀^2 - ω^2 + 2im*ξ*ω₀*ω)
    else
        return @. (ω₀^2 + 2im*ξ*ω₀*ω)/(ω₀^2 - ω^2 + 2im*ξ*ω₀*ω)
    end
end

"""
    freq_response(s::SDOF, F::Vector{Float64}, freq; type = :force)

Compute the frequency response of a single degree of freedom (SDOF) system due to an external force or a base motion

# Inputs
- s: Structure containing the parameters of the SDOF system
- F: Vector of the force [N] or base motion amplitude [m]
- freq: Vector of frequencies [Hz]
- type: Type of excitation
    - :force: External force (default)
    - :base: Base motion

# Output
- y: Response of the system at the given frequencies
"""
function freq_response(s::SDOF, F::Vector{Float64}, freq; type = :force)
    return frf(s, freq, type = type).*F
end