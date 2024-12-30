"""
    Rectangle(F₀, ts, duration)

Struct to define a rectangular excitation signal

# Parameters
* F₀ : Amplitude of the force [N]
* ts : Starting time of the excitation [s]
* duration : Duration of the excitation [s]
"""
@with_kw struct Rectangle
    F₀::Float64
    ts::Float64
    duration::Float64
end

"""
    Triangle(F₀, ts, duration)

Struct to define a triangular excitation signal

# Parameters
* F₀ : Amplitude of the force [N]
* ts : Starting time of the excitation [s]
* duration : Duration of the excitation [s]
"""
@with_kw struct Triangle
    F₀::Float64
    ts::Float64
    duration::Float64
end

"""
    RandomExc(F₀, ts, duration, σ)

Struct to define a random excitation signal

# Parameters
* F₀ : Amplitude of the force [N]
* ts : Starting time of the excitation [s]
* duration : Duration of the excitation [s]
* σ : Standard deviation of the random excitation
"""
@with_kw struct RandomExc
    F₀::Float64
    ts::Float64
    duration::Float64
    σ::Float64
end

"""
    Hammer(F₀, ts, k, θ)

Struct to define a hammer impact excitation signal

# Parameters
* F₀ : Amplitude of the force [N]
* ts : Starting time of the excitation [s]
* k : Shape parameter
* θ : Intensity parameter [s]
"""
@with_kw struct Hammer
    F₀::Float64
    ts::Float64
    k::Float64
    θ::Float64
end

"""
    SmoothRect(F₀, ts, tr, duration)

Struct to define a smooth rectangular excitation signal

# Parameters
* F₀ : Amplitude of the force [N]
* ts : Starting time of the excitation [s]
* tr : Rise time from 0 to F₀ [s]
* duration : Duration of the excitation [s]
"""
@with_kw struct SmoothRect
    F₀::Float64
    ts::Float64
    tr::Float64
    duration::Float64
end

"""
    excitation(type, t)

Computes different types of excitation signals

# Parameters
* type : Struct of excitation type
    1. Triangle
    2. Rectangle
    3. Hammer
    4. SmoothRect
    5. Random

# Output
* F : Vector of excitation evolution over time [N]
"""
# Rectangle excitation
function excitation(type::Rectangle, t)

    (; F₀, ts, duration) = type
    Ft = zeros(length(t))

    pos_start = argmin((t .- ts).^2.)
    pos_end = argmin((t .- ts .- duration).^2.)

    pos_exc_t = findall(t[pos_start] .≤ t .≤ t[pos_end])

    Ft[pos_exc_t] .= F₀

    return Ft
end

# Triangle excitation
function excitation(type::Triangle, t)

    (; F₀, ts, duration) = type

    Ft = zeros(length(t))

    tr = (2ts + duration)/2.
    pos_start = argmin((t .- ts).^2.)
    pos_middle = argmin((t .- tr).^2.)
    pos_end = argmin((t .- ts .- duration).^2.)
    amp = 2F₀/duration

    Ft[pos_start:pos_middle] = amp*(t[pos_start:pos_middle] .- type.ts)
    Ft[pos_middle + 1:pos_end] = F₀ .- amp*(t[pos_middle + 1:pos_end] .- tr)

    return Ft
end

# Random excitation
function excitation(type::RandomExc, t)

    (; F₀, ts, duration, σ) = type

    Ft = zeros(length(t))

    pos_start = argmin((t .- ts).^2.)
    pos_end = argmin((t .- ts .- duration).^2.)

    pos_exc_t = findall(t[pos_start] .≤ t .≤ t[pos_end])

    Ft[pos_exc_t] .= F₀ .+ σ*randn(length(pos_exc_t))

    return Ft
end

# Hammer excitation
function excitation(type::Hammer, t)

    (; F₀, ts, k, θ) = type

    Ft = zeros(length(t))

    if !isa(t, Array)
        t = collect(t);
    end

    pos_start = argmin((t .- ts).^2.)

    t_hammer = t[pos_start:end] .- ts

    Ft[pos_start:end] = F₀*t_hammer.^(k - 1.).*exp.(-t_hammer./θ)/((k - 1.)*θ)^(k - 1.)/exp(1. - k)

    return Ft
end

# Smooth rectangular excitation
function excitation(type::SmoothRect, t)

    (; F₀, ts, tr, duration) = type

    Ft = zeros(length(t))

    pos_start = argmin((t .- ts).^2.)
    pos_end = argmin((t .- ts .- duration).^2.)

    Trect = duration - 2tr
    if isless(Trect, 0.)
        error("Il faut que duration >= 2tr")
    end

    pos_rect_start = argmin((t .- ts .- tr).^2.)
    pos_rect_end = argmin((t .- ts .- duration .+ tr).^2.)

    α = 2tr/duration

    t_rise = t[pos_start:pos_rect_start] .- ts
    Frise = 0.5*F₀*(1. .- cos.(2π*t_rise/α/duration))

    t_desc = t[pos_rect_end:pos_end] .- ts .- duration
    Fdesc = 0.5*F₀*(1. .- cos.(2π*t_desc/α/duration))

    Ft[pos_start:pos_rect_start] = Frise
    Ft[pos_rect_start:pos_rect_end] .= F₀
    Ft[pos_rect_end:pos_end] = Fdesc

    return Ft
end