abstract type ArbitraryExc end

"""
    Rectangle(F₀, tstart, duration)

Struct to define a rectangular excitation signal

**Fields**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
"""
@with_kw struct Rectangle <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    duration::Float64
end

"""
    Triangle(F₀, tstart, duration)

Struct to define a triangular excitation signal

**Fields**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
"""
@with_kw struct Triangle <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    duration::Float64
end

"""
    Hammer(F₀, tstart, k, θ)

Struct to define a hammer impact excitation signal

**Fields**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `k`: Shape parameter
* `θ`: Intensity parameter [s]
"""
@with_kw struct Hammer <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    k::Float64
    θ::Float64
end

"""
    SmoothRect(F₀, tstart, tr, duration)

Struct to define a smooth rectangular excitation signal

**Fields**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
* `trise`: Rise time from 0 to F₀ [s]
"""
@with_kw struct SmoothRect <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    duration::Float64
    trise::Float64
end

"""
    SineWave(F₀, tstart, duration, freq; zero_end = true)

Struct to define a sine wave excitation signal

**Constructor parameters**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
* `freq`: Frequency of the excitation [Hz]
* `zero_end`: Boolean to set the excitation to 0 at the end of the duration (default = true)

**Fields**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
* `ω`: Frequency of the excitation [Hz]
* `zero_end`: Boolean to set the excitation to 0 at the end of the duration (default = true)
"""
@with_kw struct SineWave <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    duration::Float64
    ω::Float64
    ϕ::Float64
    zero_end::Bool

    SineWave(F₀, tstart, duration, freq, ϕ = 0.; zero_end = true) = new(F₀, tstart, duration, 2π*freq, ϕ, zero_end)
end

"""
    HalfSine(F₀, tstart, duration)

Struct to define a half sine excitation signal

**Fields**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
"""
@with_kw struct HalfSine <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    duration::Float64
end

"""
    HaverSine(F₀, tstart, duration)

Struct to define a Haversine (or versed sine) excitation signal

**Fields**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
"""
@with_kw struct HaverSine <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    duration::Float64
end

"""
    SweptSine(F₀, tstart, duration, fstart, fend, type)

Struct to define a swept sine excitation signal

**Fields**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
* `fstart`: Starting frequency [Hz]
* `fend`: Ending frequency [Hz]
* `type`: Type of sweep
    * `:lin` - linear (default)
    * `:quad` - quadratic
    * `:log` - logarithmic
* `zero_end`: Boolean to set the excitation to 0 at the end of the duration (default = true)
"""
@with_kw struct SweptSine <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    duration::Float64
    fstart::Float64
    fend::Float64
    type_swept::Symbol
    zero_end::Bool

    SweptSine(F₀, tstart, duration, fstart, fend, type_swept = :lin; zero_end = true) = new(F₀, tstart, duration, fstart, fend, type_swept, zero_end)
end

"""
    GaussianPulse(F₀, tstart, duration, fc)

Struct to define a Gaussian pulse excitation signal

**Fields**
* `F₀`: Amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
* `fc`: Center frequency of the pulse [Hz]
* `precision`: Precision of the pulse (default = 4)

*Note on the parameter `precision`*

The precision parameter calibrates the standard deviation of the pulse, so that the duration = n x σ, within some precision. If n = 1.96 then the confidence interval of the Gaussian distribution is 95%. To do so, we compute n so that at t = duration = n x σ , the amplitude of F₀*exp(-0.5*(t - duration/2)^2/sigma^2) = 10^(-precision)
"""
@with_kw struct GaussianPulse <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    duration::Float64
    fc::Float64
    precision::Float64

    GaussianPulse(F₀, tstart, duration, fc; precision = 4) = new(F₀, tstart, duration, fc, precision)
end

"""
    ColoredNoise(F₀, color)

Struct to define a colored noise excitation signal

**Fields**
* `F₀`: Mean amplitude of the force [N]
* `tstart`: Starting time of the excitation [s]
* `duration`: Duration of the excitation [s]
* `σ`: Target standard deviation of the colored noise
* `color`: Color of the noise
    * `:white` (default)
    * `:pink`
    * `:blue`
    * `:brown`
    * `:purple`
* `band_freq`: Frequencies used to defined the bandpass filter applied to the colored noise
"""
@with_kw struct ColoredNoise <: ArbitraryExc
    F₀::Float64
    tstart::Float64
    duration::Float64
    σ::Float64
    color::Symbol
    band_freq::Vector{Float64}

    ColoredNoise(F₀, tstart, duration, σ = 1.; color = :white, band_freq = Float64[]) = new(F₀, tstart, duration, σ, color, band_freq)
end

"""
    excitation(type, t)

Computes different types of excitation signals

**Inputs**
* type : Excitation type
    1. `Triangle`
    2. `Rectangle`
    3. `Hammer`
    4. `SmoothRect`
    5. `SineWave`
    6. `HalfSine`
    7. `HaverSine`
    8. `SweptSine`
    9. `GaussianPulse`
    10. `ColoredNoise`
* `t`: Time vector

**Output**
* `F`: Excitation signal
"""
function excitation(type::Rectangle, t)

    (; F₀, tstart, duration) = type
    Ft = zeros(length(t))

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)

    pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])

    Ft[pos_exc_t] .= F₀

    return Ft
end

# Triangle excitation
function excitation(type::Triangle, t)

    (; F₀, tstart, duration) = type

    Ft = zeros(length(t))

    trise = (2tstart + duration)/2.
    pos_start = argmin(@. (t .- tstart)^2.)
    pos_middle = argmin(@. (t - trise)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)
    amp = 2F₀/duration

    @. Ft[pos_start:pos_middle] = amp*(t[pos_start:pos_middle] - type.tstart)
    @. Ft[pos_middle + 1:pos_end] = F₀ - amp*(t[pos_middle + 1:pos_end] - trise)

    return Ft
end

# Hammer excitation
function excitation(type::Hammer, t)

    (; F₀, tstart, k, θ) = type

    Ft = zeros(length(t))

    # Check the type of t
    !isa(t, Array) ? t = collect(t) : nothing

    pos_start = argmin(@. (t - tstart)^2.)

    t_hammer = @. t[pos_start:end] - tstart

    @. Ft[pos_start:end] = F₀*t_hammer^(k - 1.)*exp(-t_hammer/θ)/((k - 1.)*θ)^(k - 1.)/exp(1. - k)

    return Ft
end

# Smooth rectangular excitation
function excitation(type::SmoothRect, t)

    (; F₀, tstart, trise, duration) = type

    Ft = zeros(length(t))

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)

    # Check duration
    Trect = duration - 2trise
    Trect < 0. ? error("duration must >= 2trise") : nothing

    pos_rect_start = argmin(@. (t - tstart - trise)^2.)
    pos_rect_end = argmin(@. (t - tstart - duration + trise).^2.)

    α = 2trise/duration

    t_rise = t[pos_start:pos_rect_start] .- tstart
    Frise = @. 0.5*F₀*(1. - cos(2π*t_rise/α/duration))

    t_desc = @. t[pos_rect_end:pos_end] - tstart - duration
    Fdesc = @. 0.5*F₀*(1. - cos(2π*t_desc/α/duration))

    Ft[pos_start:pos_rect_start] .= Frise
    Ft[pos_rect_start:pos_rect_end] .= F₀
    Ft[pos_rect_end:pos_end] .= Fdesc

    return Ft
end

# Sine wave excitation
function excitation(type::SineWave, t)

    (; F₀, tstart, duration, ω, ϕ, zero_end) = type

    nt = length(t)
    Ft = zeros(nt)

    pos_start = argmin(@. (t - tstart)^2.)

    tsw = tstart + duration

    if tsw ≥ t[end]
        pos_end = nt
    else
        if zero_end
            n = round(Int, (ω*duration + ϕ)/2π)
            duration = (2π*n - ϕ)/ω
            pos_end = argmin(@. (t - tstart - duration)^2.)

            # Check duration
            tstart + duration ≥ t[end] ? Warn("The duration of the sine wave is too long to performed zero-end operation.") : nothing
        else
            pos_end = argmin(@. (t - tsw)^2.)
        end
    end

    pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])

    @. Ft[pos_exc_t] = F₀*sin(ω*(t[pos_exc_t] - tstart) + ϕ)

    return Ft
end

# Half sine excitation
function excitation(type::HalfSine, t)

    (; F₀, tstart, duration) = type

    Ft = zeros(length(t))

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)
    pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])

    @. Ft[pos_exc_t] = F₀*sin(π*(t[pos_exc_t] - tstart)/duration)

    return Ft
end

# Haversine excitation
function excitation(type::HaverSine, t)

    (; F₀, tstart, duration) = type

    Ft = zeros(length(t))

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)
    pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])

    @. Ft[pos_exc_t] = F₀*(1. - cos(2π*(t[pos_exc_t] - tstart)/duration))/2

    return Ft
end

# Swept sine excitation
function excitation(type::SweptSine, t)

    (; F₀, tstart, duration, fstart, fend, type_swept, zero_end) = type

    Ft = zeros(length(t))

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)
    pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])

    if fstart == fend
        ϕ = @. 2π*fstart*(t[pos_exc_t] - tstart)
    else
        if type_swept == :lin
            if zero_end
                n = round((fstart + fend)*duration/2)
                duration = 2n/(fstart + fend)

                pos_end = argmin(@. (t - tstart - duration)^2.)
                pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])
            end

            β = (fend - fstart)/duration
            ϕ = @. 2π*(fstart*(t[pos_exc_t] - tstart) + 0.5*β*(t[pos_exc_t] - tstart)^2)

        elseif type_swept == :quad
            if zero_end
                n = round((2fstart + fend)*duration/3)
                duration = 3n/(2fstart + fend)

                pos_end = argmin(@. (t - tstart - duration)^2.)
                pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])
            end

            β = (fend - fstart)/duration^2
            ϕ = @. 2π*(fstart*(t[pos_exc_t] - tstart) + β*(t[pos_exc_t] - tstart)^3/3)

        elseif type_swept == :log
            # Check frequency condition
            fstart*fend ≤ 0. ? error("fstart and fend must have the same sign and be different from 0.") : nothing

            if zero_end
                n = round((fend - fstart)*duration/log(fend/fstart))
                duration = n*log(fend/fstart)/(fend - fstart)

                # Check duration
                tstart + duration ≥ t[end] ? warning("The duration of the swept sine is too long to performed zero-end operation.") : nothing

                pos_end = argmin(@. (t - tstart - duration)^2.)
                pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])
            end

            β = (fend/fstart)^(1/duration)

            ϕ = @. 2π*fstart*(β^(t[pos_exc_t] - tstart) - 1)/log(β)
        else
            error("The type of swept sine is not defined.")
        end
    end

    @. Ft[pos_exc_t] = F₀*sin(ϕ)

    return Ft
end

# Gaussian pulse excitation
function excitation(type::GaussianPulse, t)
    (; F₀, tstart, duration, fc, precision) = type

    Ft = zeros(length(t))

    pos_start = argmin((t .- tstart).^2.)
    tpulse = t[pos_start:end] .- tstart .- duration/2

    # The precision parameter calibrates the standard deviation of the pulse, so that the duration = n*σ, within some precision. If n = 1.96 then the confidence interval of the Gaussian distribution is 95%. To do so, we compute n so that at t = duration = n*σ , the amplitude of F₀*exp(-0.5*(t - duration/2)^2/sigma^2) = 10^(-precision)
    n = 2sqrt(precision*2log(10) + 2log(F₀))
    σ = duration/n

    @. Ft[pos_start:end] = F₀*exp(-0.5*(tpulse/σ)^2)*cos(2π*fc*tpulse)

    return Ft
end

# Colored noise excitation
function excitation(type::ColoredNoise, t)

    (; F₀, tstart, duration, σ, color, band_freq) = type
    N = length(t)
    Ft = zeros(N)

    # Sampling frequency
    fs = 1/(t[2] - t[1])

    # Generate the fft of a white noise
    if color == :white
        Ft .= F₀ .+ σ*randn(N)
    else
        white_fft = rfft(randn(N))
        freq = rfftfreq(N, fs)

        scale = ones(length(freq))
        if color == :pink
            @. scale[2:end] = 1/sqrt(freq[2:end])
        elseif color == :blue
            @. scale = sqrt(freq)
        elseif color == :brown
            @. scale[2:end] = 1/freq[2:end]
        elseif color == :purple
            scale .= freq
        end

        # Energy preservation of the white noise
        scale ./= sqrt(mean(scale.^2))

        colored_noise_fft = white_fft.*scale

        # Colored noise with scaled variance
        colored_noise =  irfft(colored_noise_fft, N)
        colored_noise .*= σ/std(colored_noise)

        @. Ft = F₀ + colored_noise
    end

    if length(band_freq) > 0
        if band_freq[1] > freq[1] && band_freq[2] < freq[end]
            # Band-pass filter
            filter_type = Bandpass(band_freq[1], band_freq[2])
        elseif band_freq[1] < freq[1] && band_freq[2] < freq[end]
            # Low-pass filter
            filter_type = Lowpass(band_freq[2])
        elseif band_freq[1] > freq[1] && band_freq[2] > freq[end]
            # High-pass filter
            filter_type = Highpass(band_freq[1])
        else
            return x .+ colored_noise
        end

        df = digitalfilter(filter_type, Butterworth(4); fs)
        Ft .= filtfilt(df, Ft)
    end

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration).^2.)

    Ft[1:pos_start] .= 0.
    Ft[pos_end:end] .= 0.

    return Ft
end