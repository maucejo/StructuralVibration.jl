abstract type ArbitraryExc end

"""
    Rectangle(F, tstart, duration)

Struct to define a rectangular excitation signal

**Fields**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
"""
@show_data struct Rectangle{Tf <: Real, Tt <: Real, Td <: Real} <: ArbitraryExc
    F::Tf
    tstart::Tt
    duration::Td
end

"""
    Triangle(F, tstart, duration)

Struct to define a triangular excitation signal

**Fields**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
"""
@show_data struct Triangle{Tf <: Real, Tt <: Real, Td <: Real} <: ArbitraryExc
    F::Tf
    tstart::Tt
    duration::Td
end

"""
    Hammer(F, tstart, k, θ)

Struct to define a hammer impact excitation signal

**Fields**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `p::Real`: Shape parameter
* `θ::Real`: Intensity parameter [s]
"""
@show_data struct Hammer{Tf <: Real, Tt <: Real, Tp <: Real} <: ArbitraryExc
    F::Tf
    tstart::Tt
    p::Tp
    θ::Tp
end

"""
    SmoothRect(F, tstart, tr, duration)

Struct to define a smooth rectangular excitation signal

**Fields**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
* `trise::Real`: Rise time from 0 to F [s]

*Note: `SmoothRect` is actually a custom Tukey window for which the coefficient α is computed to satisfy the `trise` given by the user*
"""
@show_data struct SmoothRect{Tf <: Real, Tt <: Real, Td <: Real, Tr <: Real} <: ArbitraryExc
    F::Tf
    tstart::Tt
    duration::Td
    trise::Tr
end

"""
    SineWave(F, tstart, duration, freq; zero_end = true)

Struct to define a sine wave excitation signal

**Constructor parameters**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
* `freq::Real`: Frequency of the excitation [Hz]
* `zero_end::Bool`: Boolean to set the excitation to 0 at the end of the duration (default = true)

**Fields**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
* `ω::Real`: Frequency of the excitation [Hz]
* `zero_end::Bool`: Boolean to set the excitation to 0 at the end of the duration (default = true)
"""
@show_data struct SineWave{Tf <: Real, Tt <: Real, Td <: Real, Tw <: Real} <: ArbitraryExc
    F::Tf
    tstart::Tt
    duration::Td
    ω::Tw
    ϕ::Tw
    zero_end::Bool

    SineWave(F::Tf, tstart::Tt, duration::Td, freq::Tw, ϕ::Tw = 0.; zero_end = true) where {Tf, Tt, Td, Tw} = new{Tf, Tt, Td, Tw}(F, tstart, duration, 2π*freq, ϕ, zero_end)
end

"""
    HalfSine(F, tstart, duration)

Struct to define a half sine excitation signal

**Fields**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
"""
@show_data struct HalfSine{Tf <: Real, Tt <: Real, Td <: Real} <: ArbitraryExc
    F::Tf
    tstart::Tt
    duration::Td
end

"""
    HaverSine(F, tstart, duration)

Struct to define a Haversine (or versed sine) excitation signal

**Fields**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
"""
@show_data struct HaverSine{Tf <: Real, Tt <: Real, Td <: Real} <: ArbitraryExc
    F::Tf
    tstart::Tt
    duration::Td
end

"""
    SweptSine(F, tstart, duration, fstart, fend, type)

Struct to define a swept sine excitation signal

**Fields**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
* `fstart::Real`: Starting frequency [Hz]
* `fend::Real`: Ending frequency [Hz]
* `type::Symbol`: Type of sweep
    * `:lin` - linear (default)
    * `:quad` - quadratic
    * `:log` - logarithmic
* `zero_end::Bool`: Boolean to set the excitation to 0 at the end of the duration (default = true)
"""
@show_data struct SweptSine{Tf <: Real, Tt <: Real, Td <: Real, Tr <: Real} <: ArbitraryExc
    F::Tf
    tstart::Tt
    duration::Td
    fstart::Tr
    fend::Tr
    type_swept::Symbol
    zero_end::Bool

    SweptSine(F::Tf, tstart::Tt, duration::Td, fstart::Tr, fend::Tr, type_swept::Symbol = :lin; zero_end::Bool = true) where {Tf, Tt, Td, Tr} = new{Tf, Tt, Td, Tr}(F, tstart, duration, fstart, fend, type_swept, zero_end)
end

"""
    GaussianPulse(F, tstart, duration, fc)

Struct to define a Gaussian pulse excitation signal

**Fields**
* `F::Real`: Amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
* `fc::Real`: Center frequency of the pulse [Hz]
* `precision::Real`: Precision of the pulse (default = 4.)

**Note**

The precision parameter calibrates the standard deviation of the pulse, so that the duration = n x σ, within some precision. If n = 1.96 then the confidence interval of the Gaussian distribution is 95%. To do so, we compute n so that at t = duration = n x σ , the amplitude of F x exp(-0.5*(t - duration/2)^2/sigma^2) = 10^(-precision)
"""
@show_data struct GaussianPulse{Tf <: Real, Tt <: Real, Td <: Real, Tc <: Real, Tp <: Real} <: ArbitraryExc
    F::Tf
    tstart::Tt
    duration::Td
    fc::Tc
    precision::Tp

    GaussianPulse(F::Tf, tstart::Tt, duration::Td, fc::Tc; precision::Tp = 4.) where {Tf, Tt, Td, Tc, Tp} = new{Tf, Tt, Td, Tc, Tp}(F, tstart, duration, fc, precision)
end

"""
    ColoredNoise(F, color)

Struct to define a colored noise excitation signal

**Fields**
* `F::Real`: Mean amplitude of the force [N]
* `tstart::Real`: Starting time of the excitation [s]
* `duration::Real`: Duration of the excitation [s]
* `σ::Real`: Target standard deviation of the colored noise
* `color::Symbol`: Color of the noise
    * `:white` (default)
    * `:pink`
    * `:blue`
    * `:brown`
    * `:purple`
* `band_freq::AbstractVector`: Frequencies used to defined the bandpass filter applied to the colored noise
"""
@show_data struct ColoredNoise{Tf <: Real, Tt <: Real, Td <: Real, Ts <: Real, Tb <: AbstractVector} <: ArbitraryExc
    F::Tf
    tstart::Tt
    duration::Td
    σ::Ts
    color::Symbol
    band_freq::Tb

    ColoredNoise(F::Tf, tstart::Tt, duration::Td, σ::Ts = 1.; color::Symbol = :white, band_freq::Tb = Tf[]) where {Tf, Tt, Td, Ts, Tb} = new{Tf, Tt, Td, Ts, Tb}(F, tstart, duration, σ, color, band_freq)
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

    (; F, tstart, duration) = type
    Ft = zeros(typeof(F), length(t))

    tb = tstart
    te = tstart + duration

    Ft[@. tb ≤ t ≤ te] .= F

    return Ft
end

# Triangle excitation
function excitation(type::Triangle, t)

    (; F, tstart, duration) = type

    Ft = zeros(typeof(F), length(t))

    # trise = (tbegin + tend)/2
    trise = (2tstart + duration)/2.

    pos_start = argmin(@. (t .- tstart)^2.)
    pos_middle = argmin(@. (t - trise)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)
    amp = 2F/duration

    @. Ft[pos_start:pos_middle] = amp*(t[pos_start:pos_middle] - tstart)
    @. Ft[pos_middle + 1:pos_end] = F - amp*(t[pos_middle + 1:pos_end] - trise)

    return Ft
end

# Hammer excitation
function excitation(type::Hammer, t)

    (; F, tstart, p, θ) = type

    Ft = zeros(typeof(F), length(t))

    # Check the type of t
    !isa(t, Array) ? t = collect(t) : nothing

    pos_start = argmin(@. (t - tstart)^2.)

    t_hammer = @. t[pos_start:end] - tstart

    @. Ft[pos_start:end] = F*(t_hammer/p/θ)^p*exp(-t_hammer/θ + p)

    return Ft
end

# Smooth rectangular excitation
function excitation(type::SmoothRect, t)

    (; F, tstart, trise, duration) = type

    Ft = zeros(typeof(F), length(t))

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)

    # Check duration
    Trect = duration - 2trise
    Trect < 0. ? error("duration must >= 2trise") : nothing

    pos_rect_start = argmin(@. (t - tstart - trise)^2.)
    pos_rect_end = argmin(@. (t - tstart - duration + trise).^2.)

    α = 2trise/duration

    t_rise = t[pos_start:pos_rect_start] .- tstart
    Frise = @. 0.5*F*(1. - cos(2π*t_rise/α/duration))

    t_desc = @. t[pos_rect_end:pos_end] - tstart - duration
    Fdesc = @. 0.5*F*(1. - cos(2π*t_desc/α/duration))

    Ft[pos_start:pos_rect_start] .= Frise
    Ft[pos_rect_start:pos_rect_end] .= F
    Ft[pos_rect_end:pos_end] .= Fdesc

    return Ft
end

# Sine wave excitation
function excitation(type::SineWave, t)

    (; F, tstart, duration, ω, ϕ, zero_end) = type

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

    @. Ft[pos_exc_t] = F*sin(ω*(t[pos_exc_t] - tstart) + ϕ)

    return Ft
end

# Half sine excitation
function excitation(type::HalfSine, t)

    (; F, tstart, duration) = type

    Ft = zeros(typeof(F), length(t))

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)
    pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])

    @. Ft[pos_exc_t] = F*sin(π*(t[pos_exc_t] - tstart)/duration)

    return Ft
end

# Haversine excitation
function excitation(type::HaverSine, t)

    (; F, tstart, duration) = type

    Ft = zeros(typeof(F), length(t))

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration)^2.)
    pos_exc_t = findall(@. t[pos_start] ≤ t ≤ t[pos_end])

    @. Ft[pos_exc_t] = F*(1. - cos(2π*(t[pos_exc_t] - tstart)/duration))/2

    return Ft
end

# Swept sine excitation
function excitation(type::SweptSine, t)

    (; F, tstart, duration, fstart, fend, type_swept, zero_end) = type

    Ft = zeros(typeof(F), length(t))

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

    @. Ft[pos_exc_t] = F*sin(ϕ)

    return Ft
end

# Gaussian pulse excitation
function excitation(type::GaussianPulse, t)
    (; F, tstart, duration, fc, precision) = type

    Ft = zeros(typeof(F), length(t))

    pos_start = argmin((t .- tstart).^2.)
    tpulse = t[pos_start:end] .- tstart .- duration/2

    # The precision parameter calibrates the standard deviation of the pulse, so that the duration = n*σ, within some precision. If n = 1.96 then the confidence interval of the Gaussian distribution is 95%. To do so, we compute n so that at t = duration = n*σ , the amplitude of F*exp(-0.5*(t - duration/2)^2/sigma^2) = 10^(-precision)
    n = 2sqrt(2(precision*log(10) + log(F)))
    σ = duration/n

    @. Ft[pos_start:end] = F*exp(-0.5*(tpulse/σ)^2)*cos(2π*fc*tpulse)

    return Ft
end

# Colored noise excitation
function excitation(type::ColoredNoise, t)

    (; F, tstart, duration, σ, color, band_freq) = type
    N = length(t)
    Ft = zeros(N)

    # Sampling frequency
    fs = 1/(t[2] - t[1])

    # Generate the fft of a white noise
    if color == :white
        Ft .= F .+ σ*randn(N)
    else
        white_fft = rfft(randn(N))
        freq = rfftfreq(N, fs)

        scale = zeros(length(freq))
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

        @. Ft = F + colored_noise
    end

    if length(band_freq) > 0
        if band_freq[1] > freq[1] && band_freq[2] < freq[end]
            # Band-pass filter
            filter_type = DSP.Bandpass(band_freq[1], band_freq[2])
        elseif band_freq[1] < freq[1] && band_freq[2] < freq[end]
            # Low-pass filter
            filter_type = DSP.Lowpass(band_freq[2])
        elseif band_freq[1] > freq[1] && band_freq[2] > freq[end]
            # High-pass filter
            filter_type = DSP.Highpass(band_freq[1])
        else
            return x .+ colored_noise
        end

        df = DSP.digitalfilter(filter_type, DSP.Butterworth(4); fs)
        Ft .= DSP.filtfilt(df, Ft)
    end

    pos_start = argmin(@. (t - tstart)^2.)
    pos_end = argmin(@. (t - tstart - duration).^2.)

    Ft[1:pos_start] .= 0.
    Ft[pos_end:end] .= 0.

    return Ft
end