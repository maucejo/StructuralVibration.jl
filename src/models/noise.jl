"""
    agwn(x, snr_dB; rst = true)

Adds a Gaussian White Noise (AGWN) to a signal `x` with a given SNR.

**Inputs**
* `x`: Signal (Real or Complex)
* `snr_dB`: signal to noise ratio [dB]
* `rst`: Reset the random number generator - Bool

**Output**
* `y`: Noisy signal
"""
function agwn(x, snr_dB; rst = true)

    # Reset the RNG if required
    if rst
        rng = MersenneTwister(1000)
    end

    SNR = 10^(snr_dB/10.)                     # SNR in linear scale
    En = vec(mean(abs2, x, dims = ndims(x)))  # Signal energy
    V = En/SNR                                # Noise variance

    σ = sqrt.(V)                              # Standard deviation
    n = σ.*randn(rng, eltype(x), size(x))     # Gaussian noise

    return x .+ n
end

"""
    acn(x, snr_dB, freq, color = :pink; rst = true)

Adds a complex Random Colored Noise (ACN) to a signal `x` with a given SNR

**Inputs**
* `x::VecOrMat{Real}`: Signal
* `snr_dB`: Signal to noise ratio [dB]
* `freq::AbstractVector`: Frequency range of interest
* `color`: Color of the noise
    * `:white`
    * `:pink` (default)
    * `:blue`
    * `:brown`
    * `:purple`
* `rst::Bool`: Reset the random number generator

**Output**
* `y`: Noisy signal
"""
function acn(x::VecOrMat{T}, snr_dB, freq::AbstractVector, color = :pink; rst = true) where {T <: Complex}

    # Reset the RNG if required
    if rst
        rng = MersenneTwister(1000)
    end

    ndx = ndims(x)
    SNR = 10^(snr_dB/10.)                     # SNR in linear scale
    En = vec(mean(abs2, x, dims = ndx))  # Signal energy
    V = En/SNR                                # Noise variance

    σ = sqrt.(V)                              # Standard deviation

    white_noise = randn(rng, eltype(x), size(x))
    scale = ones(eltype(x), size(x, ndx))
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
    colored_noise = white_noise.*scale'

    # Colored noise with scaled variance
    colored_noise .*= σ/std(colored_noise, dims = ndx)

    return x .+ colored_noise
end

"""
    acn(x, snr_dB, fs, color = :pink; band_freq = Float64[], rst = true)

Adds a complex Colored Noise (ACN) to a signal `x` with a given SNR

**Inputs**
* `x`: Signal
* `snr_dB`: Signal to noise ratio [dB]
* `fs`: Sampling frequency [Hz]
* `color`: Color of the noise
    * `:white`
    * `:pink` (default)
    * `:blue`
    * `:brown`
    * `:purple`
* `band_freq`: Frequencies used to defined the bandpass filter applied to the colored noise
* `rst`: Reset the random number generator

**Output**
* `y`: Noisy signal
"""
function acn(x::VecOrMat{T}, snr_dB, fs, color = :pink; band_freq = T[], rst = true) where {T <: Real}

    # Reset the RNG if required
    if rst
        rng = MersenneTwister(1000)
    end

    ndx = ndims(x)
    L = size(x, ndx)                        # Data dimensions
    SNR = 10^(snr_dB/10.)                   # SNR in linear scale
    En = vec(mean(abs2, x, dims = ndx))     # Signal energy
    V = En/SNR                              # Noise variance

    σ = sqrt.(V)                            # Standard deviation

    white_fft = rfft(randn(rng, size(x)), ndx)
    freq = rfftfreq(L, fs)

    scale = ones(eltype(x), length(freq))
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
    if ndx == 1
        colored_fft = white_fft.*scale
    else
        colored_fft = white_fft.*scale'
    end

    # Colored noise with scaled variance
    colored_noise = irfft(colored_fft, L, ndx)
    colored_noise .*= σ./std(colored_noise, dims = ndx)

    if length(band_freq) > 0
        flag = true
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
            flag = false
        end

        if flag
            df = DSP.digitalfilter(filter_type, DSP.Butterworth(4); fs)
            if ndx == 1
                colored_noise .= DSP.filtfilt(df, colored_noise)
            else
                N = size(x, 1)
                for i in 1:N
                    colored_noise[i, :] .= DSP.filtfilt(df, colored_noise[i, :])
                end
            end
        end
    end

    return x .+ colored_noise
end

"""
    mgwn(x, snr_dB; rst = true)

Adds a multiplicative Gaussian White Noise (MGWN) to a signal `x` with a given SNR

**Inputs**
* `x`: Signal
* `snr_dB`: Signal to noise ratio [dB]
* `rst`: Reset the random number generator

**Output**
* `y`: Noisy signal
"""
function mgwn(x, snr_dB; rst = true)

    # Reset the RNG if required
    if rst
        rng = MersenneTwister(1000)
    end

    SNR = 10^(snr_dB/10.)
    n = randn(rng, eltype(x), size(x))/sqrt(SNR)

    return @. (1. + n)*x
end

"""
    mixed_noise(x, snr_dB; rst = true)

Adds both additive and multiplicative Gaussian White Noise to a signal `x` with a given SNR

**Inputs**
* `x`: Signal
* `snr_dB`: Signal to noise ratio [dB]
* `rst`: Reset the random number generator

**Output**
* `y`: Noisy signal
"""
function mixed_noise(x, snr_dB; rst = true)

    # Reset the RNG if required
    if rst
        rng = MersenneTwister(1000)
    end

    SNR = 10^(snr_dB/10.)                      # SNR in linear scale
    En = mean(abs2, x, dims = ndims(x))[:]     # Signal energy
    V = En/SNR                                 # Noise variance

    σ = sqrt.(V)                               #Standard deviation
    addn = σ.*randn(rng, eltype(x), size(x))   # Gaussian noise

    muln = randn(rng, eltype(x), size(x))/sqrt(SNR)

    return @. (1. + muln)*x + addn
end