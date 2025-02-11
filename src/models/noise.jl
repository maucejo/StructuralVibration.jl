"""
    agwn(x, snr_dB; rst = true)

Adds a Gaussian White Noise (AWGN) to a signal `x` with a given SNR.

# Inputs
* `x`: signal
* `snr_dB`: signal to noise ratio [dB]
* `rst`: reset the random number generator - Bool

# Output
* `y`: noisy signal - Matrix{ComplexF64}
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

# Inputs
* `x`: Signal
* `snr_dB`: Signal to noise ratio [dB]
* `freq`: Frequency range of interest
* `color`: Color of the noise
    * `:white`
    * `:pink` (default)
    * `:blue`
    * `:brown`
    * `:purple`
* `band_freq`: Frequencies used to defined the bandpass filter applied to the colored noise
* `rst`: reset the random number generator

# Output
* `y`: noisy signal
"""
function acn(x::VecOrMat{Complex{Float64}}, snr_dB, freq::AbstractArray, color = :pink; rst = true)

    # Reset the RNG if required
    if rst
        rng = MersenneTwister(1000)
    end

    SNR = 10^(snr_dB/10.)                     # SNR in linear scale
    En = vec(mean(abs2, x, dims = ndims(x)))  # Signal energy
    V = En/SNR                                # Noise variance

    σ = sqrt.(V)                              # Standard deviation

    white_noise = randn(rng, eltype(x), size(x))
    scale = ones(size(x, ndims(x)))
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
    colored_noise = white_noise.*scale

    # Colored noise with scaled variance
    colored_noise .*= σ/std(colored_noise)

    return x .+ colored_noise
end

"""
    acn(x, snr_dB, color = :pink; band_freq = Float64[], rst = true)

Adds a complex Colored Noise (ACN) to a signal `x` with a given SNR

# Inputs
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
* `rst`: reset the random number generator

# Output
* `y`: noisy signal
"""
function acn(x::VecOrMat{Float64}, snr_dB, fs::Float64, color = :pink; band_freq = Float64[], rst = true)

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
    if ndx == 1
        colored_fft = white_fft.*scale
    else
        colored_fft = white_fft.*scale'
    end

    # Colored noise with scaled variance
    colored_noise = irfft(colored_fft, L, ndx)
    colored_noise .*= σ/std(colored_noise)

    if length(band_freq) > 0
        flag = true
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
            flag = false
        end

        if flag
            df = digitalfilter(filter_type, Butterworth(4); fs)
            if ndx == 1
                colored_noise .= filtfilt(df, colored_noise)
            else
                N = size(x, 1)
                for i in 1:N
                    colored_noise[i, :] .= filtfilt(df, colored_noise[i, :])
                end
            end
        end
    end

    return x .+ colored_noise
end

"""
    mult_noise(x, snr_dB; rst = true)

Adds a multiplicative Gaussian White Noise (AWGN) to a signal `x` with a given SNR

# Inputs
* `x`: Signal
* `snr_dB`: Signal to noise ratio [dB]
* `rst`: reset the random number generator

# Output
* `y`: noisy signal
"""
function mult_noise(x, snr_dB; rst = true)

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

# Inputs
* `x`: Signal
* `snr_dB`: Signal to noise ratio [dB]
* `rst`: reset the random number generator

# Output
* `y`: noisy signal
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