"""
    FFTParameters

Structure to store the time and frequency parameters for the FFT analysis

# Constructor parameters
* `fs`: Sampling rate of the signal
* `bs`: Block size of the signal
* `pow2`: Flag for finding the next power of 2

# Inputs
* `t`: Time vector
* `dt`: Time step
* `freq`: Frequency vector
* `df`: Frequency resolution
"""
@with_kw struct FFTParameters
    t
    dt::Float64
    tspan
    freq
    df::Float64
    freq_span
    fs::Int
    bs::Int

    function FFTParameters(fs::Int, bs::Int; pow2 = true)
        # Sample rate & block size
        if pow2
            fs = nextpow(2, fs)
            bs = nextpow(2, bs)
        end

        useful_bandwidth = fs/2.56

        # Time vector
        dt = 1/fs
        tmax = bs*dt
        t = 0.:dt:(tmax - dt)
        tspan = (0., tmax)

        # Frequency vector
        df = fs/bs
        fmax = useful_bandwidth
        freq = 0.:df:(fmax - df)
        freq_span = (0., fmax)

        return new(t, dt, tspan, freq, df, freq_span, fs, bs)
    end
end

"""
    tfestimate(input_signal::Vector{Float64}, output_signal::Vector{Float64}, fft_params::FFTParameters, window_input = hanning(fft_params.bs), window_output = window_input; overlap::Float64 = 0., type::Symbol = :h1)

Estimation of the one-sided transfer function between two signals

# Inputs
* `input_signal`: Input signal
* `output_signal`: Output signal
* `bs`: Size of a block
* `window_input`: Window function for the input signal
* `window_output`: Window function for the output signal
* `fs`: Sampling rate of the signal
* `overlap`: Overlap ratio between the segments
* `type`: Type of transfer function to estimate
    * :h1 (default)
    * :h2
    * :hv

# Outputs
* `H`: Transfer function
* `coh`: Coherence

# Available window functions

## From DSP.jl
* rect
* hanning
* hamming
* tukey
* cosine
* lanczos
* triang
* bartlett
* gaussian
* bartlett_hann
* blackman
* kaiser
* dpss

## From StructuralVibration.jl
* exponential
* force
* flattop
* nutall
* blackman_nutall
* parzen
* planck
"""
function tfestimation(input_signal::Vector{Float64}, output_signal::Vector{Float64}, bs, window_input = hanning(bs), window_output = window_input; fs = 1., overlap::Float64 = 0., type::Symbol = :h1)

    # FFT Parameters
    (; freq) = FFTParameters(fs, bs, pow2 = false)

    # Initialization
    if iseven(bs)
        N = bs ÷ 2 + 1
    else
        N = bs ÷ 2
    end
    Gxx = zeros(ComplexF64, N)
    Gyy = zeros(ComplexF64, N)
    Gxy = zeros(ComplexF64, N)
    H = undefs(ComplexF64, N)
    coh = undefs(Float64, N)

    # Signal segmentation
    input_segments, n_segments = signal_segmentation(input_signal, bs, window_input, overlap)
    output_segments = signal_segmentation(output_signal, bs, window_output, overlap)[1]

    # Auto and Cross-spectral densities estimation
    for (iseg, oseg) in zip(input_segments, output_segments)

        # FFT of the segments
        X = rfft(iseg)
        Y = rfft(oseg)

        # Spectral densities
        @. Gxx += abs2(X)
        @. Gyy += abs2(Y)
        @. Gxy += X*conj(Y)
    end

    # Averaging correction
    Gxx ./= n_segments
    Gyy ./= n_segments
    Gxy ./= n_segments

    # Correcting factors for energy due to windowing
    energy_correction_input = 1/rms(window_input)
    energy_correction_output = 1/rms(window_output)

    # Energy correction
    Gxx .*= energy_correction_input^2
    Gyy .*= energy_correction_output^2
    Gxy .*= energy_correction_input*energy_correction_output

    # Tranfer function estimation
    if type == :h1
        @. H = Gxy/Gxx
    elseif type == :h2
        @. H = Gyy/conj(Gxy)
    elseif type == :hv
        # Hv = sqrt(H1*H2)
        @. H = Gxy*sqrt(Gyy/Gxx)/abs(Gxy)
    end

    # Coherence calculation - coh = H1/H2
    @. coh = abs2(Gxy)/(Gxx*Gyy)

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = rfftfreq(bs, fs)

    useful_freqs = findall(freqs .<= freq[end])

    return H[useful_freqs], freqs[useful_freqs], coh[useful_freqs]
end

"""
    welch(input_signal::Vector{Float64}, fft_params::FFTParameters, window = hanning(fft_params.bs); overlap = 0.5, scale = :psd)

Estimation of the one-sided Power Spectral Density (PSD) of a signal using the Welch method

# Inputs
* `input_signal`: Input signal
* `bs`: Size of a block
* `window`: Window function
* `fs`: Sampling rate of the signal
* `overlap`: Overlap ratio between the segments
* `scaling`: Scale of the PSD - see https://community.sw.siemens.com/s/article/the-autopower-function-demystified for more information
    * :psd (default) - Autopower Power Spectral Density
    * :esd - Autopower Energy Spectral Density
    * :spectrum - Autopower spectrum
    * :linear - Autopower linear

# Output
* `pxx`: Autopower
* `freq`: Frequency vector

# Available window functions

## From DSP.jl
* rect
* hanning
* hamming
* tukey
* cosine
* lanczos
* triang
* bartlett
* gaussian
* bartlett_hann
* blackman
* kaiser
* dpss

## From StructuralVibration.jl
* exponential
* force
* flattop
* nutall
* blackman_nutall
* parzen
* planck
"""
function welch(input_signal, bs, window = hanning(bs); fs = 1, overlap = 0.5, scaling = :psd)
    # FFT Parameters
    (; tspan, freq) = FFTParameters(fs, bs, pow2 = false)

    # Signal segmentation
    signal, n_segments = signal_segmentation(input_signal, bs, window, overlap)

    # Initialize storage for periodograms of each segment
    if iseven(bs)
        N = bs ÷ 2 + 1
    else
        N = bs ÷ 2
    end

    Pxx = zeros(N)
    for segment in signal
        Pxx .+= abs2.(rfft(segment))
    end

    # Window energy - see https://community.sw.siemens.com/s/article/window-correction-factors
    scale = 1.
    if scaling == :psd || scaling == :esd
        # Explanation: By definition, the double-sided PSD is PSD = Pxx/Δf, where Pxx = abs2.(fft(x))/N^2, where N is the bs and Δf = fs/N. So, PSD = abs2.(fft(x))/(fs*N). Finally, this result must take into account the energy correction factor (ECF) to compensate for the windowing effect and preserve the signal energy. By definition, ECF = 1/rms(window) = sqrt(N/sum(abs2, window)). Therefore, PSD_corrected = PSD*ECF^2 = abs2.(fft(x))/(fs*sum(abs2, window)). Hence, the value of win_corr and scale...

        win_corr = sum(abs2, window)
        scale *= win_corr*fs

        if scaling == :esd
            # Explanation - ESD = T*PSD where T is the the acquisition block time
            scale /= tspan[2]
        end
    elseif scaling == :spectrum || scaling == :linear
        # Explanation: By definition, the double-sided autopower is Pxx = abs2.(fft(x))/N^2, where N is the bs. This result must take into account the amplitude correction factor (ECF) to compensate for the windowing effect and preserve the signal amplitude. By definition, ACF = 1/mean(window) = N/sum(window)). Therefore, Pxx_corrected = Pxx*ACF^2 = abs2.(fft(x))/sum(window)^2. Hence, the value of scale...
        scale = sum(window)^2
    end

    # Normalize the periodograms by the window energy and the sampling rate
    Pxx ./= (scale*n_segments)

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = rfftfreq(bs, fs)
    if iseven(bs)
        pxx = [Pxx[1]; 2Pxx[2:end-1]; Pxx[end]]
    else
        pxx = [Pxx[1]; 2Pxx[2:end]]
    end

    useful_freqs = findall(freqs .<= freq[end])

    if scaling == :linear
        @. pxx = sqrt(pxx)
    end

    return pxx[useful_freqs], freqs[useful_freqs]
end

"""
    spectrum(input_signal::Vector{Float64}, fft_params::FFTParameters, window = hanning(fft_params.bs); overlap = 0.5)

Estimation of the spectrum of a signal

# Inputs
* `input_signal`: Input signal
* `bs`: Size of a block
* `window`: Window function
* `fs`: Sampling rate of the signal
* `overlap`: Overlap ratio between the segments

# Output
* `y`: Spectrum
* `freq`: Frequency vector

# Available window functions

## From DSP.jl
* rect
* hanning
* hamming
* tukey
* cosine
* lanczos
* triang
* bartlett
* gaussian
* bartlett_hann
* blackman
* kaiser
* dpss

## From StructuralVibration.jl
* exponential
* force
* flattop
* nutall
* blackman_nutall
* parzen
* planck
"""
function spectrum(input_signal, bs, window = hanning(bs); fs = 1., overlap = 0.5)

    # FFT Parameters
    (; freq) = FFTParameters(fs, bs, pow2 = false)

    # Signal segmentation
    signal, n_segments = signal_segmentation(input_signal, bs, window, overlap)

    if iseven(bs)
        N = bs ÷ 2 + 1
    else
        N = bs ÷ 2
    end

    y = zeros(ComplexF64, N)
    for segment in signal
        y .+= rfft(segment)
    end

    # Averaging + Amplitude normalization + Amplitude correction factor (ACF) due to windowing - ACF = 1/mean(window)
    y ./= n_segments*bs*mean(window)

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = rfftfreq(bs, fs)
    if iseven(bs)
        y = [y[1]; 2y[2:end-1]; y[end]]
    else
        y = [y[1]; 2y[2:end]]
    end

    useful_freqs = findall(freqs .<= freq[end])

    return y[useful_freqs], freqs[useful_freqs]
end

"""
    signal_segmentation(signal, bs, window, overlap)

Segmentation of a signal into blocks

# Inputs
* `signal`: Signal to be segmented
* `bs`: Size of a block
* `window`: Window function
* `overlap`: Overlap ratio between the segments

# Output
* `segments`: Segmented signal
* `n_segments`: Number of segments
"""
function signal_segmentation(signal, bs, window, overlap)
    # Signal length
    N = length(signal)

    # bsength of the overlap segment
    noverlap = floor(Int, bs*overlap)
    step = bs - noverlap

    # Number of signal segments
    n_segments = step == 0 ? 1 : floor(Int, (N - bs)/step) + 1

    return [signal[(i - 1)*step .+ (1:bs)].*window for i in 1:n_segments], n_segments
end