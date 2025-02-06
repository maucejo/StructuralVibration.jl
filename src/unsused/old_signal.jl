"""
    FFTParameters

Structure to store the time and frequency parameters for the FFT analysis

# Constructor parameters
* `sample_rate`: Sample rate of the signal
* `block_size`: Block size of the signal

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
    sample_rate::Int
    block_size::Int

    function FFTParameters(sample_rate::Int, block_size::Int)
        # Sample rate
        sample_rate = nextpow(2, sample_rate)
        useful_bandwidth = sample_rate/2.56

        # Block size
        block_size = nextpow(2, block_size)

        # Time vector
        dt = 1/sample_rate
        tmax = block_size*dt
        t = 0.:dt:(tmax - dt)
        tspan = (0., tmax)

        # Frequency vector
        df = sample_rate/block_size
        fmax = useful_bandwidth
        freq = 0.:df:(fmax - df)
        freq_span = (0., fmax)

        return new(t, dt, tspan, freq, df, freq_span, sample_rate, block_size)
    end
end

"""
    tfestimate(input_signal::Vector{Float64}, output_signal::Vector{Float64}, fft_params::FFTParameters, window_input = hanning(fft_params.block_size), window_output = window_input; overlap_ratio::Float64 = 0., type_frf::Symbol = :h1)

Estimation of the transfer function between two signals

# Inputs
* `input_signal`: Input signal
* `output_signal`: Output signal
* `fft_params`: FFT parameters
* `window_input`: Window function for the input signal
* `window_output`: Window function for the output signal
* `overlap_ratio`: Overlap ratio between the segments
* `type_frf`: Type of transfer function to estimate
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
function tfestimation(input_signal::Vector{Float64}, output_signal::Vector{Float64}, fft_params::FFTParameters, window_input = hanning(fft_params.block_size), window_output = window_input; overlap_ratio::Float64 = 0., type_frf::Symbol = :h1)

    # FFT Parameters
    (; freq, block_size, sample_rate) = fft_params

    # Initialization
    Gxx = zeros(ComplexF64, block_size)
    Gyy = zeros(ComplexF64, block_size)
    Gxy = zeros(ComplexF64, block_size)
    H = undefs(ComplexF64, block_size)
    coh = undefs(Float64, block_size)

    # Signal segmentation
    input_segments, n_segments = signal_segmentation(input_signal, block_size, window_input, overlap_ratio)
    output_segments = signal_segmentation(output_signal, block_size, window_output, overlap_ratio)[1]

    # Auto and Cross-spectral densities estimation
    for (iseg, oseg) in zip(input_segments, output_segments)

        # FFT of the segments
        X = fft(iseg)
        Y = fft(oseg)

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
    if type_frf == :h1
        @. H = Gxy/Gxx
    elseif type_frf == :h2
        @. H = Gyy/conj(Gxy)
    elseif type_frf == :hv
        # Hv = sqrt(H1*H2)
        @. H = Gxy*sqrt(Gyy/Gxx)/abs(Gxy)
    end

    # Coherence calculation - coh = H1/H2
    @. coh = abs2(Gxy)/(Gxx*Gyy)

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = fftfreq(block_size, sample_rate)
    if iseven(block_size)
        idx_pos = 1:(block_size ÷ 2)
    else
        idx_pos = 1:(block_size ÷ 2 + 1)
    end

    useful_freqs = findall(freqs[idx_pos] .<= freq[end])

    return H[useful_freqs], freqs[useful_freqs], coh[useful_freqs]
end

"""
    welch(input_signal::Vector{Float64}, fft_params::FFTParameters, window = hanning(fft_params.block_size); overlap_ratio = 0.5, scale = :psd)

Estimation of the one-sided Power Spectral Density (PSD) of a signal using the Welch method

# Inputs
* `input_signal`: Input signal
* `fft_params`: FFT parameters
* `window`: Window function
* `overlap_ratio`: Overlap ratio between the segments
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
function welch(input_signal, fft_params::FFTParameters, window = hanning(fft_params.block_size); overlap_ratio = 0.5, scaling = :psd)
    # FFT Parameters
    (; tspan, freq, block_size, sample_rate) = fft_params

    # Signal segmentation
    signal, n_segments = signal_segmentation(input_signal, block_size, window, overlap_ratio)

    # Initialize storage for periodograms of each segment
    Pxx = zeros(block_size)
    for segment in signal
        Pxx .+= abs2.(fft(segment))
    end

    # Window energy - see https://community.sw.siemens.com/s/article/window-correction-factors
    scale = 1.
    if scaling == :psd || scaling == :esd
        # Explanation: By definition, the double-sided PSD is PSD = Pxx/Δf, where Pxx = abs2.(fft(x))/N^2, where N is the block_size and Δf = fs/N. So, PSD = abs2.(fft(x))/(fs*N). Finally, this result must take into account the energy correction factor (ECF) to compensate for the windowing effect and preserve the signal energy. By definition, ECF = 1/rms(window) = sqrt(N/sum(abs2, window)). Therefore, PSD_corrected = PSD*ECF^2 = abs2.(fft(x))/(fs*sum(abs2, window)). Hence, the value of win_corr and scale...

        win_corr = sum(abs2, window)
        scale *= win_corr*sample_rate

        if scaling == :esd
            # Explanation - ESD = T*PSD where T is the the acquisition block time
            scale /= tspan[2]
        end
    elseif scaling == :spectrum || scaling == :linear
        # Explanation: By definition, the double-sided autopower is Pxx = abs2.(fft(x))/N^2, where N is the block_size. This result must take into account the amplitude correction factor (ECF) to compensate for the windowing effect and preserve the signal amplitude. By definition, ACF = 1/mean(window) = N/sum(window)). Therefore, Pxx_corrected = Pxx*ACF^2 = abs2.(fft(x))/sum(window)^2. Hence, the value of scale...
        scale = sum(window)^2
    end

    # Normalize the periodograms by the window energy and the sampling rate
    Pxx ./= (scale*n_segments)

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = fftfreq(block_size, sample_rate)
    if iseven(block_size)
        idx_pos = 1:(block_size ÷ 2)
        pxx = [Pxx[1]; 2Pxx[2:(block_size ÷ 2)]; Pxx[block_size ÷ 2 + 1]]
    else
        idx_pos = 1:(block_size ÷ 2 + 1)
        pxx = [Pxx[1]; 2Pxx[2:((block_size + 1) ÷ 2)]]
    end

    useful_freqs = findall(freqs[idx_pos] .<= freq[end])

    if scaling == :linear
        @. pxx = sqrt(pxx)
    end

    return pxx[useful_freqs], freqs[useful_freqs]
end

"""
    spectrum(input_signal::Vector{Float64}, fft_params::FFTParameters, window = hanning(fft_params.block_size); overlap_ratio = 0.5)

Estimation of the spectrum of a signal

# Inputs
* `input_signal`: Input signal
* `fft_params`: FFT parameters
* `window`: Window function
* `overlap_ratio`: Overlap ratio between the segments

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
function spectrum(input_signal, fft_params::FFTParameters, window = hanning(fft_params.block_size); overlap_ratio = 0.5)

    # FFT Parameters
    (; freq, block_size, sample_rate) = fft_params

    # Signal segmentation
    signal, n_segments = signal_segmentation(input_signal, block_size, window, overlap_ratio)

    y = zeros(ComplexF64, block_size)
    for segment in signal
        y .+= fft(segment)
    end

    # Averaging + Amplitude normalization + Amplitude correction factor (ACF) due to windowing - ACF = 1/mean(window)
    y ./= n_segments*block_size*mean(window)

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = fftfreq(block_size, sample_rate)
    if iseven(block_size)
        idx_pos = 1:(block_size ÷ 2)
        y = [y[1]; 2y[2:(block_size ÷ 2)]; y[block_size ÷ 2 + 1]]
    else
        idx_pos = 1:(block_size ÷ 2 + 1)
        y = [y[1]; 2y[2:((block_size + 1) ÷ 2)]]
    end

    useful_freqs = findall(freqs[idx_pos] .<= freq[end])

    return y[useful_freqs], freqs[useful_freqs]
end

"""
    signal_segmentation(signal, block_size, window, overlap_ratio)

Segmentation of a signal into blocks

# Inputs
* `signal`: Signal to be segmented
* `block_size`: Size of the block
* `window`: Window function
* `overlap_ratio`: Overlap ratio between the segments

# Output
* `segments`: Segmented signal
* `n_segments`: Number of segments
"""
function signal_segmentation(signal, block_size, window, overlap_ratio)
    # Signal length
    N = length(signal)

    # block_sizeength of the overlap segment
    noverlap = floor(Int, block_size*overlap_ratio)
    step = block_size - noverlap

    # Number of signal segments
    n_segments = step == 0 ? 1 : floor(Int, (N - block_size)/step) + 1

    return [signal[(i - 1)*step .+ (1:block_size)].*window for i in 1:n_segments], n_segments
end