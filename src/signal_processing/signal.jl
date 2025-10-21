"""
    FFTParameters

Structure to store the time and frequency parameters for the FFT analysis

**Constructor**
* `fs::Real`: Sampling rate of the signal
* `bs::Real`: Block size of the signal
* `pow2::Bool`: Flag for finding the next power of 2

**Fields**
* `t::AbstractRange`: Time vector
* `dt::Float64`: Time step
* `tspan::Tuple{Float64, Float64}`: Time span
* `freq::AbstractRange`: Frequency vector
* `df::Float64`: Frequency resolution
* `freq_span::Tuple{Float64, Float64}`: Frequency span
* `fs::Int`: Sampling rate
* `bs::Int`: Block size
"""
@show_data struct FFTParameters{Tt <: AbstractRange, Tf <: AbstractRange, T <: Real}
    t::Tt
    dt::T
    tspan::Tuple{T, T}
    freq::Tf
    df::T
    freq_span::Tuple{T, T}
    fs::Int
    bs::Int

    function FFTParameters(fs::Int, bs::Int; pow2 = false)
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
        tspan = (0., tmax - dt)

        # Frequency vector
        df = fs/bs
        fmax = useful_bandwidth
        freq = 0.:df:(fmax - df)
        freq_span = (0., fmax - df)

        return new{typeof(t), typeof(freq), typeof(dt)}(t, dt, tspan, freq, df, freq_span, fs, bs)
    end
end

"""
    tfestimate(input_signal::AbstractMatrix, output_signal::AbstractMatrix, bs::Int, window_input = hanning, window_output = window_input; fs::Int = 1, overlap = 0., type = :h1)

    tfestimate(input_signal::AbstractVector, output_signal::AbstractVector, bs::Int, window_input = hanning, window_output = window_input; fs::Int = 1, overlap = 0., type = :h1)

Estimation of the one-sided transfer function between two signals

**Inputs**
* `input_signal::AbstractVector` or `AbstractMatrix`: Input signal
* `output_signal::AbstractVector` or `AbstractMatrix`: Output signal
* `bs::Int` Block size
* `window_input`: Window function for the input signal (default: hanning)
* `window_output`: Window function for the output signal (default: same as window_input)
* `fs::Int`: Sampling rate (default: 1)
* `overlap::Real`: Overlap ratio between the segments (default: 0.)
* `nfft::Int`: Number of FFT points (default: bs)
* `type::Symbol`: Type of transfer function to estimate
    * `:h1` (default)
    * `:h2`
    * `:h3` - h3 = (h1 + h2)/2
    * `:hv` - hv = sqrt(h1*h2)

**Outputs**
* `H`: Transfer function
* `freq`: Frequency range
* `coh`: Coherence
"""
function tfestimate(input_signal::AbstractMatrix, output_signal::AbstractMatrix, bs::Int, window_input = hanning, window_output = window_input; fs::Int = 1, overlap = 0., nfft::Int = bs, type = :h1)

    ni, nti = size(input_signal)
    no, nto = size(output_signal)

    if nti != nto
        throw(DimensionMismatch("Input and output signals must have the same number of time samples"))
    end

    # Estimate the number of fft points
    (; freq) = FFTParameters(fs, bs)
    m = nfft > bs ? nfft : bs
    freqs = rfftfreq(m, fs)
    useful_freqs = findall(freqs .<= freq[end])
    nf = length(useful_freqs)

    H = similar(complex.([1.], [1.]), no, ni, nf)
    coh = similar(real(H))
    for (i, xi) in enumerate(eachrow(input_signal))
        for (j, yj) in enumerate(eachrow(output_signal))
            H[j, i, :], _, coh[j, i, :] = tfestimate(xi, yj, bs, window_input, window_output; fs = fs, overlap = overlap, nfft = nfft, type = type)
        end
    end

    return H, freqs[useful_freqs], coh
end

function tfestimate(input_signal::AbstractVector, output_signal::AbstractVector, bs::Int, window_input = hanning, window_output = window_input; fs::Int = 1, overlap = 0., nfft::Int = bs, type = :h1)

    # FFT Parameters
    (; freq) = FFTParameters(fs, bs)

    # Check if the overlap ratio is between 0 and 1
    0. ≤ overlap ≤ 1. ? nothing : throw(DomainError("overlap ratio must be between 0 and 1"))

    m = nfft > bs ? nfft : bs

    # Check if the number of fft points is even or odd
    n = iseven(m) ? m ÷ 2 + 1 : m ÷ 2

    Gxx = zeros(Complex{eltype(freq)}, n)
    Gyy = zeros(Complex{eltype(freq)}, n)
    Gxy = zeros(Complex{eltype(freq)}, n)
    X = similar(Gxx)
    Y = similar(Gxx)
    H = similar(Gxx)
    coh = similar(real(Gxx))

    # Signal segmentation
    input_segments, n_segments = signal_segmentation(input_signal, bs, window_input(bs), overlap)
    output_segments = signal_segmentation(output_signal, bs, window_output(bs), overlap)[1]

    # Auto and Cross-spectral densities estimation
    for (iseg, oseg) in zip(input_segments, output_segments)

        # FFT of the segments
        if nfft > bs
            X .= rfft(zpad(iseg, nfft))
            Y .= rfft(zpad(oseg, nfft))
        else
            X .= rfft(iseg)
            Y .= rfft(oseg)
        end

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
    energy_correction_input = 1/rms(window_input(bs))
    energy_correction_output = 1/rms(window_output(bs))

    # Energy correction
    Gxx .*= energy_correction_input^2
    Gyy .*= energy_correction_output^2
    Gxy .*= energy_correction_input*energy_correction_output

    # Tranfer function estimation
    if type == :h1
        @. H = Gxy/Gxx
    elseif type == :h2
        @. H = Gyy/conj(Gxy)
    elseif type == :h3
        # H3 = (H1 + H2)/2 - Arithmetic mean
        @. H = (abs2(Gxy) + Gyy*Gxx)/2(Gxx*conj(Gxy))
    elseif type == :hv
        # Hv = sqrt(H1*H2) - Geometric mean
        @. H = Gxy*sqrt(Gyy/Gxx)/abs(Gxy)
    end

    # Coherence calculation - coh = H1/H2
    @. coh = abs2(Gxy)/(Gxx*Gyy)

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = rfftfreq(m, fs)

    useful_freqs = findall(freqs .<= freq[end])

    return H[useful_freqs], freqs[useful_freqs], coh[useful_freqs]
end

"""
    welch(input_signal, bs, window = hanning; fs = 1,
          overlap = 0.5, nfft = bs, scaling = :psd)

Estimation of one-sided Autopower functions of a signal using the Welch method

**Inputs**
* `input_signal::AbstractVector`: Input signal
* `bs::Int`: Block size
* `window`: Window function (default: hanning)
* `fs::Int`: Sampling rate (default: 1)
* `overlap`: Overlap ratio between the segments (default: 0.5)
* `nfft::Int`: Number of FFT points (default: bs)
* `scaling`: Scale of the PSD - see https://community.sw.siemens.com/s/article/the-autopower-function-demystified for more information
    * `:psd` (default) - Power Spectral Density
    * `:esd` - Autopower Energy Spectral Density
    * `:spectrum` - Autopower spectrum
    * `:linear` - Autopower linear

**Outputs**
* `pxx`: Autopower
* `freq`: Frequency range

**Notes**
- The `welch`` function is already implemented in `DSP.jl` under the name `welch_pgram`. The function `welch` is implemented here for pedagogical purposes.
- The `welch` function is equivalent to calling `csd(x, x, ...)`.
"""
function welch(input_signal::AbstractVector, bs::Int, window = hanning; fs::Int = 1, overlap = 0.5, nfft::Int = bs, scaling = :psd)

    # FFT Parameters
    (; tspan, freq) = FFTParameters(fs, bs, pow2 = false)

    # Check if the overlap ratio is between 0 and 1
    0. ≤ overlap ≤ 1. ? nothing : throw(DomainError("overlap ratio must be between 0 and 1"))

    m = nfft > bs ? nfft : bs

    # Check if the number of fft points is even or odd
    n = iseven(m) ? m ÷ 2 + 1 : m ÷ 2

    # Signal segmentation
    signal, n_segments = signal_segmentation(input_signal, bs, window(bs), overlap)

    Pxx = zeros(eltype(input_signal), n)
    for segment in signal
        if nfft > bs
            Pxx .+= abs2.(rfft(zpad(segment, nfft)))
        else
            Pxx .+= abs2.(rfft(segment))
        end
    end

    # Window energy - see https://community.sw.siemens.com/s/article/window-correction-factors
    scale = 1.
    if scaling == :psd || scaling == :esd
        # Explanation: By definition, the double-sided PSD is PSD = Pxx/Δf, where Pxx = abs2.(fft(x))/n^2, where n is the bs and Δf = fs/n. So, PSD = abs2.(fft(x))/(fs*n). Finally, this result must take into account the energy correction factor (ECF) to compensate for the windowing effect and preserve the signal energy. By definition, ECF = 1/rms(window) = sqrt(n/sum(abs2, window)). Therefore, PSD_corrected = PSD*ECF^2 = abs2.(fft(x))/(fs*sum(abs2, window)). Hence, the value of win_corr and scale...

        win_corr = sum(abs2, window(bs))
        scale *= win_corr*fs

        # Explanation - ESD = T*PSD where T is the the acquisition block time
        scaling == :esd ? scale /= tspan[2] : nothing

    elseif scaling == :spectrum || scaling == :linear
        # Explanation: By definition, the double-sided autopower is Pxx = abs2.(fft(x))/n^2, where n is the bs. This result must take into account the amplitude correction factor (ACF) to compensate for the windowing effect and preserve the signal amplitude. By definition, ACF = 1/mean(window) = n/sum(window)). Therefore, Pxx_corrected = Pxx*ACF^2 = abs2.(fft(x))/sum(window)^2. Hence, the value of scale...
        scale = sum(window(bs))^2
    end

    # normalize the periodograms by the window energy and the sampling rate
    Pxx ./= (scale*n_segments)

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = rfftfreq(m, fs)

    # Correcting the amplitude of the spectrum
    pxx = iseven(m) ? [Pxx[1]; 2Pxx[2:end-1]; Pxx[end]] : [Pxx[1]; 2Pxx[2:end]]

    useful_freqs = findall(freqs .<= freq[end])

    # Linear scaling
    if scaling == :linear
        @. pxx = sqrt(pxx)
    end

    return pxx[useful_freqs], freqs[useful_freqs]
end

"""
    csd(x, y, bs, window_x = hanning, window_y = window_x; fs = 1,
        overlap = 0.5, nfft = bs, scaling = :psd)

Estimation of the one-sided Cross-spectral density between two signals

**Inputs**
* `x::AbstractVector`: First signal
* `y::AbstractVector`: Second signal
* `bs::Int`: Block size
* `window_x`: Window function for the first signal (default: hanning)
* `window_y`: Window function for the second signal (default: same as window_x)
* `fs::Int`: Sampling rate (default: 1)
* `overlap`: Overlap ratio between the segments
* `nfft::Int`: Number of FFT points (default: bs)
* `scaling`: Scale of the CSD
    * `:psd` (default) - Power Spectral Density
    * `:esd` - Energy Spectral Density
    * `:spectrum` - Spectrum
    * `:linear` - Linear

**Outputs**
* `pxy`: Cross-spectral density
* `freq`: Frequency range

**Notes**
- The `csd(x, x, ...) = welch(x, ...)`.
"""
function csd(x::AbstractMatrix, y::AbstractMatrix, bs::Int, window_x = hanning, window_y = window_x; fs::Int = 1, overlap = 0.5, nfft::Int = bs, scaling = :psd)

    ni, nti = size(x)
    no, nto = size(y)

    if nti != nto
        throw(DimensionMismatch("Input and output signals must have the same number of time samples"))
    end

    # Estimate the number of fft points
    (; freq) = FFTParameters(fs, bs)
    m = nfft > bs ? nfft : bs
    freqs = rfftfreq(m, fs)
    useful_freqs = findall(freqs .<= freq[end])
    nf = length(useful_freqs)

    Pxy = similar(complex.([1.], [1.]), no, ni, nf)
    for (i, xi) in enumerate(eachrow(x))
        for (j, yj) in enumerate(eachrow(y))
            Pxy[j, i, :] = csd(xi, yj, bs, window_x, window_y; fs = fs, overlap = overlap, nfft = nfft, scaling = scaling)[1]
        end
    end

    return Pxy, freqs[useful_freqs]
end

function csd(x::AbstractVector, y::AbstractVector, bs::Int, window_x = hanning, window_y = window_x; fs::Int = 1, overlap = 0.5, nfft::Int = bs, scaling = :psd)

    # FFT Parameters
    (; tspan, freq) = FFTParameters(fs, bs)

    # Check if the overlap ratio is between 0 and 1
    0. ≤ overlap ≤ 1. ? nothing : throw(DomainError("overlap ratio must be between 0 and 1"))

    m = nfft > bs ? nfft : bs

    # Check if the block size is even or odd
    n = iseven(m) ? m ÷ 2 + 1 : m ÷ 2

    # Signal segmentation
    xsig, n_segments = signal_segmentation(x, bs, window_x(bs), overlap)
    ysig = signal_segmentation(y, bs, window_y(bs), overlap)[1]

    Pxy = zeros(Complex{eltype(freq)}, n)
    X = similar(Pxy)
    Y = similar(Pxy)
    # Auto and Cross-spectral densities estimation
    for (xs, ys) in zip(xsig, ysig)

        # FFT of the segments
        if nfft > bs
            X .= rfft(zpad(xs, nfft))
            Y .= rfft(zpad(ys, nfft))
        else
            X .= rfft(xs)
            Y .= rfft(ys)
        end

        # Spectral densities
        @. Pxy += X*conj(Y)
    end

    # Window energy - see https://community.sw.siemens.com/s/article/window-correction-factors
    scale = 1.
    if scaling == :psd || scaling == :esd
        # Explanation: By definition, the double-sided PSD is PSD = Pxx/Δf, where Pxx = abs2.(fft(x))/n^2, where n is the bs and Δf = fs/n. So, PSD = abs2.(fft(x))/(fs*n). Finally, this result must take into account the energy correction factor (ECF) to compensate for the windowing effect and preserve the signal energy. By definition, ECF = 1/rms(window) = sqrt(n/sum(abs2, window)). Therefore, PSD_corrected = PSD*ECF^2 = abs2.(fft(x))/(fs*sum(abs2, window)). Hence, the value of win_corr and scale...

        win_corr = bs*rms(window_x(bs))*rms(window_y(bs))
        scale *= win_corr*fs

        # Explanation - ESD = T*PSD where T is the the acquisition block time
        scaling == :esd ? scale /= tspan[2] : nothing

    elseif scaling == :spectrum || scaling == :linear
        # Explanation: By definition, the double-sided autopower is Pxx = abs2.(fft(x))/n^2, where n is the bs. This result must take into account the amplitude correction factor (ACF) to compensate for the windowing effect and preserve the signal amplitude. By definition, ACF = 1/mean(window) = n/sum(window)). Therefore, Pxx_corrected = Pxx*ACF^2 = abs2.(fft(x))/sum(window)^2. Hence, the value of scale...
        scale = sum(window_x(bs))*sum(window_y(bs))
    end

    # normalize the periodograms by the window energy and the sampling rate
    Pxy ./= (scale*n_segments)

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = rfftfreq(m, fs)

    # Correcting the amplitude of the spectrum
    pxy = iseven(m) ? [Pxy[1]; 2Pxy[2:end-1]; Pxy[end]] : [Pxy[1]; 2Pxy[2:end]]

    useful_freqs = findall(freqs .<= freq[end])

    # Linear scaling
    if scaling == :linear
        @. pxy = sqrt(pxy)
    end

    return pxy[useful_freqs], freqs[useful_freqs]
end

"""
    spectrum(input_signal, bs, window = hanning; fs = 1, overlap = 0.5,
             nfft = bs)

Estimation of the spectrum of a signal

**Inputs**
* `input_signal::Vector{Real}`: Input signal
* `bs::Int`: Block size
* `window`: Window function (default: hanning)
* `fs::Int`: Sampling rate
* `overlap`: Overlap ratio between the segments
* `nfft::Int`: Number of FFT points (default: bs)

**Outputs**
* `y`: Signal spectrum
* `freq`: Frequency range
"""
function spectrum(input_signal::Vector{T}, bs::Int, window = hanning; fs::Int = 1, overlap = 0.5, nfft::Int = bs) where {T <: Real}

    # FFT Parameters
    (; freq) = FFTParameters(fs, bs, pow2 = false)

    # Check if the overlap ratio is between 0 and 1
    0. ≤ overlap ≤ 1. ? nothing : throw(DomainError("overlap ratio must be between 0 and 1"))

    # Signal segmentation
    signal, n_segments = signal_segmentation(input_signal, bs, window(bs), overlap)

    m = nfft > bs ? nfft : bs

    # Check if the block size is even or odd
    n = iseven(m) ? m ÷ 2 + 1 : m ÷ 2

    y = zeros(Complex{eltype(input_signal)}, n)
    for segment in signal
        if nfft > bs
            y .+= rfft(zpad(segment, nfft))
        else
            y .+= rfft(segment)
        end
    end

    # Averaging + Amplitude normalization + Amplitude correction factor (ACF) due to windowing - ACF = 1/mean(window)
    y ./= n_segments*bs*mean(window(bs))

    # Convert full-spectrum to one-sided spectrum
    # Extracting the positive frequencies
    freqs = rfftfreq(m, fs)

    # Correcting the amplitude of the spectrum
    y .= iseven(m) ? [y[1]; 2y[2:end-1]; y[end]] : [y[1]; 2y[2:end]]

    useful_freqs = findall(freqs .<= freq[end])

    return y[useful_freqs], freqs[useful_freqs]
end

"""
    signal_segmentation(signal, bs, window, overlap)

Segmentation of a signal into blocks

**Inputs**
* `signal`: Signal to be segmented
* `bs::Int`: Size of a block
* `window`: Window function
* `overlap`: Overlap ratio between the segments

**Outputs**
* `segments`: Segmented signal
* `n_segments`: number of segments
"""
function signal_segmentation(signal, bs::Int, window, overlap)
    # Signal length
    n = length(signal)

    # bsength of the overlap segment
    noverlap = floor(Int, bs*overlap)
    step = bs - noverlap

    # number of signal segments
    n_segments = step == 0 ? 1 : floor(Int, (n - bs)/step) + 1

    return [signal[(i - 1)*step .+ (1:bs)].*window for i in 1:n_segments], n_segments
end

"""
    anti_alias(signal, fc; fs = 2fc)

Apply an anti-aliasing filter to a signal

**Inputs**
* `signal`: Signal to be filtered
* `fc`: Cut-off frequency
* `fs`: Sampling rate

**Output**
* `signal`: Filtered signal
"""
function anti_alias(signal, fc; fs = 2fc)
    # nyquist frequency
    fn = fs/2

    # Filter design
    order = 200        # Filter order
    fb = 0.05*fc       # Filter bandwith = 2*fb
    freq_filt = [(0., fc - fb) => 1, (fc, fn) => 0]
    filt_coeff = remez(order, freq_filt, Hz = fs, maxiter = 50)

    return filtfilt(filt_coeff, signal)
end

"""
    zpad(x::AbstractVector, N)
    zpad(x::AbstractMatrix, N)

Zero-pad a vector to a specified length N for usage in RFFT/IRFFT

**Inputs**
- `x`: Input vector or matrix
- `N`: Desired length after zero-padding
- `type`: Padding type
    * `:dir` (default): Zero-padding for direct RFFT
    * `:inv`: Zero-padding for inverse RFFT

**Output**
- `xpad`: Zero-padded vector or matrix

**Notes**
- If the length of `x` is already greater than or equal to `N`, the original `x` is returned.
- The function supports both vectors and matrices, padding along the last dimension.
"""
function zpad(x::AbstractVector, N, type = :dir)
    nx = length(x)

    if N > nx
        if type == :inv
            if N%2 == 0
                M = Int(round(N/2 + 1))
            else
                M = Int(round((N + 1)/2))
            end

            xpad = [x; zeros(eltype(x), M - nx)]
        else
            xpad = [x; zeros(eltype(x), N - nx)]
        end
    else
        return x
    end

    return xpad
end

function zpad(x::AbstractMatrix, N, type = :dir)
    nx, ny = size(x)

    if N > ny
        if type == :inv
            if N%2 == 0
                M = Int(round(N/2 + 1))
            else
                M = Int(round((N+1)/2))
            end

            xpad = [x zeros(eltype(x), nx, M - ny)]
        else
            xpad = [x zeros(eltype(x), nx, N - ny)]
        end
    else
        return x
    end

    return xpad
end