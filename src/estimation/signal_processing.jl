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
* flattop
* nutall
* blackman_nutall
* parzen
* planck
"""
function tfestimate(input_signal::Vector{Float64}, output_signal::Vector{Float64}, fft_params::FFTParameters, window_input = hanning(fft_params.block_size), window_output = window_input; overlap_ratio::Float64 = 0., type_frf::Symbol = :h1)

    # FFT Parameters
    (; freq, block_size, sample_rate) = fft_params

    # Signal length
    N = length(input_signal)

    # Length of the overlap segment
    noverlap = floor(Int, block_size*overlap_ratio)
    step = block_size - noverlap
    step == 0 ? step = 1 : nothing

    # Number of signal segments
    n_segments = floor(Int, (N - block_size)/step) + 1

    # Initialization
    Gxx = zeros(ComplexF64, block_size)
    Gyy = zeros(ComplexF64, block_size)
    Gxy = zeros(ComplexF64, block_size)
    H = undefs(ComplexF64, block_size)
    coh = undefs(Float64, block_size)

    # Correcting factors for energy due to windowing
    energy_correction_input = 1/rms(window_input)
    energy_correction_output = 1/rms(window_output)

    # Auto and Cross-spectral densities estimation
    for i in 1:n_segments
        # Indices defining the segment
        start_idx = (i-1)*step + 1
        end_idx = start_idx + block_size - 1

        # Segments extraction
        input_segment = input_signal[start_idx:end_idx].*window_input
        output_segment = output_signal[start_idx:end_idx].*window_output

        # FFT of the segments
        X = fft(input_segment)
        Y = fft(output_segment)

        # Spectral densities
        @. Gxx += abs2(X)
        @. Gyy += abs2(Y)
        @. Gxy += X*conj(Y)
    end

    # Averaging correction
    Gxx ./= n_segments
    Gyy ./= n_segments
    Gxy ./= n_segments

    # Energy correction
    Gxx .*= energy_correction_input
    Gyy .*= energy_correction_output
    Gxy .*= sqrt(energy_correction_input * energy_correction_output)

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

    # Extracting the positive frequencies
    freqs = fftfreq(block_size, sample_rate)
    if iseven(block_size)
        idx_pos = 1:(block_size ÷ 2)
    else
        idx_pos = 1:(block_size ÷ 2 + 1)
    end

    pos_useful_freqs = findall(freqs[idx_pos] .<= freq[end])

    return H[pos_useful_freqs], coh[pos_useful_freqs]
end