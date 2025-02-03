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
    freq
    df::Float64
    sample_rate::Int
    block_size::Int

    function FFTParameters(sample_rate::Int, block_size::Int)
        # Sample rate
        sample_rate = nextpow(2, sample_rate)
        useful_bandwidth = sample_rate/2.56

        # Block size
        block_size = nextpow(2, block_size)
        useful_spectral_lines = block_size/2.56

        # Time vector
        dt = 1/sample_rate
        tmax = block_size*dt
        t = 0:dt:(tmax - dt)

        # Frequency vector
        df = useful_bandwidth/useful_spectral_lines
        fmax = useful_bandwidth
        freq = 0:df:(fmax - df)

        return new(t, dt, freq, df, sample_rate, block_size)
    end
end

function tfestimate(input_signal::Vector{Float64}, output_signal::Vector{Float64}, fft_params::FFTParameters, window_input, window_output, overlap_ratio::Float64; type_frf::Symbol = :h1)
    # FFT Parameters
    (; freq, block_size, sample_rate) = fft_params

    # Signal length
    N = length(input_signal)

    # Segment length
    segment_length = block_size

    # Window generation
    win = window_input(segment_length)
    wout = window_output(segment_length)

    # Length of the overlap segment
    overlap = round(Int, overlap_ratio*segment_length)

    # Number of signal segments
    n_segments = floor(Int, (N - overlap)/(segment_length - overlap))

    # Initialization
    Gxx = zeros(ComplexF64, block_size)
    Gyy = zeros(ComplexF64, block_size)
    Gxy = zeros(ComplexF64, block_size)
    H = undefs(ComplexF64, block_size)
    coh = undefs(Float64, block_size)

    # Correcting factors for energy due to windowing
    energy_correction_input = 1/rms(win)
    energy_correction_output = 1/rms(wout)

    # Auto and Cross-spectral densities estimation
    for i in 1:n_segments
        # Indices defining the segment
        start_idx = (i-1)*(segment_length - overlap) + 1
        end_idx = start_idx + segment_length - 1

        # Segments extraction
        input_segment = input_signal[start_idx:end_idx].*win
        output_segment = output_signal[start_idx:end_idx].*wout

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
        @. H = Gxy*sqrt(Gxx/Gyy)/abs(Gxy)
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

    return H[pos_useful_freqs], coherence[pos_useful_freqs]
end