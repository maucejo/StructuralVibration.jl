# Export windows from DSP.jl
import DSP

# Rectangular window
"""
    rect

Rectangular window

Exported from DSP.jl
"""
rect = DSP.rect

# Hann window
"""
    hann

Hann(ing) window

Exported from DSP.jl
"""
hann = DSP.hanning

# Hamming window
"""
    hamming

Hamming window

Exported from DSP.jl
"""
hamming = DSP.hamming

# Tukey window
"""
    tukey

Tukey window

Exported from DSP.jl
"""
tukey = DSP.tukey

# Cosine window
"""
    cosine

Cosine window

Exported from DSP.jl
"""
cosine = DSP.cosine

# Lanczos window
"""
    lanczos

Lanczos window

Exported from DSP.jl
"""
lanczos = DSP.lanczos

# Triangle window
"""
    triang

Triangle window

Exported from DSP.jl
"""
triang = DSP.triang

# Bartlett window
"""
    bartlett

Bartlett window

Exported from DSP.jl
"""
bartlett = DSP.bartlett

# Gaussian window
"""
    gaussian

Gaussian window

Exported from DSP.jl
"""
gaussian = DSP.gaussian

# Bartlett-Hann window
"""
    bartlett_hann

Bartlett-Hann window

Exported from DSP.jl
"""
bartlett_hann = DSP.bartlett_hann

# Blackman window
"""
    blackman

Blackman window

Exported from DSP.jl
"""
blackman = DSP.blackman

# Kaiser window
"""
    kaiser

Kaiser window

Exported from DSP.jl
"""
kaiser = DSP.kaiser

# DPSS window
"""
    dpss

DPSS window

Exported from DSP.jl
"""


"""
    exponential(N, exponential_end = 0.01)

Create an exponential window

# Inputs
* `N`: Number of points
* `exponential_end`: End value of the exponential window

# Output
* `w`: Exponential window
"""
function exponential(N, exponential_end = 0.01)
    return exp.(log(exponential_end)*(0:N-1)/(N-1))
end

"""
    force(N, width = 0.1)

Create a force window

# Inputs
* `N`: Number of points
* `width`: Width of the force window

# Output
* `w`: Force window
"""
function force(N, width = 0.1)
    w = zeros(typeof(width), N)
    w[1:round(Int, width*N)] .= 1

    return w
end

"""
    flattop(N)

Create a flat top window

# Input
* `N`: Number of points

# Output
* `w`: Flattop window
"""
function flattop(N)
    n = 0:N-1
    a0 = 0.21557895
    a1 = 0.41663158
    a2 = 0.277263158
    a3 = 0.083578947
    a4 = 0.006947368

    return @. a0 - a1*cos(2π*n/(N-1)) + a2*cos(4π*n/(N-1)) .- a3*cos(6π*n/(N-1)) + a4*cos(8π*n/(N-1))
end

"""
    nutall(N)

Create a Nutall window

# Input
* `N`: Number of points

# Output
* `w`: Nutall window
"""
function nutall(N)
    n = 0:N-1
    a0 = 0.355768
    a1 = 0.487396
    a2 = 0.144232
    a3 = 0.012604

    return @. a0 - a1*cos(2π*n/(N-1)) + a2*cos(4π*n/(N-1)) .- a3*cos(6π*n/(N-1))
end

"""
    blackman_nutall(N)

Create a Blackman-Nutall window

# Input
* `N`: Number of points

# Output
* `w`: Blackman-Nutall window
"""
function blackman_nutall(N)
    n = 0:N-1
    a0 = 0.3635819
    a1 = 0.4891775
    a2 = 0.1365995
    a3 = 0.0106411

    return @. a0 - a1*cos(2π*n/(N-1)) + a2*cos(4π*n/(N-1)) - a3*cos(6π*n/(N-1))
end

"""
    parzen(N)

Create a Parzen window

# Input
* `N`: Number of points

# Output
* `w`: Parzen window
"""
function parzen(N)
    L = N + 1
    n = round.(Int, (0:N-1) .- N/2)

    w = similar(n, Float64)
    pos1 = @. 0 ≤ abs(n) ≤ L/4
    pos2 = @. L/4 < abs(n) ≤ L/2
    @. w[pos1] = 1 - 6*(2n[pos1]/L)^2*(1 - 2abs(n[pos1])/L)
    @. w[pos2] = 2*(1 - 2abs(n[pos2])/L)^3

    return w
end

"""
    planck(N, ϵ = 0.25)

Create a Planck-taper window

# Inputs
* `N`: Number of points
* `ϵ`: Parameter controlling the tapering

# Output
* `w`: Planck-taper window
"""
function planck(N, ϵ = 0.25)
    n = 0:(N-1)
    δ = ϵ*N
    w = similar(n, typeof(δ))

    pos1 = @. 1 ≤ n < δ
    pos2 = @. δ ≤ n < Int(N/2)
    pos3 = @. 0 ≤ n < Int(N/2)
    @. w[pos1] = 1/(1 + exp(δ/n[pos1] - δ/(δ - n[pos1])))
    @. w[pos2] = 1.
    @. w[N .- n[pos3]] = w[pos3]

    return w
end