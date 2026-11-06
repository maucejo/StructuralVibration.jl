using StructuralVibration

## Detrending
# Signal 1
x1 = -0.5:0.01:0.5
y1 = @. sin(π*x1) + 0.25
y1_const = detrend(x1, y1, 0)
y1_lin = detrend(x1, y1, 1)

# Signal 2
x2 = 0:0.1:20
y2 = @. 3sin(x2) + x2
y2_const = detrend(x2, y2, 0)
y2_lin = detrend(x2, y2, 1)

## Gradient
# Signal
x = LinRange(0., 3π, 100)

# N = 30
# idx = round.(Int, (N - 1)rand(N) .+ 1) |> sort |> unique
# x = x[idx]

y = sin.(x)

# True gradient
dy = cos.(x)

# Estimated gradient
dy_approx = gradient(y, x)
## Decibel conversion
# Acceleration amplitude of 1 mm/s² (ref = 1 µm/s²) -> 60 dB
L_amp = to_dB(1e-3, 1e-6)

# Same quantity expressed as a power quantity (ref = (1 µm/s²)²) -> 60 dB
L_pow = to_dB((1e-3)^2, (1e-6)^2, type = :power)

# Complex FRF
L_frf = to_dB([1 + 1im, 10., 0.])

# Using the usual reference values
L_ref_amp = to_dB(1e-3, ref_dB(:acceleration))
L_ref_pow = to_dB((1e-3)^2, ref_dB(:acceleration, type = :power), type = :power)

## Amplitude of a signal
fs = 1000
t = (0:fs - 1)/fs
x = @. 2sin(2π*5t) + 0.5

x_peak = peak_value(x)  # 2.5
x_pp = peak2peak(x)     # 4.0
x_rms = rms_value(x)    # sqrt(2 + 0.25) = 1.5

# Multichannel signal
X = [x 2x]
X_peak = peak_value(X, dims = 1)
X_pp = peak2peak(X, dims = 1)
X_rms = rms_value(X, dims = 1)

## Phase unwrapping
freq = 0:0.5:100
τ = 0.02
H = @. exp(-2im*π*freq*τ)

ϕ_rad = unwrap_phase(H)              # -2πfτ
ϕ_deg = unwrap_phase(H, deg = true)  # -360fτ

# Multichannel signal (one channel per column)
Hm = [H H.^2]
ϕm = unwrap_phase(Hm)
