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
y = sin.(x)

# True gradient
dy = cos.(x)

# Estimated gradient
dy_approx = gradient(y, x)