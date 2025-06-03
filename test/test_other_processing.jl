using StructuralVibration
@usingany GLMakie

## Detrend
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

fig = Figure()
ax1 = Axis(fig[1, 1], ylabel = "signal", title = "sin(πt) + 0.25")
ax2 = Axis(fig[2, 1], xlabel = "x", ylabel = "signal", title = "3sin(x) + x")
hlines!(ax1, 0, color = :black, linewidth = 0.5)
vlines!(ax1, 0, color = :black, linewidth = 0.5)
lines!(ax1, x1, y1, label = "Original")
lines!(ax1, x1, y1_const, label = "Constant trend")
lines!(ax1, x1, y1_lin, label = "Linear trend")
xlims!(ax1, x1[1], x1[end])

hlines!(ax2, 0, color = :black, linewidth = 0.5)
lines!(ax2, x2, y2, label = "Original")
lines!(ax2, x2, y2_const, label = "Constant trend")
lines!(ax2, x2, y2_lin, label = "Linear trend")
xlims!(ax2, x2[1], x2[end])

Legend(fig[:, 2], ax1)

## Gradient
x = LinRange(0., 3π, 100)
y = sin.(x)

dy = cos.(x)
dy_approx = gradient(y, x)

fig1 = Figure()
ax1 = Axis(fig1[1, 1], xlabel = "x", ylabel = "signal")
lines!(ax1, x, dy, label = "cos(x)")
lines!(ax1, x, dy_approx, linestyle = :dash, label = "gradient")
xlims!(ax1, x[1], x[end])
axislegend(ax1, position = :rt)