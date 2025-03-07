#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
#| output: false
using StructuralVibration, ShareAdd
@usingany CairoMakie
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| echo: false
@doc SdofFreeTimeProblem
```
#
#
#
#
#
#
#
#
#
#| echo: false
@doc solve(prob::SdofFreeTimeProblem)
```
#
#
#
#
#
#| ouput: false

# Structural parameters
m = 1.0
f₀ = 1.

# Time vector
t = 0.:0.01:10.

# Initial conditions
u₀ = [0.1, -2]

# Undamped system
sdof_nd = Sdof(m, f₀, 0.)
prob_nd = SdofFreeTimeProblem(sdof_nd, u₀, t)
x_nd = solve(prob_nd).u

# Underdamped system
sdof_ud = Sdof(m, f₀, 0.1)
prob_ud = SdofFreeTimeProblem(sdof_ud, u₀, t)
x_ud = solve(prob_ud).u

# Critically damped system
sdof_cd = Sdof(m, f₀, 1.)
prob_cd = SdofFreeTimeProblem(sdof_cd, u₀, t)
x_cd = solve(prob_cd).u

# Overdamped system
sdof_od = Sdof(m, f₀, 2.)
prob_od = SdofFreeTimeProblem(sdof_od, u₀, t)
x_od = solve(prob_od).u;
#
#
#
#| echo: false

fig = Figure()
ax = Axis(fig[1, 1], ylabel="Displacement (m)", title = "Undamped system")
lines!(ax, t, x_nd)

ax1 = Axis(fig[1, 2], title = "Underdamped system")
lines!(ax1, t, x_ud)

ax2 = Axis(fig[2, 1], xlabel="Time (s)", ylabel="Displacement (m)", title = "Critically damped system")
lines!(ax2, t, x_cd)

ax3 = Axis(fig[2, 2], xlabel="Time (s)", title = "Overdamped system")
lines!(ax3, t, x_od)

display(fig);
#
#
#
