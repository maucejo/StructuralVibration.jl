using Parameters, DSP, LinearAlgebra, Interpolations
@usingany GLMakie

includet("../../../src/models/sdof.jl")

## SDOF system
m = 1.
ω₀ = 2π*10.
ξ = 0.
sdof = SDOF(m, ω₀, ξ)

#Time vector
Δt = 1e-3
t = 0.:Δt:1.

## Harmonic response
freq_exc = 9.
x = forced_response(sdof, 1., 2π*freq_exc, t)[1]

## Plot
pos_begin = 0. .≤ t .≤ 0.342
pos_bm = 0.342 .≤ t .≤ 0.395
pos_middle = 0.395 .≤ t .≤ 0.604
pos_me = 0.604 .≤ t .≤ 0.659
pos_end = 0.659 .≤ t .≤ 1.

tbegin = t[pos_begin]
xbegin = x[pos_begin]

# For color gradient
tbm = t[pos_bm]
xbm = x[pos_bm]
nbm = length(tbm)
drbm = (0.796 - 0.22)/(nbm - 1)
dgbm = (0.596 - 0.235)/(nbm - 1)
dbbm = (0.2 - 0.149)/(nbm - 1)
colbm = RGBf.(0.796:-drbm:0.22, 0.235:dgbm:0.596, 0.2:-dbbm:0.149)

tmiddle = t[pos_middle]
xmiddle = x[pos_middle]

# For color gradient
tme = t[pos_me]
xme = x[pos_me]
nme = length(tme)
drme = (0.584 - 0.22)/(nme - 1)
dgme = (0.596 - 0.345)/(nme - 1)
dbme = (0.698 - 0.149)/(nme - 1)
colme = RGBf.(0.22:drme:0.584, 0.596:-dgme:0.345, 0.149:dbme:0.698)

tend = t[pos_end]
xend = x[pos_end]

# Plot
lw = 7.
fig = Figure()
ax = Axis(fig[1, 1])
hidespines!(ax)
hidedecorations!(ax)
lines!(ax, tbegin, xbegin, color = RGBf(0.796, 0.235, 0.2), linewidth = lw)
lines!(ax, tbm, xbm, color = colbm, linewidth = lw, linecap = :round)
lines!(ax, tmiddle, xmiddle, color = RGBf(0.22, 0.596, 0.149), linewidth = lw)
lines!(ax, tme, xme, color = colme, linewidth = lw, linecap = :round)
lines!(ax, tend, xend, color = RGBf(0.584, 0.345, 0.698), linewidth = lw)
fig