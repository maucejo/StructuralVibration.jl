using Parameters, LinearAlgebra

includet("../src/models/bar_rod.jl")
includet("../src/models/beam.jl")
includet("../src/models/FEmodel.jl")

bar = Bar(1., 2.1e11, 7800., 0.3)
mesh = Mesh(bar, 0., 3, :CF)

K, M = assembly(bar, mesh)