using StructuralVibration, LinearAlgebra, Test

## Second-order convergence of the Newmark family
# System parameters
M = [2. 0.; 0. 1.]
K = [3000. -1000.; -1000. 1000.]
C = 0.01M + 1e-4K

# Initial conditions
x0 = [1e-2, -2e-2]
v0 = [0.1, 0.3]
u0 = (x0, v0)

# External forces - F(t) = f0*sin(Ω*t)
f0 = [10., 0.]
Ω = 15.

# Final time
T = 1.

# Exact solution - The harmonic force is added to the state vector (x, v, s, c)
A = zeros(6, 6)
A[1:2, 3:4] .= I(2)
A[3:4, 1:2] .= -M\K
A[3:4, 3:4] .= -M\C
A[3:4, 5] .= M\f0
A[5, 6] = Ω
A[6, 5] = -Ω
x_exact = (exp(A*T)*[x0; v0; 0.; 1.])[1:2]

# Relative error at the final time
function relative_error(alg, h)
    t = 0.:h:T
    F = f0*sin.(Ω*t)'
    x = solve(DirectTimeProblem(K, M, C, F, u0, t), alg, progress = false).u

    return norm(x[:, end] - x_exact)/norm(x_exact)
end

algs = [
    "Newmark" => Newmark(),
    "FoxGoodwin" => FoxGoodwin(),
    "LinearAcceleration" => LinearAcceleration(),
    "MidPoint" => MidPoint(),
    "HHT(ρ = 1)" => HHT(),
    "HHT(ρ = 0.8)" => HHT(ρ = 0.8),
    "HHT(ρ = 0.5)" => HHT(ρ = 0.5),
    "HHT(αf = 0.1)" => HHT(αf = 0.1),
    "WBZ(ρ = 1)" => WBZ(),
    "WBZ(ρ = 0.8)" => WBZ(ρ = 0.8),
    "WBZ(ρ = 0.5)" => WBZ(ρ = 0.5),
    "WBZ(αm = 0.2)" => WBZ(αm = 0.2),
    "GeneralizedAlpha(ρ = 1)" => GeneralizedAlpha(),
    "GeneralizedAlpha(ρ = 0.8)" => GeneralizedAlpha(ρ = 0.8),
    "GeneralizedAlpha(ρ = 0.5)" => GeneralizedAlpha(ρ = 0.5),
    "GeneralizedAlpha(αf = 0.4, αm = 0.3)" => GeneralizedAlpha(αf = 0.4, αm = 0.3)
]

@testset "Second-order convergence" begin
    @testset "$name" for (name, alg) in algs
        order = log2(relative_error(alg, 5e-4)/relative_error(alg, 2.5e-4))

        @test order ≈ 2 atol = 0.05
    end
end
