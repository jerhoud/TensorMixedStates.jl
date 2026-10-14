# A domain wall, the left half of an XX chain up and the right half down, melting after a quench.
# The XX chain is free fermions in disguise: its magnetization is known exactly, from the
# exponential of the single particle hopping matrix. The evolution is computed by Tdvp and by
# ApproxW, and again with half the time step, which is how a time step is checked when nothing
# exact is known.

using TensorMixedStates, .Spins, LinearAlgebra

n = 20
t = 4.0
xx = sum(Sx(i) * Sx(i + 1) + Sy(i) * Sy(i + 1) for i in 1:n-1)
wall = State{Pure}(System(n, Spin(1/2)), [i ≤ n ÷ 2 ? "1/2" : "-1/2" for i in 1:n])
lim = Limits(maxdim = 64)

magnetization(st) = real(expect1(st, Sz))
tdvp = magnetization(evolve(-im * xx, t, wall; nsteps = 40, limits = lim))
tdvp_half = magnetization(evolve(-im * xx, t, wall; nsteps = 80, limits = lim))
approx = magnetization(evolve(-im * xx, t, wall; algo = ApproxW(order = 4), nsteps = 40,
                              limits = lim))

# exactly: an up spin is a fermion hopping by 1/2, ⟨Sz_i⟩ = Σ_j |U_ij|² - 1/2 over the j up at
# time 0, U = exp(-i h t)
h = SymTridiagonal(zeros(n), fill(0.5, n - 1))
u = exp(-im * t * Matrix(h))
exact = [ sum(abs2(u[i, j]) for j in 1:n÷2) - 0.5 for i in 1:n ]

println("magnetization at t = $t:  ", round.(tdvp; digits = 4))
println("largest error of Tdvp ", maximum(abs.(tdvp - exact)), ", with half the step ",
        maximum(abs.(tdvp_half - exact)), ", of ApproxW ", maximum(abs.(approx - exact)))
println("Tdvp against itself with half the step ", maximum(abs.(tdvp - tdvp_half)),
        ", against ApproxW ", maximum(abs.(tdvp - approx)))
@assert maximum(abs.(tdvp - exact)) < 1e-4
@assert maximum(abs.(approx - exact)) < 1e-4
