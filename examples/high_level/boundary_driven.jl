# A chain of free fermions between two baths, the left one filling site 1 at the rate b1 and
# emptying it at the rate a1, the right one site L at the rates bL and aL: its steady state,
# computed exactly by GaussianState at L = 100, against its closed form, a density flat inside the
# chain and a current constant along it. The structure is that of Karevski and Platini (2009),
# derived again for these rates: with the fillings ν = b / (a + b) the baths impose and their
# strengths g = (a + b) / 2t, the correlation ⟨c_j† c_{j+1}⟩ = i y with
# y = (νL - ν1) / (g1 + 1/g1 + gL + 1/gL), and the current is -2t y.

using TensorMixedStates, .Fermions

L, t = 100, 1.0
a1, b1, aL, bL = 0.3, 0.7, 0.9, 0.1
hop = t * sum(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j) for j in 1:L-1)
baths = a1 * Dissipator(C)(1) + b1 * Dissipator(dag(C))(1) +
        aL * Dissipator(C)(L) + bL * Dissipator(dag(C))(L)
current(j) = im * t * (dag(C)(j) * C(j + 1) - dag(C)(j + 1) * C(j))

sim = runTMS(SimData(
    name = "boundary_driven",
    description = "free fermions between two baths, L = $L, steady state",
    backend = GaussianState,
    phases = [
        CreateState{Mixed}(System(L, Fermion()), "FullyMixed"),
        # algo, which a search on MPS needs, has no use on a Gaussian state, as the log notes
        SteadyState(lindbladian = -im * hop + baths, algo = Dmrg(nsweeps = 10),
                    final_measurements = "profile.dat" => [N, "current" => current(L ÷ 2)]),
    ],
))

ν1, νL = b1 / (a1 + b1), bL / (aL + bL)
g1, gL = (a1 + b1) / 2t, (aL + bL) / 2t
y = (νL - ν1) / (g1 + 1 / g1 + gL + 1 / gL)
density = [ν1 + y / g1; fill(ν1 + y * (g1 + 1 / g1), L - 2); νL - y / gL]
println("density inside $(density[2]), current $(-2t * y)")
@assert maximum(abs(expect(sim.state, N(j)) - density[j]) for j in 1:L) < 1e-12
@assert abs(expect(sim.state, current(L ÷ 2)) + 2t * y) < 1e-12
