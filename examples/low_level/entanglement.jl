# Entanglement, measured on states where it is known. A GHZ state has the entropy log 2 across
# every cut. A Bell pair depolarized with probability q is the Werner state
# (1 - q)|Φ⟩⟨Φ| + q I/4, whose reduced density matrix, logarithmic negativity log(2 - 3q/2) and
# Rényi-2 mutual information are known in closed form: on a mixed state the entropy of a part no
# longer measures its entanglement with the rest, which the negativity does.

using TensorMixedStates, .Qubits, .Circuits, LinearAlgebra

# a GHZ state of 6 qubits: the entanglement entropy across each of its cuts
ghz = ghz_state(System(6, Qubit()), "Up", "Dn")
entropies = [ first(entanglement_entropy(ghz, cut)) for cut in 1:5 ]   # the entropy, then its spectrum
println("GHZ, entropy across each cut ", entropies, ", log 2 = ", log(2))

# a Bell pair, then the depolarization of the two qubits together
q = 0.2
bell = apply(CX(1, 2) * H(1), State{Mixed}(System(2, Qubit()), "Up"))
werner = apply(depolarizing_gate(q, 2)(1, 2), bell)
values = measure(werner, [ReducedDensityMatrix([1, 2]), VonNeumannEntropy([1]),
                          LogNegativity([1], [2]), MutualInfoRenyi2([1])])
foreach(println, values)

Φ = [1, 0, 0, 1] / √2
exact_rdm = (1 - q) * Φ * Φ' + q * I(4) / 4
exact_purity = sum(abs2, exact_rdm)
@assert all(abs.(entropies .- log(2)) .< 1e-12)
@assert norm(last(values[1]) - exact_rdm) < 1e-12
@assert abs(last(values[2]) - log(2)) < 1e-12
@assert abs(last(values[3]) - log(2 - 3q / 2)) < 1e-12
@assert abs(last(values[4]) - (2log(2) + log(exact_purity))) < 1e-12
