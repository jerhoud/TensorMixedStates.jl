# A GHZ state made by a circuit, a Hadamard gate and a ladder of CNOT gates, then each qubit
# depolarized with probability p. Its fidelity to the GHZ state and the statistics of the parity
# X⊗X⊗…⊗X, which a perfect GHZ state has equal to 1, are known in closed form, and checked
# here. The parity measured, keeping the result +1, post-selects a state closer to GHZ.

using TensorMixedStates, .Qubits, Random

n = 5
p = 0.1
sys = System(n, Qubit())

# the rightmost gate acts first: H on the first qubit, then CNOT(1, 2), CNOT(2, 3), ...
circuit = prod(controlled(X)(i, i + 1) for i in n-1:-1:1) * H(1)
noise = prod(depolarizing_gate(p)(i) for i in 1:n)
ρ = apply(noise, apply(circuit, State{Mixed}(sys, "Up")))
ghz = ghz_state(sys, "Up", "Dn")

# a qubit depolarized loses its coherence with the others: the fidelity keeps (1 - p)^n and what
# the incoherent mixtures still share with GHZ
exact_fidelity = (1 - p)^n + sum(binomial(n, k) * p^k * (1 - p)^(n - k) * 2.0^-(k + 1) for k in 1:n-1) +
                 p^n * 2.0^-n
parity = prod(X(i) for i in 1:n)
exact_plus = (1 + (1 - p)^n) / 2              # ⟨X⊗…⊗X⟩ = (1 - p)^n

println("fidelity to GHZ ", fidelity(ρ, ghz), "   exact ", exact_fidelity)
println("parity: ", probabilities(ρ, parity), "   exact +1 => ", exact_plus)
println("ten draws of Z(1): ", [ sample(ρ, Z(1); rng = Xoshiro(k)) for k in 1:10 ])
result, probability, kept = collapse(ρ, parity => 1)
println("post-selected on the parity +1, probability ", probability, ", fidelity ",
        fidelity(kept, ghz))
@assert abs(fidelity(ρ, ghz) - exact_fidelity) < 1e-10
@assert abs(last(last(probabilities(ρ, parity))) - exact_plus) < 1e-10
@assert fidelity(kept, ghz) > fidelity(ρ, ghz)
