# Qubits driven by a field and jumping from "Up" to "Dn": their density matrix relaxes to a steady
# state, reached here twice, by a long evolution of the Lindblad equation and by steady_state,
# which searches the zero of the Lindbladian directly. The qubits being independent, the steady
# state is known: ⟨Z⟩ = -γ²/(γ² + 2Ω²).

using TensorMixedStates, .Qubits

n = 4
Ω = 0.8                    # the drive
γ = 0.5                    # the rate of the jumps from "Up" to "Dn" under Sm
lindbladian = sum(-im * (Ω / 2) * X(i) + Dissipator(sqrt(γ) * Sm)(i) for i in 1:n)
start = State{Mixed}(System(n, Qubit()), "Up")

relaxed = evolve(lindbladian, 60.0, start; nsteps = 120, limits = Limits(maxdim = 16))
residual, steady = steady_state(lindbladian, start; algo = Dmrg(nsweeps = 10),
                                limits = Limits(maxdim = 16))

exact = -γ^2 / (γ^2 + 2Ω^2)
println("⟨Z⟩ after a long evolution ", expect1(relaxed, Z))
println("⟨Z⟩ of the steady state    ", expect1(steady, Z), ", residual ", residual)
println("exact                      ", exact)
println("purity of the steady state ", purity(steady))
@assert maximum(abs.(expect1(relaxed, Z) .- exact)) < 1e-6
@assert maximum(abs.(expect1(steady, Z) .- exact)) < 1e-8
