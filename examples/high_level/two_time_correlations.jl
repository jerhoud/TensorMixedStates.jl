# The correlations at two times ⟨X_4(0) X_j(t)⟩ in the ground state of a transverse field Ising
# chain, for every site j: the ground state |ψ⟩ is kept under a name, the operator X_4 applied
# to a copy |φ⟩ = X_4|ψ⟩, and the two states evolved together by one Evolve acting on both,
# Overlap(:φ, X) measuring ⟨φ(t)|X_j|ψ(t)⟩ = ⟨ψ|X_4 e^{iHt} X_j e^{-iHt}|ψ⟩ on every site j.
# The last one is checked against the exact evolution of the chain, small enough here to be
# written as a dense matrix.

using TensorMixedStates, .Qubits, LinearAlgebra

n = 8
ising = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - 1.2 * sum(X(i) for i in 1:n)
limits = Limits(maxdim = 64)
t = 2.0

sim = runTMS(SimData(
    name = "correlations_$n",
    description = "the correlations <X_4(0) X_j(t)> in the ground state of an Ising chain",
    seed = 1,
    phases = [
        CreateState(type = Pure, system = System(n, Qubit()), randomize = 8),
        GroundState(hamiltonian = ising, algo = Dmrg(nsweeps = 10, tol = 1e-12), limits = limits),
        CopyState(:state => :φ),
        Gates(gates = X(4), on = :φ),
        Evolve(evolver = -im * ising, duration = t, time_step = 0.05, limits = limits,
               on = [:state, :φ], measurements = "correlations" => Overlap(:φ, X)),
    ],
))

found = measure(sim.state, Overlap(:φ, X); φ = sim.φ)

# exactly, from the dense matrices of the operators on the chain
sys = System(n, Qubit())
dense = matrix(ising, sys)
ψ0 = eigvecs(Hermitian(dense))[:, 1]
u = exp(-im * t * dense)
exact = [ ψ0' * matrix(X(4), sys) * u' * matrix(X(j), sys) * u * ψ0 for j in 1:n ]

println("⟨X_4(0) X_j($t)⟩: ", round.(found; digits = 4))
println("largest error ", maximum(abs.(found - exact)))
@assert maximum(abs.(found - exact)) < 1e-6
