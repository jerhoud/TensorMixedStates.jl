using TensorMixedStates, .Qubits

limits = Limits(
    cutoff = 1e-16,
    maxdim = 100,
)

output = "data" => [X, Y, Z, (X, X), (Y, Y), (Z, Z), Purity]

phases(n) = [
  CreateState(
    name = "Building my very special state",
    final_measures = output,
    type = Pure(),
    system = System(n, Qubit()),
    state = "Up",
  ),
  Gates(
    name = "Applying gate X to all qubits",
    final_measures = output,
    gates = prod(X(i) for i in 1:n),
    limits = limits,
  ),
]

sim_data(n) = SimData(
    name = "my amazing simulation with $n qubits",
    description = "A simulation with $n qubits",
    phases = phases(n)
)

# do it with 5 sites
runTMS(sim_data(5))
