using TensorMixedStates, .Qubits

limits = Limits(
    cutoff = 1e-16,
    maxdim = 100,
)

measurements = "data" => [X, Y, Z, (X, X), (Y, Y), (Z, Z), Purity]

phases(n) = [
  CreateState(
    name = "Building my very special state",
    final_measures = measurements,
    type = Pure(),
    system = System(n, Qubit()),
    state = "Up",
  ),
  Gates(
    name = "Applying gate X to all qubits",
    final_measures = measurements,
    gates = prod(X(i) for i in 1:n),
    limits = limits,
  ),
]

sim_data(n) = SimData(
    name = "my_amazing_simulation_with_$(n)_qubits",
    description = "A simulation with $n qubits",
    phases = phases(n)
)

# do it with 5 sites
runTMS(sim_data(5))
