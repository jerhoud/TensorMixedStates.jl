using TensorMixedStates, .Qubits

measurements =
    "data" => [X, Y, Z, Trace, Purity, :step, :energy]

sim_data(n) = SimData(
    name = "A simple demo of Dmrg with $n qubits",
    description = "A simulation showing the use of dmrg in a simple case.",
    # the random state the search starts from, drawn again the same by every run
    seed = 1,
    phases = [
        CreateState(
            system = System(n, Qubit()),
            type = Pure,
            randomize = 10,
        ),
        GroundState(
            hamiltonian = sum(-Z(i) for i in 1:n) + 0.5 * sum(Z(i)Z(i + 1) for i in 1:n-1),
            limits = Limits(cutoff = 1e-10, maxdim = [5, 10, 25, 50, 100]),
            algo = Dmrg(nsweeps = 6, tol = 1e-8),
            measurements = measurements,
        ),
    ]
)

runTMS(sim_data(8))
