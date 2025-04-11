using TensorMixedStates, .Qubits

output =
    "data" => [X, Y, Z, Trace, Purity, :sweep, :energy]

sim_data(n) = SimData(
    name = "A simple demo of Dmrg with $n qubits",
    description = "A simulation showing the use of dmrg in a simple case.", 
    phases = [
        CreateState(
            system = System(n, Qubit()),
            type = Pure(),
            randomize = 10,
        ),
        Dmrg(
            hamiltonian = sum(-Z(i) for i in 1:n) + 0.5 * sum(Z(i)Z(i + 1) for i in 1:n-1),
            limits = Limits(cutoff = 1e-10, maxdim = [5, 10, 25, 50, 100]),
            nsweeps = 10,
            measures = output,
            tolerance = 1e-8,       
        ),
    ]
)

runTMS(sim_data(10))
