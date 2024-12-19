using TensorMixedStates

limits = Limits(
    cutoff = 1e-16,
    maxdim = 100,
)

output = [
    "data" => [X, Y, Z, Trace, Purity, :sweep, :energy],
    "log" => :sweep
]

sim_data(n) = SimData(
    name = "A simple demo of Dmrg with $n qubits",
    description = """
            A simulation showing the use of dmrg in a simple case.
            On pure and mixed representation.""", 
    phases = [
        CreateState(
            system = System(n, "Qubit"),
            type = Pure,
            randomize = 10,
        ),
        Dmrg(
            hamiltonian = sum(-Z(i) for i in 1:n) + 0.5 * sum(Z(i)Z(i + 1) for i in 1:n-1),
            cutoff = 1e-10,
            maxdim = [5, 10, 25, 50, 100],
            nsweeps = 10,
            measures = output,
            tolerance = 1e-8,       
        ),
        CreateState(
            system = System(n, "Qubit"),
            type = Mixed,
            randomize = 10,
        ),
        Dmrg(
            hamiltonian = sum(-Z(i) for i in 1:n) + 0.5 * sum(Z(i)Z(i + 1) for i in 1:n-1),
            cutoff = 1e-10,
            maxdim = [5, 10, 25, 50, 100],
            nsweeps = 10,
            measures = output,
            tolerance = 1e-8,       
        ),
    ]
)

runTMS(sim_data(10))
