using TensorMixedStates

limits = Limits(
    cutoff = 1e-16,
    maxdim = 100,
)

output(n) = [
    "data" => [X, Y, Z, Y(1)Y(2), Z(1)Y(2)Y(3), EE(n ÷ 2, 4), Purity, Trace],
    "log" => "sim time",
]

sim_data(n) = SimData(
    name = "my simulation with $n qubits",
    description = """
            A simulation with $n qubits
            starting from a complete graph state
            With dissipation toward Up using tdvp""",
    phases = [
        CreateState(
            name = "Creating complete graph state",
            state = graph_state(Pure, complete_graph(n)),
            final_measures = output(n),
        ),
        ToMixed(
            final_measures = output(n),
            limits = limits,
        ),
        Evolve(
            limits = limits,
            duration = 3,
            time_step = 0.025,
            algo = Tdvp,
            evolver = sum(DUp(i) for i in 1:n),
            measures = output(n),
            measures_period = 4,
        ),
    ],
)

runTMS(sim_data(10))