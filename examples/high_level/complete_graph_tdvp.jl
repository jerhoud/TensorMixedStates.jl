using TensorMixedStates, .Qubits

limits = Limits(
    cutoff = 1e-16,
    maxdim = 100,
)

measurements(n) =
    "data" => [X, Y, Z, Y(1)Y(2), Z(1)Y(2)Y(3), EntanglementEntropy(n ÷ 2, 4), Purity, Trace]

sim_data(n) = SimData(
    name = "my_simulation_with_$(n)_qubits",
    description = """
            A simulation with $n qubits
            starting from a complete graph state
            With dissipation toward Up using tdvp""",
    phases = [
        create_graph_state(complete_graph(n); limits),
        ToMixed(
            final_measures = measurements(n),
            limits = limits,
        ),
        Evolve(
            limits = limits,
            duration = 3,
            time_step = 0.025,
            algo = Tdvp(),
            evolver = sum(Dissipator(Sp)(i) for i in 1:n),
            measures = measurements(n),
            measures_period = 4,
        ),
    ],
)

runTMS(sim_data(8))