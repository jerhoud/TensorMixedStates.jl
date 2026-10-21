using TensorMixedStates, .Qubits

limits = Limits(
    cutoff = 1e-16,
    maxdim = 100,
)

# the first qubit, starting along x, precesses as cos(2t), which the check compares with
measurements = "data" => [X, Y, Z,
    "check_cos" => c -> Check(expect(c.state, X(1)), cos(2c.time); tol = 1e-8)]


phases(n) = [
    CreateState(
        name = "Create $n qubits pointing in the 6 directions",
        final_measurements = measurements,
        type = Pure,
        system = System(n, Qubit()),
        state = repeat(["X+", "Y+", "Z+", "X-", "Y-", "Z-"], (n + 5) ÷ 6)[1:n],
      ),
    Evolve(
        name = "Evolve using tdvp with hamiltonian Z on all sites",
        measurements = measurements,
        algo = Tdvp(),
        duration = 7,
        time_step = 0.1,
        evolver = -im * sum(Z(i) for i in 1:n),
        limits = limits,
    ),
]

simdata(n) = SimData(
    name = "precession",
    description = "no description",
    phases = phases(n),
)

runTMS(simdata(8))
