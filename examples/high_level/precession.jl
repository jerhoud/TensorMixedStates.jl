using TensorMixedStates

limits = Limits(
    cutoff = 1e-16,
    maxdim = 100,
)

output = "data" => [ X, Y, Z, Check("check_cos", X(1), t->cos(2t))]
    

phases(n) = [
    CreateState(
        name = "Create $n qubits pointing in the 6 directions",
        final_measures = output,
        type = Pure,
        system = System(n, "Qubit"),
        state = repeat(["X+", "Y+", "Z+", "X-", "Y-", "Z-"], (n + 5) ÷ 6)[1:n],
      ),
    Evolve(
        name = "Evolve using tdvp with hamiltonian Z on all sites",
        measures = output,
        algo = Tdvp(),
        duration = 7,
        time_step = 0.1,
        evolver = -im * sum(Z(i) for i in 1:n),
        limits = limits,
    ),
]

simdata(n) = SimData(
    name = "A simulation showing precession",
    description = "no description",
    phases = phases(n),
)

runTMS(simdata(10))
