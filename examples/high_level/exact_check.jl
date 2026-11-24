# A program checked on exact states: an Ising chain quenched with damping, written once as a
# SimData and run with backend = ExactState, which holds the density matrix whole and evolves it
# without truncation nor error of the time step, then on MPS of growing link dimensions. Against
# the exact run, the MPS runs show the error of the truncation, which vanishes once the link
# dimension holds the whole state, 4^(n/2) for a density matrix: Tdvp then projects on the whole
# space, and its steps are exact for an evolver that does not depend on time. Checked so on a
# small system, the same program runs at the sizes no exact computation reaches.

using TensorMixedStates, .Qubits

n = 6
ham = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - 0.7 * sum(X(i) for i in 1:n)
damping = sum(Dissipator(sqrt(0.1) * Sm)(i) for i in 1:n)

sim_data(backend, maxdim) = SimData(
    name = "exact_check_$(nameof(backend))_$maxdim",
    description = "an Ising chain quenched with damping, backend $backend, maxdim $maxdim",
    backend = backend,
    phases = [
        CreateState{Mixed}(System(n, Qubit()), "+"),
        Evolve(evolver = -im * ham + damping, duration = 2.0, time_step = 0.05,
               limits = Limits(maxdim = maxdim),
               measurements = Data("z") => [Z(1), Z(n ÷ 2), "ZZ" => Z(1) * Z(n), Purity]),
    ],
)

largest_difference(rows, exact) = maximum(abs(row[c] - ref[c]) for (row, ref) in zip(rows, exact)
                                          for c in ("Z(1)", "Z($(n ÷ 2))", "ZZ", "Purity"))

exact = runTMS(sim_data(ExactState, 4^(n ÷ 2))).data["z"]
errors = Dict(maxdim => largest_difference(runTMS(sim_data(State, maxdim)).data["z"], exact)
              for maxdim in (4, 8, 16, 4^(n ÷ 2)))
for (maxdim, e) in sort(collect(errors))
    println("maxdim = $maxdim: largest difference with the exact run ", e)
end
@assert errors[4^(n ÷ 2)] < 1e-8
