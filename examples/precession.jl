using TensorMixedStates

limits = Limits(
    # cutoff to keep singular values in svd
    cutoff = 1e-16,
    # maximum bond dimension
    maxdim = 100,
)

output = Output(
    # write expectation values of these observables for all sites to data file
    expect = ["X", "Y", "Z"],
    # whether to write information on the state to log () (trace, bond dimension, memory usage...)
    state_info = true,
    # C style output format for data numbers
    data_format = "%7.3f",
    # a list of checks that show the difference between an observable and an expected value depending on simulation time
    # for each check it write the simulation time, the absolute difference and both the observable measured and the expected value
    # here this will work for the first Tdvp but not the second
    checks = [X(1)=>t->cos(2t)],
)

phases(n) = [
    CreateState(
        name = "Create $n qubits pointing in the 6 directions",
        type = Pure,
        output = output,
        system = fill(site("Qubit"), n),
        state = repeat(["X+", "Y+", "Z+", "X-", "Y-", "Z-"], (n + 5) ÷ 6)[1:n],
    ),
    Tdvp(
        name = "Evolve using tdvp with hamiltonian Z on all sites",
        output = output,
        duration = 7,
        tau = 0.1,
        hamiltonian = sum(Z(i) for i in 1:n),
        limits = limits,
        data_output = output,
    ),
    CreateState(
        name = "Create $n qubits pointing in the 6 directions",
        type = Mixed,
        output = output,
        system = fill(site("Qubit"), n),
        state = repeat(["X+", "Y+", "Z+", "X-", "Y-", "Z-"], (n + 5) ÷ 6)[1:n],
    ),
    Tdvp(
        name = "Evolve using tdvp with hamiltonian Z on all sites",
        output = output,
        duration = 7,
        tau = 0.1,
        hamiltonian = sum(Z(i) for i in 1:n),
        limits = limits,
        data_output = output,
    ),
]

simdata(n) = SimData(
    name = "A simulation showing precession",
    description = "no description",
    debug = true,
    phases = phases(n),
)

runNL(simdata(10))
