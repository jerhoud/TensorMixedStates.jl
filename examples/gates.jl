using TensorMixedStates

limits = Limits(
    # cutoff to keep singular values in svd
    cutoff = 1e-16,
    # maximum bond dimension
    maxdim = 100,
)

output = Output(
    # write expectation values of these observables for all sites to data file
    expect = [X, Y, Z],
    # write correlations values for each couple of observables to data file
    correl = [(X, X), (Y, Y), (Z, Z)],
    # whether to write information on the state to log () (trace, bond dimension, memory usage...)
    state_info = true,
    # C style output format for data numbers
    data_format = "%6.2f",
)

phases(n) = [
  CreateState(
    # this gets written in the log
    name = "Building my very special state",
    # this is written when phase is over (here when state is created)
    output = output,
    # creating a pure state (other possible value: Mixed)
    type = Pure,
    # defining the size and nature of the system, here n qubits
    # mixed systems are possible as well as conserved quantum numbers
    system = fill(site("Qubit"), n),
    # define state here all "Up" (or equivalently "0", "↑" or "Z+")
    state = "Up",
  ),
  Gates(
    # this gets written in the log
    name = "Applying gate X to all qubits",
    # this is written when phase is over (here after the gates are applied)
    output = output,
    # the gates to apply
    gates = prod(X(i) for i in 1:n),
    # the limits to apply during the process
    limits = limits,
  ),
]

sim_data(n) = SimData(
    # this the name of the folder created to save simultion results
    name = "my amazing simulation with $n qubits",
    # this is written to a file so you remember what this simulation is about
    description = "A simulation with $n qubits",
    # in debug mode no files a created, it is all written to stdout
    debug = true,
    # the description of the simulation
    phases = phases(n)
)

# do it with 5 sites
runTMS(sim_data(5))
