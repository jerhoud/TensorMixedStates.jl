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
    # whether to write information on the state to log () (trace, bond dimension, memory usage...)
    state_info = true,
    # C style output format for data numbers
    data_format = "%7.3f",
)

phases(n) = [
    CreateState(
        # this gets written in the log
        name = "Creating a random state with $n qubits",
        # this is written when phase is over (here when state is created)
        output = output,
        # creating a pure state (other possible value: Mixed)
        type = Pure,
        # defining the size and nature of the system, here n qubits
        # mixed systems are possible as well as conserved quantum numbers
        system = fill(site("Qubit"), n),
        # the bond dimension of the random state created
        randomize = 3,
      ),
    Dmrg(
        hamiltonian = sum(Z(i) for i in 1:n) + 0.5 *  sum(Z(i)Z(i + 1) for i in 1:n-1),
        cutoff = 1e-10,
        maxdim = [5, 10, 25, 50, 100],
        nsweep = 10,
        output_period = 1,
        data_output = output,      
    ),
]

sim_data(n) = SimData(
    # this is the name of the folder created to save simultion results
    name = "A simple demo of Dmrg with $n qubits",
    # this is written to a file so you remember what this simulation is about
    description = "A simulation with $n qubits",
    # in debug mode no files are created, it is all written to stdout
    debug = true,
    # the description of the simulation
    phases = phases(n),
)

runTMS(sim_data(5))