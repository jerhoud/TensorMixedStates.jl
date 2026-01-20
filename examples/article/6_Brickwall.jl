# importing TMS and qubits related identifiers
using TensorMixedStates, .Qubits

#Maximum bond dimension:
MAXDIM = 600
#SVD truncation level:
CUTOFF = 1e-30

# A Limit object to pass to Gates
limits = Limits(
    cutoff = CUTOFF,
    maxdim = MAXDIM,
)

# Depolarization channel
DPL(p) = (1 - 0.75p) * Gate(Id) + 0.25p * Gate(X) + 0.25p * Gate(Y) + 0.25p * Gate(Z)
# Ising coupling gates
Rxx(ϕ) = exp(-im * ϕ * X ⊗ X)
Rzz(ϕ) = exp(-im * ϕ * Z ⊗ Z)


# A function depending on the system size producing a vector of output specifications file => data
# each entry in the data file have the format: data_name simulation_time data1 data2 ... 
output(n) = [
    "OSEE.dat" => EE(n ÷ 2, 4), # OSEE at mid point and first 4 eigenvalues
    "purity.dat" => Purity, # purity of the density matrix
    "trace.dat" => Trace, # trace to observe eventual deviations from 1
    "bond_dim.dat" => Linkdim # maximum bon dimension
]

# A function depending on the simulation parameters
#   producing the SimData object describing the simulation
sim_data(n, ϕ, steps, noise) = SimData(
    # the name of the file directory created (different parameters give different names)
    name = "Brickwall_$(n)_$(ϕ)_$(steps)_$noise",
    # the content of the description file (for later reference)
    description = """
        $n qubits
        cutoff = $CUTOFF
        maxdim = $MAXDIM
    """,
    # A vector of Phase objects describing the steps of the simulation
    # subvectors of phases are allowed
    phases = [
        # a Phase to create the initial state
        CreateState(
            name = "Initialization", # for the log file
            type = Mixed(), # type of representation
            system = System(n, Qubit()), # system description (n qubits)
            state =  "0",  # initial all 0 state
            final_measures = output(n), # save measurements for inital state
        ),
        # A vector of vector of Gates phases. Subvectors to any depth are ok in phases
        # one Gates phase has a name for the log, gates to apply, and limits on the MPS to enforce
        # note that `prod` is a standard Julia function not specific to TMS
        [[Gates(
            name = "Applying Rxx(ϕ) gates on qubits [1,2], [3,4],...",
            gates = prod(Rxx(ϕ)(2*i-1, 2*i) for i in 1:n÷2),
            limits = limits,
        ),
        Gates(
            name = "Applying Rzz(ϕ) gates on qubits [N,1], [2,3], [4,5],...",
            gates = Rzz(ϕ)(1, n)*prod(Rzz(ϕ)(2*i, 2*i+1) for i in 1:(n÷2)-1),
            limits = limits,
        ),
        Gates(
            name = "Depolarization channel on all qubits",
            final_measures = output(n),
            gates = prod(DPL(noise)(i) for i in 1:n),
            limits = limits,
        )
        ] for step in 1:steps]
    ]
)

# parameters of this simulation
n = 20
ϕ = 0.5
steps = 50
noise = 0.02

# run the simulation
# this creates a file directory whose name was given above 
# the "restart" parameter erase the previous directory with the same name if present 
# creates a description file as specified before
# copy this script to "prog.jl"
# while running a "running" file is present
# a "stamp" file gives versions of julia, TMS and the current date
# a "log" file describing the progress of the simulation is created
# results file are created ("OSEE.dat", "purity.dat", "trace.dat", "bond_dim.dat")
runTMS(sim_data(n, ϕ, steps, noise); restart = true)

# we could call runTMS several times in one script with different parameters
