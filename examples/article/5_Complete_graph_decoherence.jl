# importing TMS and qubit related identifiers
using TensorMixedStates, .Qubits

# A function depending on the simulation parameters
#   producing the SimData object describing the simulation
simdata(n, c) = SimData(
    # the name of the file directory created (different parameters give different names)
    name = "graph_$(n)_$c",
    # A vector of Phase objects describing the steps of the simulation

    phases = [
        # a call to a phase creator to build a complete graph state
        create_graph_state(complete_graph(n); limits = Limits(cutoff=1e-14)),
        # a phase to switch to mixed representation
        ToMixed(),
        # a phase for Lindblad evolution
        Evolve(
            # total evolution time
            duration = 5,
            # time step
            time_step = 0.2,
            # evolution algorithm
            algo = ApproxW(order = 1, w = 2),
            # the limits on the MPS 
            limits = Limits(maxdim=20, cutoff=1e-20),
            # evolver must be -im * hamiltonian + dissipators
            # here no Hamiltonian
            evolver = c * sum(Dissipator(Sp)(i) for i in 1:n),
            # save measurements each step
            measures = [
                "sanity.dat" => [Trace, Purity, Linkdim],
                "data.dat" => Y(1)Y(2)Z(3)
            ]
        )
    ]
)

# run the simulation
# this creates a file directory whose name was given above 
# the "restart" parameter erase the previous directory with the same name if present 
# creates a description file as specified before
# copy this script to "prog.jl"
# while running a "running" file is present
# a "stamp" file gives versions of julia, TMS and the current date
# a "log" file describing the progress of the simulation is created
# results file are created ("sanity.dat" and "data.dat")

runTMS(simdata(32, 1))

# we could call runTMS several times in one script with different parameters
