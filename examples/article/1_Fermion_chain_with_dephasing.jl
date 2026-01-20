# importing TMS and fermion related identifiers
using TensorMixedStates, .Fermions

# SVD truncation level:
CUTOFF = 1e-30

# Maximum bond dimension:
MAXDIM = 200

# A Limit object to pass to Evolve
limits = Limits(
    cutoff = CUTOFF,
    maxdim = MAXDIM,
)

# A function depending on the system size producing a vector of output specifications file => data
# each entry in the data files have the format: data_name simulation_time data1 data2 ... 
output(n) = [
    "density.dat" => N, # occupation number for all sites
    "OSEE.dat" => EE(n ÷ 2, 4), # OSEE at midpoint and 4 eignevalues 
    "purity.dat" => Purity, # purity of the density matrix
    "trace.dat" => Trace, # trace of the density matrix to check
    "bond_dim.dat" => Linkdim, # maximum bond dimension
]

# A function depending on the simulation parameters
#   producing the SimData object describing the simulation
sim_data(gamma, n, step, alg) = SimData(
    # the name of the file directory created (different parameters give different names)
    name = "Fermion_chain_$(n)_$(gamma)_$(step)_$alg",
    # the content of the description file (for later reference)
    description = """
        Spinless fermion chain (1D) quench with $n sites and open boundary conditions 
        cutoff = $CUTOFF
        maxdim = $MAXDIM
        algo #$alg
        time_step = $step
        dephasing gamma= $gamma
        alternating initial state |101010...>
        model: see https://arxiv.org/pdf/2501.07095
    """,
    # A vector of Phase objects describing the steps of the simulation
    phases = [
        # a Phase to create the initial state
        CreateState(
            name = "Initialization in state |101010...>", # for the log file
            type = Mixed(), # type of representation
            system = System(n, Fermion()), # system description (n fermion sites)
            state =  [i % 2 == 0 ? "Occ" : "Emp" for i in 1:n], # initial state
            final_measures = output(n), # save measurements for inital state
        ),
        # a phase for Lindblad evolution
        Evolve(
            # evolution algorithm TDVP (alg=0) or ApproxW 1 (alg=1) or ApproxW 2 (alg=2)
            algo= (
                (alg==0) ? Tdvp() : (
                (alg==1) ? ApproxW(order = 4, w = 1)
                         : ApproxW(order = 4, w = 2)
                )),
            # the limits on the MPS    
            limits = limits,
            # total evolution time
            duration = 4,
            # time step
            time_step = step,
            # evolver must be -im * hamiltonian + dissipators
            evolver = 
                im*sum(dag(C)(i)C(i+1)+dag(C)(i+1)C(i) for i in 1:n-1)
                + sum(Dissipator(sqrt(4 * gamma)N)(i) for i in 1:n),
            # save measurements
            measures = output(n),
            # every two time steps
            measures_period = 2,
        ),
    ]
)

# the parameters of the current simulation
n = 40
gamma = 0.75
#gamma = 1.0
#gamma = 1.25

alg = 0 # TDVP
#alg = 1 # W1 order 4
#alg = 2 # W2 order 4

time_step = 0.05


# run the simulation
# this creates a file directory whose name was given above 
# the "restart" parameter erase the previous directory with the same name if present 
# creates a description file as specified before
# copy this script to "prog.jl"
# while running a "running" file is present
# a "stamp" file gives versions of julia, TMS and the current date
# a "log" file describing the progress of the simulation is created
# results file are created ("density.dat", "OSEE.dat", "purity.dat",
#   "trace.dat", "bond_dim.dat")
runTMS(sim_data(gamma, n, time_step, alg); restart = true)

# we could call runTMS several times in one script with different parameters