# importing TMS and fermion related identifiers
using TensorMixedStates, .Fermions

#Maximum bond dimension:
MAXDIM = 100

#SVD truncation level:
CUTOFF = 1e-30

# A Limit object to pass to Evolve
limits = Limits(
    cutoff = CUTOFF,
    maxdim = MAXDIM,
)

# A function depending on the system size producing a vector of output specifications file => data
# each entry in the data files have the format: data_name simulation_time data1 data2 ... 
output(n) = [
    "density.dat" => N, # occupation number for all sites
    "total_number.dat" => sum(N(i) for i in 1:n), # the total occupation number
    "OSEE.dat" => EE(n ÷ 2, 4), # OSEE at mid point and first 4 eigenvalues
    "purity.dat" => Purity, # purity of the density matrix
    "trace.dat" => Trace, # trace to observe eventuals deviation from 1
    "bond_dim.dat" => Linkdim, # maximum bond dimension
]

# A function depending on the simulation parameters
#   producing the SimData object describing the simulation
sim_data(n, step, duration, alg, Gamma) = SimData(
    # the name of the file directory created (different parameters give different names)
    name = "Fermions_with_source_$(n)_$(alg)_$Gamma",
    # the content of the description file (for later reference)
    description = """
        Free fermion chain with $n sites and particle source in the center
        cutoff = $CUTOFF
        maxdim = $MAXDIM
        algo #$alg W2o2
        time_step = $step
        duration = $duration
        fermion injection rate at the center Gamma = $Gamma
        empty initial state: |00...0>
        model and notations: see P L Krapivsky et al J. Stat. Mech. (2019) 113108 http://doi.org/10.1088/1742-5468/ab4e8e
    """,
    # A vector of Phase objects describing the steps of the simulation
    phases = [
        # a Phase to create the initial state
        CreateState(
            name = "Initialization - creating an empty (spinless) fermion chain", # for the log file
            type = Mixed(), # type of representation
            system = System(n, Fermion()), # system description (n fermion sites)
            state =  "0", # initial empty state
            final_measures = output(n), # save measurements for inital state
        ),
        # a phase for Lindblad evolution
        Evolve(
            # evolution algorithm TDVP (alg=0) or ApproxW 1 (alg=1) or ApproxW 2 (alg=2)
            # n_hermitianize = 5 makes the density matrix Hermitian again every 5 steps
            algo= (
                (alg==0) ? Tdvp() : (
                (alg==1) ? ApproxW(order = 4, w = 1, n_hermitianize = 5)
                         : ApproxW(order = 4, w = 2, n_hermitianize = 5)
                )),
            # the limits on the MPS
            limits = limits,
            # total time of evolution
            duration = duration,
            # time step
            time_step = step,
            # evolver must be -im * hamiltonian + dissipators
            evolver =
                -im * sum(dag(C)(i)*C(i+1)+dag(C)(i+1)*C(i) for i in 1:n-1)
                + Dissipator(sqrt(2*Gamma)*dag(C))(n ÷ 2), #dissipative term injecting fermions in the center of the chain
            # save measurements
            measures = output(n),
            # every two time steps
            measures_period = 2,
        ),
 
    ]
)

# the parameters of the current simulation
n = 50
step = 0.1
duration = 10
alg = 0 # TDVP
#alg = 1 # W1 order 4
#alg = 2 # W2 order 4
Gamma = 0.2

# run the simulation
# this creates a file directory whose name was given above 
# the "restart" parameter erase the previous directory with the same name if present 
# creates a description file as specified before
# copy this script to "prog.jl"
# while running a "running" file is present
# a "stamp" file gives versions of julia, TMS and the current date
# a "log" file describing the progress of the simulation is created
# results file are created ("density.dat", "total_number.dat", "OSEE.dat", "purity.dat",
#   "trace.dat", "bond_dim.dat")
runTMS(sim_data(n, step, duration, alg, Gamma); restart = true)

# we could call runTMS several times in one script with different parameters
