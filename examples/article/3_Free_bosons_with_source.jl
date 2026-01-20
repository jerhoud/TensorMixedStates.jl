# importing TMS and boson related identifiers
using TensorMixedStates, .Bosons

# Maximum number of bosons per lattice site:
MAXOCC = 3

# Maximum bond dimension:
MAXDIM = 100

# SVD truncation level:
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
sim_data(n, step, duration, alg, gamma) = SimData(
    # the name of the file directory created (different parameters give different names)
    name = "Boson_central_injection_$(n)_$(alg)_$gamma",
    # the content of the description file (for later reference)
    description = """
        Free boson chain with $n sites and open boundary dissipation
        Maximum boson occupancy = $MAXOCC
        cutoff = $CUTOFF
        maxdim = $MAXDIM
        algo #$alg
        time_step = $step
        duration = $duration
        boson injection rate at the center: $gamma
        empty initial state: |00...0>
        model and notations: see P L Krapivsky et al. https://doi.org/10.1088/1742-5468/ab8118
    """,
    # A vector of Phase objects describing the steps of the simulation
    phases = [
        # a Phase to create the initial state
        CreateState(
            name = "Initialization - creating an empty boson chain", # for the log file
            type = Mixed(), # type of representation
            system = System(n, Boson(MAXOCC + 1)), # system description (n boson sites)
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
                -im * sum(A(i)dag(A)(i+1) + dag(A)(i)A(i+1) for i in 1:n-1)
                + Dissipator(2 * sqrt(gamma)dag(A))(n ÷ 2),
            # save measurements
            measures = output(n),
            # every two time steps
            measures_period = 2,
        ),
    ]
)

# the parameters of the current simulation
n = 20
step = 0.1
duration = 15
#alg = 0 # TDVP
#alg = 1 # W1 order 4
alg = 2 # W2 order 4
gamma = 0.05

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
runTMS(sim_data(n, step, duration, alg, gamma); restart = true)

# we could call runTMS several times in one script with different parameters
