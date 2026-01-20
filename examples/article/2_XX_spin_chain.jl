# importing TMS and qubit related identifiers
using TensorMixedStates, .Qubits

# SVD truncation level:
CUTOFF = 1e-30

# Maximum bond dimension:
MAXDIM = 300

limits = Limits(
    cutoff = CUTOFF,
    maxdim = MAXDIM,
)

# A function depending on the system size producing a vector of output specifications file => data
# each entry in the data files have the format: data_name simulation_time data1 data2 ... 
output(n) = [
    "mag.dat" => Sz, # magnetization for all sites
    "current.dat" => 2(X(1)Y(2)-Y(1)X(2)), # measure of "current" between site 1 and 2
    "OSEE.dat" => EE(n ÷ 2, 4), # OSEE at midpoint and 4 eignevalues 
    "purity.dat" => Purity, # purity of the density matrix
    "trace.dat" => Trace, # trace of the density matrix to check
    "bond_dim.dat" => Linkdim, # maximum bond dimension
]

# A function depending on the simulation parameters
#   producing the SimData object describing the simulation
sim_data(n, step, duration, alg, eL, eR, muL, muR) = SimData(
    # the name of the file directory created (different parameters give different names)
    name = "XX_chain_$(n)_$(alg)_$(eL)_$(eR)_$(muL)_$muR",
    # the content of the description file (for later reference)
    description = """
        XX spin chain (1D) quench with $n sites and boundary dissipation
        cutoff = $CUTOFF
        maxdim = $MAXDIM
        algo #$alg
        time_step = $step
        duration = $duration
        initial state: infinite temperature product state
        Boundary bath parameters:
        epsilon_{L,R}= $eL, $eR
        mu_{L,R}= $muL, $muR
        model and notations: see https://scipost.org/10.21468/SciPostPhys.14.5.112
    """,
    # A vector of Phase objects describing the steps of the simulation
    phases = [
        # a Phase to create the initial state
        CreateState(
            name = "Initialization in the infinite temperature state (rho ~ identity)", # for the log
            type = Mixed(), # type of representation
            system = System(n, Qubit()), # system description (n qubit sites)
            state =  "FullyMixed", # initial state (infinite temperature)
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
            duration = duration,
            # time step
            time_step = step,
            # evolver must be -im * hamiltonian + dissipators
            evolver =
                -im*sum(X(i)*X(i+1)+Y(i)*Y(i+1) for i in 1:n-1)
                + Dissipator(sqrt(eL*(1+muL)*0.5)*Sp)(1)
                + Dissipator(sqrt(eL*(1-muL)*0.5)*Sm)(1)
                + Dissipator(sqrt(eR*(1+muR)*0.5)*Sp)(n)
                + Dissipator(sqrt(eR*(1-muR)*0.5)*Sm)(n),
            # save measurements
            measures = output(n),
            # every two time steps
            measures_period = 2,
        )
    ]
)

# the parameters of the current simulation
n = 30
eL = 1.0
eR = 1.0
muL = 1.0
muR = -muL
step = 0.1
duration = 40
#alg = 0 # TDVP
#alg = 1 # W1 order 4
alg = 2 # W2 order 4

# run the simulation
# this creates a file directory whose name was given above 
# the "restart" parameter erase the previous directory with the same name if present 
# creates a description file as specified before
# copy this script to "prog.jl"
# while running a "running" file is present
# a "stamp" file gives versions of julia, TMS and the current date
# a "log" file describing the progress of the simulation is created
# results file are created ("mag.dat", "current.dat", "OSEE.dat", "purity.dat",
#   "trace.dat", "bond_dim.dat")
runTMS(sim_data(n, step, duration, alg, eL, eR, muL, muR); restart = true)

# we could call runTMS several times in one script with different parameters