using TensorMixedStates, .Bosons
#Maximum number of bosons per lattice site:
MAXOCC = 4
#Maximum bond dimension:
MAXDIM = 100
#SVD truncation level:
CUTOFF = 1e-30

limits = Limits(
    cutoff = CUTOFF,
    maxdim = MAXDIM,
)

output(n) = [
    "density.dat" => N,
    "total_number.dat" => sum(N(i) for i in 1:n),
    "OSEE.dat" => EE(n ÷ 2, 4),
    "purity.dat" => Purity,
    "trace.dat" => Trace,
    "bond_dim.dat" => Linkdim,
]

sim_data(n,step,duration,alg,gamma) = SimData(
    name = "Q6",
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
    phases = [
        CreateState(
            name = "Initialization - creating an empty boson chain",
            type = Mixed(),
            system = System(n, Boson(MAXOCC)),
            state =  "0",
            final_measures = output(n),
        ),
        Evolve(
            algo= (
                (alg==0) ? Tdvp() : (
                (alg==1) ? ApproxW(order = 4, w = 1,n_symmetrize = 5) : ApproxW(order = 4, w = 2,n_symmetrize = 5)
                )),
            limits = limits,
            duration = duration,
            time_step = step,
            evolver =  -im*(
                    sum(A(i)*dag(A)(i+1)+dag(A)(i)*A(i+1) for i in 1:n-1)
                    )
                    +Dissipator(2*sqrt(gamma)*dag(A))(n ÷ 2),
            measures = output(n),
            measures_period = 2,
        ),
 
    ]
)

n=20
step=0.1
duration=15
#alg=0 # TDVP
#alg=1 # W1 order 4
alg=2 # W2 order 4
gamma=0.05
runTMS(sim_data(n,step,duration,alg,gamma),restart=true)