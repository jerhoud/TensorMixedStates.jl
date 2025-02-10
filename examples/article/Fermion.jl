using TensorMixedStates
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
    "bond_dim.dat" => Linkdim
]

sim_data(n,step,duration,alg,Gamma) = SimData(
    name = "F0",
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
    phases = [
        CreateState(
            name = "Initialization - creating an empty (spinless) fermion chain",
            type = Mixed,
            system = System(n, "Fermion"),
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
                    sum(Cdag(i)*C(i+1)+Cdag(i+1)*C(i) for i in 1:n-1)
                    )
                    +sqrt(2*Gamma)*DCdag(n ÷ 2), #dissipative term injecting fermions in the center of the chain
            measures = output(n),
            measures_period = 2,
        ),
 
    ]
)

n=50
step=0.1
duration=10
alg=0 # TDVP
#alg=1 # W1 order 4
#alg=2 # W2 order 4
Gamma=0.2
runTMS(sim_data(n,step,duration,alg,Gamma),restart=true)