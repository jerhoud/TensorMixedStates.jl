using TensorMixedStates, .Fermions

MAXDIM = 200

limits = Limits(
    cutoff = 1e-30,
    maxdim = MAXDIM,
)

output(n) = [
    "density.dat" => N,
    "OSEE.dat" => EE(n ÷ 2, 4),
    "purity.dat" => Purity,
    "trace.dat" => Trace,
    "bond_dim.dat" => Linkdim,
]

sim_data(gamma,n,step,alg) = SimData(
    name = "T6b",
    description = """
        Spinless fermionchain (1D) quench with $n sites and open boundary conditions 
        cutoff = 1e-30
        maxdim = $MAXDIM
        algo #$alg
        time_step = $step
        dephasing gamma=$gamma
        alternating initial state |101010...>
        model: see https://arxiv.org/pdf/2501.07095
    """,
    phases = [
        CreateState(
            name = "Initialization in state |101010...>",
            type = Pure(),
            system = System(n, Fermion()),
            state =  [i % 2 == 0 ? "Occ" : "Emp" for i in 1:n],
            final_measures = output(n),
        ),
        ToMixed(
            final_measures=output(n),
            limits = limits,
        ), 
        Evolve(
            algo= (
                (alg==0) ? Tdvp() : (
                (alg==1) ? ApproxW(order = 4, w = 1) : ApproxW(order = 4, w = 2)
                )),
            limits = limits,
            duration = 4,
            time_step = step,
            evolver =  -im*(
                    -sum(dag(C)(i)*C(i+1)+dag(C)(i+1)*C(i) for i in 1:n-1)
                    )
                    +sum(Dissipator(sqrt(4*gamma)*N)(i) for i in 1:n),
            measures = output(n),
            measures_period = 2,
        ),
 
    ]
)

n=40
gamma=0.75
#gamma=1.0
#gamma=1.25

alg=0 # TDVP
#alg=1 # W1 order 4
#alg=2 # W2 order 4

time_step=0.05
runTMS(sim_data(gamma,n,time_step,alg),restart=true)