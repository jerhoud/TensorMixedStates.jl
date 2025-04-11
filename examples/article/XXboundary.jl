using TensorMixedStates, .Qubits

MAXDIM = 300

limits = Limits(
    cutoff = 1e-30,
    maxdim = MAXDIM,
)

output(n) = [
    "mag.dat" => Sz,
    "current.dat" => 2(X(1)Y(2)-Y(1)X(2)),
    "OSEE.dat" => EE(n ÷ 2, 4),
    "purity.dat" => Purity,
    "trace.dat" => Trace,
    "bond_dim.dat" => Linkdim,
]

sim_data(n,step,duration,alg,eL,eR,muL,muR) = SimData(
    name = "T3",
    description = """
        XX spin chain (1D) quench with $n sites and boundary dissipation
        cutoff = 1e-30
        maxdim = $MAXDIM
        algo #$alg
        time_step = $step
        duration = $duration
        initial state: infinite temperature product state
        Boundary bath parameters:
        epsilon_{L,R}=$eL,$eR
        mu_{L,R}=$muL,$muR
        model and notations: see https://scipost.org/10.21468/SciPostPhys.14.5.112
    """,
    phases = [
        CreateState(
            name = "Initialization in the infinite temperature state (rho ~ identity)",
            type = Mixed(),
            system = System(n, Qubit()),
            state =  "FullyMixed",
            final_measures = output(n),
        ),
        Evolve(
            algo= (
                (alg==0) ? tdvp : (
                (alg==1) ? ApproxW(order = 4, w = 1) : ApproxW(order = 4, w = 2)
                )),
            limits = limits,
            duration = duration,
            time_step = step,
            evolver =  -im*(
                    sum(X(i)*X(i+1)+Y(i)*Y(i+1) for i in 1:n-1)
                    )
                    +Dissipator(sqrt(eL*(1+muL)*0.5)*Sp)(1)
                    +Dissipator(sqrt(eL*(1-muL)*0.5)*Sm)(1)
                    +Dissipator(sqrt(eR*(1+muR)*0.5)*Sp)(n)
                    +Dissipator(sqrt(eR*(1-muR)*0.5)*Sm)(n),
            measures = output(n),
            measures_period = 2,
        )
    ]
)

n=30
eL=1.0
eR=1.0
muL=1.0
muR=-muL
step=0.1
duration=40
#alg=0 # TDVP
#alg=1 # W1 order 4
alg=2 # W2 order 4
runTMS(sim_data(n,step,duration,alg,eL,eR,muL,muR),restart=true)