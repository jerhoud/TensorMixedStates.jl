using TensorMixedStates

limits = Limits(
    cutoff = 1e-30,
    maxdim = 50,
)

output(n) = [
    "X.dat" => X,
    "Y.dat" => Y,
    "Z.dat"=> Z,
    "Z1Z2.dat" => Z(1)Z(2),
    "Z1Z3.dat" => Z(1)Z(3),
    "Z1Z4.dat" => Z(1)Z(4),
    "Z1Z5.dat" => Z(1)Z(5),
    "X1X2.dat" => X(1)X(2),
    "X1X3.dat" => X(1)X(3),
    "X1X4.dat" => X(1)X(4),
    "X1X5.dat" => X(1)X(5),
    "OSEE.dat" => EE(n ÷ 2, 4),
    "purity.dat" => [Purity, Trace, Linkdim],
    "log" => "sim time"
]

sim_data(J,h,n) = SimData(
    name = "Q1",
    description = """
        Ising chain (1D) quench with $n spins and periodic boundary conditions 
        cutoff = 1e-30
        maxdim = 100
        algo = W2 order 4
        time_step = 0.04
    """,
    phases = [
        CreateState(
            name = "Initialization in state |++...+>",
            type = Pure,
            system = System(n, "Qubit"),
            state = "X+",
            final_measures = output(n),
        ),
        ToMixed(
            final_measures=output(n),
            limits = limits,
        ), 
        Evolve(
            #algo = Tdvp(n_symmetrize = 5),
            algo =  ApproxW(order = 4, n_symmetrize = 5),
            limits = limits,
            duration = 5,
            time_step = 0.04,
            evolver =  -im*(
                    J*(sum(Z(i)*Z(i+1) for i in 1:n-1)+Z(n)*Z(1))
                    -h*sum(X(i) for i in 1:n)
                    ),
            measures = output(n),
            measures_period = 5,
        ),
 
    ]
)

n=10
J=1.0
h=-1.0
runTMS(sim_data(J,h,n),restart=true)