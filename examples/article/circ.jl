using TensorMixedStates, .Qubits

#Maximum bond dimension:
MAXDIM = 600
#SVD truncation level:
CUTOFF = 1e-30

limits = Limits(
    cutoff = CUTOFF,
    maxdim = MAXDIM,
)

# Depolarization channel
DPL(p) = (1 - 0.75p) * Gate(Id) + 0.25p * Gate(X) + 0.25p * Gate(Y) + 0.25p * Gate(Z)
# Ising coupling gates
Rxx(ϕ) = exp(-im * ϕ * X ⊗ X)
Rzz(ϕ) = exp(-im * ϕ * Z ⊗ Z)


output(n) = [
    "OSEE.dat" => EE(n ÷ 2, 4),
    "purity.dat" => Purity,
    "trace.dat" => Trace,
    "bond_dim.dat" => Linkdim
]

sim_data(n,ϕ,steps,noise) = SimData(
    name = "Z0",
    description = """
        $n qubits
        cutoff = $CUTOFF
        maxdim = $MAXDIM
    """,
    phases = [
        CreateState(
            name = "Initialization",
            type = Mixed(),
            system = System(n, Qubit()),
            state =  "0",
            final_measures = output(n),
        ),
        [[Gates(
            name = "Applying exp(I*XX*ϕ) gates on qubits [1,2],[3,4],...",
            gates = prod(Rxx(ϕ)(2*i-1,2*i) for i in 1:n÷2),
            limits = limits,
        ),
        Gates(
            name = "Applying exp(I*ZZ*ϕ) gates on qubits [N,1], [2,3],[4,5],...",
            gates = Rzz(ϕ)(1, n)*prod(Rzz(ϕ)(2*i, 2*i+1) for i in 1:(n÷2)-1),
            limits = limits,
        ),
        Gates(
            name = "Depolarization channel on all qubits",
            final_measures = output(n),
            gates = prod(DPL(noise)(i) for i in 1:n),
            limits = limits,
        )
        ] for step in 1:steps]
    ]
)


n=20
ϕ=0.5
steps=50
noise=0.02
runTMS(sim_data(n,ϕ,steps,noise),restart=true)