using TensorMixedStates
import ITensors: op, @OpName_str, @SiteType_str
using ITensors

#Maximum bond dimension:
MAXDIM = 600
#SVD truncation level:
CUTOFF = 1e-30

limits = Limits(
    cutoff = CUTOFF,
    maxdim = MAXDIM,
)

#Depolarization channel
DPL(i; p) = create_mixed_gate("DP", [I, X, Y, Z], [1-0.75*p, p*0.25, p*0.25, p*0.25],i; p)


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
            type = Mixed,
            system = System(n, "Qubit"),
            state =  "0",
            final_measures = output(n),
        ),
        [[Gates(
            name = "Applying exp(I*XX*ϕ) gates on qubits [1,2],[3,4],...",
            gates = prod(Rxx(2*i-1,2*i;ϕ) for i in 1:n÷2),
            limits = limits,
        ),
        Gates(
            name = "Applying exp(I*ZZ*ϕ) gates on qubits [N,1], [2,3],[4,5],...",
            gates = Rzz(1,n;ϕ)*prod(Rzz(2*i,2*i+1;ϕ) for i in 1:(n÷2)-1),
            limits = limits,
        ),
        Gates(
            name = "Depolarization channel on all qubits",
            final_measures = output(n),
            gates = prod(DPL(i;p=noise) for i in 1:n),
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