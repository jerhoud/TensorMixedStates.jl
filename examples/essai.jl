using TensorMixedStates

n = 16

limits = Limits(
    cutoff = 1e-16,
    maxdim = 100,
)

output = Output(
#    expect = ["X", "Y", "Z"],
#    correl = [("X", "X"), ("Y", "Y"), ("Z", "Z")],
#    observables = [X(1), Y(3), Z(5), X(1)*Y(3), X(1)*Y(3)*Z(6)],
    observables = [X(1), Y(1), Z(1), Y(1)*Y(2), Y(1)*Z(2)*Y(3)],
    state_info = true,
    data_format = "%8.4f",
)

#= phases1 = [
    cregraph_state(Pure(n), complete_graph(n); output=output, limits=limits)...,
    ToMixed(
        output = output,
        limits = limits,
    ),
    Tdvp(
        duration = 3,
        tau = 0.01,
        limits = limits,
        data_output = output,
        data_periodicity = 10,
        dissipator = sum(DUp(i) for i in 1:n),
    ),
]


phases2 =[
    create_graph_state(Mixed(n), complete_graph(n); output=output, limits=limits)...,
    #=     ToMixed(
        output = output,
        limits = limits,
        ),
 =#]

 =#

n = 6

#= phases3 = [
    CreateState(
    
        system = Pure(fill(6, site("Qubit"))),
        state = ["X+", "X-","Y+", "Y-", "Z+", "Z-"],
        output = output,
    ),
    Gates(
        name = "applying X",
        output = output,
        limits = limits,
        gates = prod(X(i) for i in 1:6),
    ),
    ToMixed(
        output = output,
        limits = limits,
    ),
 ]
 =#

n = 5

phases4 = [
    graph_state(Pure, complete_graph(5); output = output),
    ToMixed(
        output = output,
        limits = limits,
    ), 
#=     LoadState(
        file = "coucou",
        object_name = "essai",
        output = output,
        limits = limits,
    ),
 =#
    Tdvp(
        duration = 3,
        tau = 0.01,
        limits = limits,
        data_output = output,
        output_periodicity = 10,
        dissipator = sum(DUp(i) for i in 1:5),
    ),
 ]

sim_data = SimData(
    name = "my simulation with $n qubits",
    description = "A simulation with $n qubits",
    debug = true,
    phases = phases4
)

runNL(sim_data)

