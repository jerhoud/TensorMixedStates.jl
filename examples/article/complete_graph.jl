using TensorMixedStates, .Qubits

dissipators(n, c) = c * sum(Dissipator(Sp)(i) for i in 1:n)

simdata(n, c) = SimData(
  name = "graph_$(n)_$c",
  phases = [
    create_graph_state(complete_graph(n); limits = Limits(cutoff=1e-14)),
    ToMixed(),
    Evolve(
      duration = 5,
      time_step = 0.2,
      algo = ApproxW(order=1, w=2),
      limits = Limits(maxdim=20, cutoff=1e-20),
      evolver = dissipators(n, c),
      measures = [
        "sanity.dat" => [Trace, Purity, Linkdim],
        "data.dat" => Y(1)Y(2)Z(3)
      ]
    )
  ]
)

runTMS(simdata(32, 1))