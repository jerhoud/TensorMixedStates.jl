using TensorMixedStates

limits = Limits(
    cutoff = 1e-16,
    maxdim = 100,
)

output = Output(
    observables = [X(1), Y(1), Z(1), Y(1)*Y(2), Y(1)*Z(2)*Y(3)],
    state_info = true,
)

phases(n) = [
  # create a graph_state here with a complete graph
  graph_state(Pure, complete_graph(n); output = output, limits = limits),
  # switch to mixed representation
  ToMixed(
      output = output,
      limits = limits,
  ),
  # use tdvp to evolve the system 
  Tdvp(
      #total time
      duration = 3,
      # time step (for a complete tdvp sweep)
      tau = 0.025,
      # limits to apply during the process
      limits = limits,
      # output written
      data_output = output,
      # every 4 tdvp steps
      output_periodicity = 4,
      # this the dissipator used
      dissipator = sum(DUp(i) for i in 1:n),
      # we could have add an hamiltonian = ...
  ),
 ]

sim_data(n) = SimData(
    name = "my simulation with $n qubits",
    description = "A simulation with $n qubits\nstarting from a complete graph state\nWith dissipation toward Up using tdvp",
    debug = true,
    phases = phases(n),
)

# do it for 10 qubits
runNL(sim_data(10))
