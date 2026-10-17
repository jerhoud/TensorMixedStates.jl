# Examples

Three folders, from the shortest scripts to those of the article:

- [`low_level`](low_level): computations written directly, building a state and calling
  `ground_state`, `evolve`, `thermalize`, `apply` or `measure` on it. Each one checks its result
  against an exact value, a dense diagonalization, a free fermion solution or a closed formula,
  and is the place to start to see how a computation is written and how it is checked.
  - `spin_one_chain.jl`: the ground state and the first excitation of the spin one Heisenberg
    chain, as an excited state and as the ground state of a sector;
  - `free_fermions.jl`: the Fermi sea of a hopping chain and the correlations of its ground
    state, conserving the number of fermions;
  - `domain_wall_quench.jl`: a domain wall melting in an XX chain, with `Tdvp`, `ApproxW` and
    half the time step;
  - `driven_damped_qubits.jl`: qubits driven and jumping from up to down, relaxing to their
    steady state, found also by `steady_state`;
  - `thermal_state.jl`: the thermal state of an Ising chain, its energy and its partition
    function;
  - `noisy_ghz.jl`: a GHZ state made by a circuit and depolarized, its fidelity, the statistics
    of its parity, and a post-selection;
  - `entanglement.jl`: entropies, reduced density matrix, negativity and mutual information on
    GHZ and Werner states;
  - `strong_then_weak.jl`: fermions dephasing under a strong symmetry, then losing particles
    under a weak one.
- [`high_level`](high_level): simulations written as a list of phases run by `runTMS`, which
  write their measurements to files and can be stopped and resumed:
  - `dmrg.jl`: a ground state search stopped by its tolerance;
  - `excited_states.jl`: the ground state and the first excited state of the frustrated chain at
    its Majumdar-Ghosh point, through a named state;
  - `ising_quench.jl`: a quench of an Ising chain with periodic boundary conditions, evolved
    with `ApproxW`;
  - `precession.jl`: qubits precessing in a field, checked against the exact solution;
  - `gates.jl`: gates applied to a state with the `Gates` phase;
  - `measured_trajectory.jl`: a chain of qubits kicked and measured at each step, a quantum
    trajectory written with an `Update` phase;
  - `thermal_phase.jl`: the thermal states of an XXZ chain over a range of temperatures, written
    as JSON Lines;
  - `two_time_correlations.jl`: the correlations at two times of an Ising chain, two states
    evolved together by one `Evolve` and measured by `Overlap`;
  - `complete_graph_tdvp.jl`: a complete graph state decaying under dissipation, evolved with
    `Tdvp`;
  - `fermion_chain_conserved.jl`: a fermion chain with dephasing conserving its number of
    particles strongly, then weakly once a loss is added.
- [`article`](article): the six examples of the article describing TMS, at the sizes it
  published, hours of computation each, see its [README](article/README.md).

Run one with `julia --threads=auto examples/low_level/free_fermions.jl`. The high level ones
write a directory named after their simulation in the working directory. All the examples of
`low_level` and `high_level` are run by `test/run_examples.jl`, which the continuous integration
runs.
