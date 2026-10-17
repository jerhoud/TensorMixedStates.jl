# Home

TensorMixedStates (TMS) is a Julia library to make simulations of closed or open quantum
systems using Matrix Product States representations.

Version 2.0 is a rewrite of the library: a program written for 1.x needs the edits listed in
[Upgrading from 1.x](@ref).

## Features

- Pure states and, for open systems, density matrices, both represented as matrix product
  states.
- A rich set of sites and operators, easily extended by the user, with a very expressive
  syntax for observables, gates, channels, Hamiltonians and Lindbladians.
- Ground states and excited states with DMRG, Hamiltonian and Lindbladian evolution with TDVP
  and with the WI and WII approximations, steady states, and thermal states by evolution in
  imaginary time.
- Prepared states: GHZ, Dicke and W states, graph states, singlets on pairs of sites, the AKLT
  state, Slater determinants and Fermi seas of free fermions, and the fully mixed states of a
  sector.
- Gates and channels, noisy gates included, and the measurement of a state with its
  probabilities, its samples and its collapse.
- Measurements of expectation values, correlations, entanglement entropies, purities, fidelities
  and reduced density matrices.
- Conserved quantities, with weak and strong symmetries for open systems.
- An optional high level interface, which writes a simulation in a few lines, with its
  measurements, its output files and checkpoints, from which an interrupted simulation resumes.

Being based on ITensor, TMS delivers high performance computations, and uses the cores of
the machine, see [Threads and performance](@ref).

## Installation

To use TMS, you need to have Julia installed on your system. Installing Julia is usually easy and fast, see [The Julia Programming Language](https://julialang.org/) for instructions. TMS requires at least Julia version 1.10.5 to run.

To install TMS in Julia, launch the Julia interface (by typing `julia` on the command line) and type

```
]add TensorMixedStates
```

The `]` switches the prompt to the package manager of Julia.

After downloading TMS, Julia compiles it, with no interaction on the user's part, the first time
and after each update of TMS or of a package it uses. This compilation runs small computations,
so that the first calls of a script do not wait for their own compilation, and takes five to six
minutes and up to 7 GB of memory. On a machine with less memory, a login node of a cluster for
instance, these computations are left out by a file `LocalPreferences.toml` next to the
`Project.toml` of the environment, written before installing TMS or before the next
`using TensorMixedStates`:

```toml
[TensorMixedStates]
precompile_workload = false
```

The compilation then takes about ten seconds and less than 1 GB, but each new session compiles
what it uses as it runs: the first simulation of [High Level Interface](@ref) takes about three
minutes rather than six seconds.

## A first example

The ground state of a transverse field Ising chain of ten qubits, and a few measurements on
it:

```julia
using TensorMixedStates, .Qubits

n = 10
mysystem = System(n, Qubit(conserve = parity(N)))
hamiltonian = -sum(X(i)X(i + 1) for i in 1:n-1) - 0.5 * sum(Z(i) for i in 1:n)
energy, ground = ground_state(hamiltonian, State{Pure}(mysystem, "Up");
                              algo = Dmrg(nsweeps = 10), limits = Limits(maxdim = 20))
measure(ground, [X(1)X(n), Z])
```

The hamiltonian flips the spins two by two, so it conserves the parity of the number of down
spins, and the qubits are told so: the state stays in the sector it starts in, that of the
ground state. `energy` is about `-9.7655`, and `measure` gives the correlation between the two
ends of the chain, about `0.748`, and the magnetization along `z` on every site. The
[Manual](@ref) explains each of these steps.

## Using TMS

To use TMS in your code you need to write

```julia
using TensorMixedStates
```

Julia script names are usually written with a .jl extension. Once you have written your script, you can execute it with

```sh
julia --threads=auto my_script.jl
```

or `julia -t auto my_script.jl` for short, which starts Julia with as many threads as the
machine has; `julia my_script.jl` runs the script too, on a single thread of Julia. What the
threads bring, and how TMS uses them, is in [Threads and performance](@ref).

You can also use TMS in the interactive Julia interpreter.

## Documentation

You are currently reading it!

You can have access to inline documentation on TMS at the Julia prompt simply by typing "?" followed by the function name or type name you are interested in. For example

```
?runTMS
```

For a function TMS shares with ITensorMPS, as `expect` or `sample`, the help shows that of
ITensorMPS first, then that of TMS.

For this to work you must have first imported TMS with

```julia
using TensorMixedStates
```

### For AI assistants

An AI asked to compute with TMS finds the library on one page in [Quick reference](@ref),
shipped with the package as `docs/src/quickref.md`, the map of the documentation in
[`llms.txt`](https://jerhoud.github.io/TensorMixedStates.jl/stable/llms.txt), and a skill for
coding agents, its workflow, conventions and traps, in the folder
[`skills/tensormixedstates`](https://github.com/jerhoud/TensorMixedStates.jl/tree/main/skills/tensormixedstates)
of the repository.

### Learning about matrix product states

For an introduction to matrix product states and tensor networks, see the lecture notes of the
course Grégoire Misguich gave at the 9th Les Houches summer school on Computational Physics:
Open Quantum Systems, in June 2026, available on [arXiv:2606.24803](https://arxiv.org/abs/2606.24803).
Their [repository](https://github.com/gregoire-misguich/Introduction-to-matrix-product-states-and-tensor-networks)
holds the sixteen Julia examples the notes use, some of them written with TMS.

## Examples

Working examples are in the folder
[`examples`](https://github.com/jerhoud/TensorMixedStates.jl/tree/main/examples) of the
repository, in three folders.

The folder [`low_level`](https://github.com/jerhoud/TensorMixedStates.jl/tree/main/examples/low_level)
holds computations written directly, each checking its result against an exact value:

- `spin_one_chain.jl`: the ground state and the first excitation of the spin one Heisenberg
  chain, as an excited state and as the ground state of a sector;
- `free_fermions.jl`: the Fermi sea of a hopping chain and the correlations of its ground state,
  conserving the number of fermions;
- `domain_wall_quench.jl`: a domain wall melting in an XX chain, with `Tdvp`, `ApproxW` and half
  the time step;
- `driven_damped_qubits.jl`: qubits driven and jumping from up to down, relaxing to their
  steady state, found also by `steady_state`;
- `thermal_state.jl`: the thermal state of an Ising chain, its energy and its partition function;
- `noisy_ghz.jl`: a GHZ state made by a circuit and depolarized, its fidelity, the statistics of
  its parity, and a post-selection;
- `entanglement.jl`: entropies, reduced density matrix, negativity and mutual information on GHZ
  and Werner states;
- `strong_then_weak.jl`: fermions dephasing under a strong symmetry, then losing particles under
  a weak one.

The folder [`high_level`](https://github.com/jerhoud/TensorMixedStates.jl/tree/main/examples/high_level)
holds simulations written with the high level interface:

- `dmrg.jl`: a ground state search with `GroundState`, stopped by its tolerance;
- `excited_states.jl`: the ground state and the first excited state of the frustrated chain at
  its Majumdar-Ghosh point, through a named state;
- `ising_quench.jl`: a quench of an Ising chain with periodic boundary conditions, evolved
  with `ApproxW`;
- `precession.jl`: qubits precessing in a field, checked against the exact solution with a
  `Check`;
- `gates.jl`: gates applied to a state with the `Gates` phase;
- `measured_trajectory.jl`: a chain of qubits kicked and measured at each step, a quantum
  trajectory written with an `Update` phase;
- `thermal_phase.jl`: the thermal states of an XXZ chain over a range of temperatures, written as
  JSON Lines;
- `two_time_correlations.jl`: the correlations at two times of an Ising chain, two states
  evolved together by one `Evolve` and measured by `Overlap`;
- `complete_graph_tdvp.jl`: a complete graph state decaying under dissipation, evolved with
  `Tdvp`;
- `fermion_chain_conserved.jl`: a fermion chain with dephasing conserving its number of
  particles strongly, then weakly once a loss is added.

The folder [`article`](https://github.com/jerhoud/TensorMixedStates.jl/tree/main/examples/article)
holds the six examples of the reference article:

- `1_Fermion_chain_with_dephasing.jl`: a spinless fermion chain with dephasing, from an
  alternating state;
- `2_XX_spin_chain.jl`: an XX spin chain with dissipation at its two ends, from the infinite
  temperature state;
- `3_Free_bosons_with_source.jl` and `4_Free_fermions_with_source.jl`: free bosons, or free
  fermions, injected at the centre of an empty chain;
- `5_Complete_graph_decoherence.jl`: the graph state of a complete graph decaying under
  dissipation;
- `6_Brickwall.jl`: a brickwall circuit of `Rxx` and `Rzz` gates with depolarizing noise.

## References

TMS is described in the following article, published in SciPost Physics Codebases. Please
cite both the article and the codebase release it documents, release 1.28, which is version
1.2.8 of TMS, as SciPost asks:

Jérôme Houdayer and Grégoire Misguich, *TensorMixedStates: A Julia library for simulating
pure and mixed quantum states using matrix product states*,
[SciPost Phys. Codebases **72** (2026)](https://doi.org/10.21468/SciPostPhysCodeb.72).

Jérôme Houdayer and Grégoire Misguich, *Codebase release 1.28 for TensorMixedStates*,
[SciPost Phys. Codebases **72-r1.28** (2026)](https://doi.org/10.21468/SciPostPhysCodeb.72-r1.28).

In BibTeX:

```bibtex
@article{TensorMixedStates,
  title     = {{TensorMixedStates}: A {Julia} library for simulating pure and mixed
               quantum states using matrix product states},
  author    = {Houdayer, J{\'e}r{\^o}me and Misguich, Gr{\'e}goire},
  journal   = {SciPost Phys. Codebases},
  pages     = {72},
  year      = {2026},
  publisher = {SciPost},
  doi       = {10.21468/SciPostPhysCodeb.72},
}

@article{TensorMixedStatesRelease,
  title     = {{Codebase} release 1.28 for {TensorMixedStates}},
  author    = {Houdayer, J{\'e}r{\^o}me and Misguich, Gr{\'e}goire},
  journal   = {SciPost Phys. Codebases},
  pages     = {72-r1.28},
  year      = {2026},
  publisher = {SciPost},
  doi       = {10.21468/SciPostPhysCodeb.72-r1.28},
}
```

## Acknowledgements and feedback

If TMS helped produce the data of a publication, a word about it in your acknowledgements
would be warmly appreciated. It costs you a single line, and it is what makes it possible
to keep the library developed and maintained.

We would also be delighted to hear from you directly, and you should not hesitate for a
moment. What are you doing with TMS? What works well, what gets in your way, what do you
need that is not there yet? The contact address is the one given in the article. Such
messages are read with care, and they are what shapes what comes next.

## Bug reports and feature requests

TMS is still in development and certainly contains bugs. If you think you have found one please report it on the
[GitHub page](https://github.com/jerhoud/TensorMixedStates.jl).

Feature requests may also be sent on the [GitHub page](https://github.com/jerhoud/TensorMixedStates.jl).

## Module

```@docs
TensorMixedStates
```

## Index

```@index
```
