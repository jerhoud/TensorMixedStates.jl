# Home

TensorMixedStates (TMS) is a Julia library to make simulations of closed or open quantum
systems using Matrix Product States representations.

## Features

- Pure states and, for open systems, density matrices, both represented as matrix product
  states.
- A rich set of sites and operators, easily extended by the user, with a very expressive
  syntax for observables, gates, Hamiltonians and Lindbladians.
- Ground states with DMRG, Hamiltonian and Lindbladian evolution with TDVP and with the WI and
  WII approximations, and steady states.
- Gates, noisy gates included.
- Conserved quantities, with weak and strong symmetries for open systems.
- An optional high level interface, which writes a simulation in a few lines, with its
  measurements, its output files and checkpoints.

Being based on ITensor, TMS delivers high performance computations and naturally runs in
parallel.

## Installation

To use TMS, you need to have Julia installed on your system. Installing Julia is usually easy and fast, see [The Julia Programming Language](https://julialang.org/) for instructions. TMS requires at least Julia version 1.10.5 to run.

To install TMS in Julia, launch the Julia interface (by typing `julia` on the command line) and type

```
]add TensorMixedStates
```

The `]` switches the prompt to the package manager of Julia.

After downloading TMS, Julia will automatically compile and install it. This process usually takes a couple of minutes and requires no interaction on the user's part.

### BLAS backend

Most of the running time of TMS is spent in the tensor contractions of ITensors, which
themselves call BLAS. Julia ships with OpenBLAS and TMS uses it as it comes: switching the
BLAS backend affects the whole Julia session, so that choice is left to you rather than
made by a library you load.

On `x86_64` machines running Linux or Windows, Intel's MKL is often noticeably faster on
these contractions. To use it, add `MKL` to your project and load it before TMS:

```julia
using MKL
using TensorMixedStates
```

MKL is not distributed for macOS nor for ARM machines, where OpenBLAS is the only option.

## A first example

The ground state of a transverse field Ising chain of ten qubits, and a few measurements on
it:

```julia
using TensorMixedStates, .Qubits

n = 10
mysystem = System(n, Qubit(conserve = parity(N)))
hamiltonian = -sum(X(i)X(i + 1) for i in 1:n-1) - 0.5 * sum(Z(i) for i in 1:n)
energy, ground = dmrg(hamiltonian, State{Pure}(mysystem, "Up"); nsweeps = 10,
                      limits = Limits(maxdim = 20))
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
machine has. A script calling the functions of TMS directly, rather than through `runTMS`,
should then start with `set_threading(:dense)`: see [Threads and performance](@ref) for what
these settings bring. `julia my_script.jl` runs the script too, on a single thread of Julia.

You can also use TMS in the interactive Julia interpreter.

## Documentation

You are currently reading it!

You can have access to inline documentation on TMS at the Julia prompt simply by typing "?" followed by the function name or type name you are interested in. For example

```
?runTMS
```

For this to work you must have first imported TMS with

```julia
using TensorMixedStates
```

### Learning about matrix product states

For an introduction to matrix product states and tensor networks, see the lecture notes of the
course Grégoire Misguich gave at the 9th Les Houches summer school on Computational Physics:
Open Quantum Systems, in June 2026, available on [arXiv:2606.24803](https://arxiv.org/abs/2606.24803).
Their [repository](https://github.com/gregoire-misguich/Introduction-to-matrix-product-states-and-tensor-networks)
holds the sixteen Julia examples the notes use, some of them written with TMS.

## Examples

Working examples are in the folder
[`examples`](https://github.com/jerhoud/TensorMixedStates.jl/tree/main/examples) of the
repository, all of them written with the high level interface.

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

The folder [`high_level`](https://github.com/jerhoud/TensorMixedStates.jl/tree/main/examples/high_level)
holds shorter ones:

- `dmrg.jl`: a ground state search with `GroundState`, stopped by its tolerance;
- `ising_quench.jl`: a quench of an Ising chain with periodic boundary conditions, evolved
  with `ApproxW`;
- `precession.jl`: qubits precessing in a field, checked against the exact solution with a
  `Check`;
- `gates.jl`: gates applied to a state with the `Gates` phase;
- `complete_graph_tdvp.jl`: a complete graph state decaying under dissipation, evolved with
  `Tdvp`;
- `fermion_chain_conserved.jl`: a fermion chain with dephasing conserving its number of
  particles strongly, then weakly once a loss is added.

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
