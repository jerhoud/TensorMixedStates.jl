# Algorithms

Each algorithm is shown at work on a simple example in the
[Algorithms](manual.md#Algorithms) section of the manual.

## Methods

An evolution or a search is given its method as `algo`: `Tdvp` or `ApproxW` for the evolutions,
`Dmrg` for the searches. Their options are their fields, and the Krylov method of their local
steps is a `Krylov`.

```@docs
TensorMixedStates.EvolutionAlgorithm
Tdvp
ApproxW
TensorMixedStates.SearchAlgorithm
Dmrg
```

## Ground states and excited states

`ground_state` searches the ground state of a Hamiltonian from a pure state, and, given the
eigenstates found before, the lowest one orthogonal to them, an excited state.

```@docs
ground_state
```

## Steady states

`steady_state` searches the steady state of a Lindbladian ``L`` by dmrg on ``L^\dagger L``,
from a mixed state.

```@docs
steady_state
```

## Time evolution

`evolve` evolves a pure state under a Hamiltonian, or a mixed state under a Lindbladian.

### Time dependent evolvers

For time evolution, it may be useful to have time dependent operators. In TMS, a time dependent
operator is a vector of terms, each an operator with its function of time, `A => f`, or an
operator alone when it does not depend on time: `[A, B => g]` is ``A + g(t) B``. For example,
to describe

```math
    h(t) = - e^{-t} \sum_{i=1}^{n-1} \sigma_x^i \sigma_x^{i+1} - \sin(t) \sum_{i=1}^n \sigma_z^i
```
we use

```@setup algorithms
using TensorMixedStates
using .Qubits
n = 6
```

```@example algorithms
hs = -im * [ -sum(X(i)X(i+1) for i in 1:n-1), -sum(Z(i) for i in 1:n)]
evolver = [ hs[1] => t -> exp(-t), hs[2] => t -> sin(t) ]
nothing # hide
```

The evolver is given to `evolve` or to an `Evolve` phase as any other evolver:

```julia
evolve(evolver, duration, initial_state; nsteps = 100, time = 0)
```

The functions are read at `time`, the time the evolution starts at, plus the time of each step;
an `Evolve` phase gives them the time of the simulation. `PreMPO(initial_state, evolver)`
prepares the evolver once for several calls, its functions kept.

`Tdvp`, and `ApproxW` of order 1 or 2, take each function at the middle of each time step,
which limits them to an error of order ``\tau^2`` in the time step ``\tau``. `ApproxW` of
order 3 or 4 keeps the error of its order: each step is made of two exponentials, each of the
evolver with its functions combined from their values at the two points of the Gauss
quadrature of the step, ``t + (1/2 \mp \sqrt{3}/6)\,\tau``, with the weights
``(3 \pm 2\sqrt{3})/12``. This is the commutator-free Magnus integrator of order 4 of Blanes
and Moan, and a step costs twice as much as with a single exponential. One of the weights is
negative: a dissipator whose rate varies by a factor of more than 13.9 between the two points
gets a negative rate in one of the two exponentials, and the step must then be shortened.

The time functions take real values. A complex one is written as its real and imaginary
parts, each with its own term, ``f(t) A = \mathrm{Re} f(t)\, A + \mathrm{Im} f(t)\, (i A)``.
A drive ``\Omega (e^{i\omega t} \sigma^+ + e^{-i\omega t} \sigma^-)`` is thus

```julia
hs = -im * [ Ω * X(1), -Ω * Y(1) ]
evolver = [ hs[1] => t -> cos(ω * t), hs[2] => t -> sin(ω * t) ]
```

```@docs
evolve
```

## Thermal states

`thermalize` takes a mixed state to ``e^{-\beta H/2} \rho \, e^{-\beta H/2}``, normalized, by
an evolution in imaginary time: from `"FullyMixed"`, the thermal state ``e^{-\beta H}/Z``.

The state can be checked without an exact reference: the energy is ``-d \log Z / d\beta``,
which the `log_trace` it returns gives apart from `expect`, a thermal state commutes with
``H``, so that a short evolution under ``-iH`` changes no measurement, and the energy falls
towards the ground state energy, which `ground_state` gives, as ``\beta`` grows. The results
must not move when the step is halved or the limits raised.

```@docs
thermalize
```

## Krylov methods

At each step of their sweeps, `Tdvp` and `Dmrg` solve a local problem by a Krylov method,
whose parameters `Krylov` holds.

```@docs
Krylov
```

## Gate application

`apply` applies an operator, gates and channels included, or an MPO, to a state.

```@docs
apply
```

## Threading

How the contractions of all these algorithms use the cores of the machine, see
[Threads and performance](@ref).

```@docs
set_threading
ThreadingState
threading_settings
```
