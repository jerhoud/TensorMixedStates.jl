# Algorithms

Each algorithm is shown at work on a simple example in the
[Algorithms](manual.md#Algorithms) section of the manual.

## Ground state computation

`dmrg` computes the ground state of a hamiltonian, on a state or a simulation.

```@docs
dmrg
DmrgObserver
```

## Steady state computation

`steady_state` computes the steady state of a Lindbladian ``L`` by dmrg on ``L^\dagger L``, on a
mixed state or a simulation.

```@docs
steady_state
```

## Time evolution

`tdvp` and `approx_W` evolve a state or a simulation under a hamiltonian or a Lindbladian.

### Time dependent evolvers

For time evolution, it may be useful to have time dependent operators. In TMS, a time dependent operator
is described in the following way: a vector of indexed operators and a vector of time functions. For example,
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
nothing # hide
```

and

```@example algorithms
coefs = [ t -> exp(-t), t -> sin(t) ]
nothing # hide
```

A term that does not depend on time takes a constant function, `t -> 1.0`. Each function is
taken at the middle of each time step. `hs` is passed to `tdvp` or `approx_W` as usual, and
`coefs` as the keyword argument `coefs`.
With a `Simulation`, `t` is the simulation time; with a `State`, the evolution starts from the
time given by the keyword argument `time_start`, 0 by default:

```julia
tdvp(hs, duration, initial_state; coefs, time_start)
```

With the high level interface, the evolver of an `Evolve` phase is made time dependent by
writing

```julia
evolver = hs => coefs
```

The time functions take real values. A complex one is written as its real and imaginary
parts, each with its own term, ``f(t) A = \mathrm{Re} f(t)\, A + \mathrm{Im} f(t)\, (i A)``.
A drive ``\Omega (e^{i\omega t} \sigma^+ + e^{-i\omega t} \sigma^-)`` is thus

```julia
hs = -im * [ Ω * X(1), -Ω * Y(1) ]
coefs = [ t -> cos(ω * t), t -> sin(ω * t) ]
```

```@docs
tdvp
approx_W
TdvpObserver
ApproxWObserver
```

## Thermal states

`thermal_state` takes a mixed state to ``e^{-\beta H/2} \rho \, e^{-\beta H/2}``, normalized,
by tdvp in imaginary time: from `"FullyMixed"`, the thermal state ``e^{-\beta H}/Z``.

The state can be checked without an exact reference: the energy is ``-d \log Z / d\beta``,
which the `log_trace` it returns gives apart from `expect`, a thermal state commutes with
``H``, so that a short evolution under ``-iH`` changes no measurement, and the energy falls
towards the ground state energy, which `dmrg` gives, as ``\beta`` grows. The results must
not move when the step is halved or the limits raised.

```@docs
thermal_state
ThermalObserver
```

## Krylov methods

At each step of their sweeps, `tdvp`, `thermal_state`, `dmrg` and `steady_state` solve a local
problem by a Krylov method, whose parameters `Krylov` holds.

```@docs
Krylov
```

## Gate application

`apply` applies gates, or an MPO, to a state or a simulation.

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
