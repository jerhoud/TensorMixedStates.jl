# Sites

## General

```@docs
AbstractSite
dim(::AbstractSite)
Index(::AbstractSite)
state(::AbstractSite, ::String)
```

The state `"FullyMixed"` represents the infinite temperature mixed state, that is a density matrix proportional to the identity matrix.

```@docs
Id
F
```

There are six predefined site types `Qubit`, `Spin`, `Boson`, `Fermion`, `Electron` and `Tj`.

## Qubit

To use `Qubit`, call

    using .Qubits

```@docs
Qubit
Phase
controlled
graph_state
create_graph_state
```

## Spins

To use `Spin`, call

    using .Spins

```@docs
Spin
```

## Boson

To use `Boson`, call

    using .Bosons

```@docs
Boson
```

## Fermion

To use `Fermion`, call

    using .Fermions

```@docs
Fermion
```

## Electron

To use `Electron`, call

    using .Electrons

```@docs
Electron
```

## Tj

To use `Tj`, call

    using .Tjs

```@docs
Tj
```

## Defining new site types

To define a new site type, you need to define a new subtype of [`AbstractSite`](@ref) and define [`dim`](@ref) and possibly `generic_state` on it (to overload do not forget to use the full name e.g. `TensorMixedStates.dim`). Then define its specific states and operators using `@def_states` and `@def_operators`.
Don't forget to define the `F` operator for fermionic sites.

```@docs
generic_state
@def_states
@def_operators
```
