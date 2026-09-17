# Sites

## General

```@docs
AbstractSite
dim(::AbstractSite)
Index(::AbstractSite)
state(::AbstractSite, ::String)
identity_operator
```

The state `"FullyMixed"` represents the infinite temperature mixed state, that is a density matrix proportional to the identity matrix.

```@docs
Id
F
```

There are eight predefined site types `Qubit`, `Qudit`, `Spin`, `Boson`, `Fermion`, `Electron`, `Tj` and `Qboson`.

## Qubit

To use `Qubit`, call

    using .Qubits

```@docs
TensorMixedStates.Qubits
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
TensorMixedStates.Spins
Spin
```

## Boson

To use `Boson`, call

    using .Bosons

```@docs
TensorMixedStates.Bosons
Boson
```

## Fermion

To use `Fermion`, call

    using .Fermions

```@docs
TensorMixedStates.Fermions
Fermion
```

## Electron

To use `Electron`, call

    using .Electrons

```@docs
TensorMixedStates.Electrons
Electron
```

## Tj

To use `Tj`, call

    using .Tjs

```@docs
TensorMixedStates.Tjs
Tj
```

## Qboson

To use `Qboson`, call

    using .Qbosons

```@docs
TensorMixedStates.Qbosons
Qboson
```

## Qudit

To use `Qudit`, call

    using .Qudits

```@docs
TensorMixedStates.Qudits
Qudit
```

## Defining new site types

To define a new site type, you need to define a new subtype of [`AbstractSite`](@ref) and define [`dim`](@ref) and possibly `string_state` on it (to overload do not forget to use the full name e.g. `TensorMixedStates.dim`). Then define its specific states and operators using `@def_states` and `@def_operators`.
Don't forget to define the `F` operator for fermionic sites.

```@docs
string_state
@def_states
@def_operators
@create_site_module
```
