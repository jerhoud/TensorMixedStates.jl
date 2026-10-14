# Sites

## General

```@docs
AbstractSite
dim(::AbstractSite)
```

The state `"FullyMixed"` represents the infinite temperature mixed state, that is a density matrix proportional to the identity matrix.

```@docs
Id
F
```

There are eight predefined site types `Qubit`, `Qudit`, `Spin`, `Boson`, `Fermion`, `Electron`, `Tj` and `Qboson`.

## Qubit

To use `Qubit`, call

```julia
using .Qubits
```

```@docs
TensorMixedStates.Qubits
Qubit
Phase
graph_state
thermal_relaxation_gate
thermal_relaxation_dissipator
```

## Spins

To use `Spin`, call

```julia
using .Spins
```

```@docs
TensorMixedStates.Spins
Spin
aklt_state
```

## Boson

To use `Boson`, call

```julia
using .Bosons
```

```@docs
TensorMixedStates.Bosons
Boson
```

## Fermion

To use `Fermion`, call

```julia
using .Fermions
```

```@docs
TensorMixedStates.Fermions
Fermion
```

## Electron

To use `Electron`, call

```julia
using .Electrons
```

```@docs
TensorMixedStates.Electrons
Electron
```

## Tj

To use `Tj`, call

```julia
using .Tjs
```

```@docs
TensorMixedStates.Tjs
Tj
```

## Qboson

To use `Qboson`, call

```julia
using .Qbosons
```

```@docs
TensorMixedStates.Qbosons
Qboson
```

## Qudit

To use `Qudit`, call

```julia
using .Qudits
```

```@docs
TensorMixedStates.Qudits
Qudit
```

## Conservation

What conserving a quantity means and costs is explained in [Conserving a quantity](@ref).

```@docs
Conservation
conserving(::Type{<:AbstractSite}, ::Vararg{Any})
strong
flux
symmetries
weaken
```

## Defining new site types

How to define a site type of your own, its states and its operators is explained in
[Site types of one's own](@ref).

```@docs
string_state
fermion_parities(::AbstractSite)
fermion_species(::AbstractSite)
show(::IO, ::AbstractSite)
@def_states
@def_operators
@create_site_module
```
