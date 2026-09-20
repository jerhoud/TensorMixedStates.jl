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

```julia
using .Qubits
```

```@docs
TensorMixedStates.Qubits
Qubit
Phase
Swap
controlled
graph_state
create_graph_state
```

## Spins

To use `Spin`, call

```julia
using .Spins
```

```@docs
TensorMixedStates.Spins
Spin
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
Sumd
```

## Defining new site types

To define a new site type, you need to define a new subtype of [`AbstractSite`](@ref) and define [`dim`](@ref) and possibly `string_state` on it (to overload do not forget to use the full name e.g. `TensorMixedStates.dim`). Then define its specific states and operators using `@def_states` and `@def_operators`.
Don't forget to define the `F` operator for fermionic sites.

### Reusing an operator name

Site types are meant to share operator names: `N` means the same thing for a `Fermion`, a
`Boson`, a `Qboson` and a `Qudit`, and each of them gives it its own matrix. Your own site
can join in, and there is nothing to do for that: name your operator `N` and
`@def_operators` will register your matrix for your site under that name.

What happens behind the scenes is that a name becomes a `const` of your module the first
time it is declared, and only then. A name that is already in scope, because you loaded a
site module exporting it or because you declared it for an earlier site of your own, is
registered for the new site and left bound as it is. So the name goes on standing for one
single operator, and the sites already using it are undisturbed.

The counterpart is that the declarations have to agree. Declaring `N` as `plain_op` when a
site already in scope declared it `selfadjoint_op` is refused, with a message saying so,
rather than quietly changing what `N` means for every site using it. If you want different
properties, you want a different name. A name already taken by something that is not an
operator at all is refused in the same way.

!!! note "Overloading `dim`"
    `dim` follows a different rule, being a function rather than a name you declare: it is
    exported by `TensorMixedStates`, so `dim(::MySite) = 2` would try to define a new
    function of your own instead of adding a method. This is why it has to be written
    `TensorMixedStates.dim(::MySite) = 2`, and the same goes for `string_state`.

```@docs
string_state
@def_states
@def_operators
@create_site_module
```
