# Systems and States

## Systems

```@docs
System
length(::System)
⊗(::System, ::System)
SysIndex
```

## States

```@docs
Representation
Pure
Mixed
Limits
State
State(::System, ::State)
AbstractState
length(::AbstractState)
maxlinkdim(::State)
mix
truncate(::State)
trace(::State)
trace2
norm(::State)
normalize(::State{Pure})
dag(::State{Pure})
hermitianize
hermiticity
inner
fidelity
hs_fidelity
RandomState
partial_trace
```

## Prepared states

Besides the product states the local states give, TMS builds the Slater determinants of free
fermions, and the Fermi sea of a quadratic hamiltonian, on sites of `Fermion` or `Electron`.
The graph states of qubits are built by [`graph_state`](@ref), and the thermal states by
[`thermal_state`](@ref).

```@docs
slater_state
fermi_sea
```

## Saving and loading

States can be written to disk in the HDF5 format and read back later. Several states may be
stored in the same file under different names. The same thing is available in the high level
interface with the `SaveState` and `LoadState` phases.

```@docs
save_state
load_state
```
