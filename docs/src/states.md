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

Besides the product states the local states give, TMS builds states with an exact MPS of small
bond dimension, as the GHZ, Dicke and W states, the states of singlets on pairs of sites, the
fully mixed state of a sector, from which [`thermal_state`](@ref) gives the canonical thermal
state, the Slater determinants of free fermions, and the Fermi sea
of a quadratic hamiltonian, on sites of `Fermion` or `Electron`, beside sites of other types if
need be, an impurity for instance.
The graph states of qubits are built by [`graph_state`](@ref), the AKLT state of spins one by
[`aklt_state`](@ref), and the thermal states by [`thermal_state`](@ref).

```@docs
ghz_state
dicke_state
w_state
dimer_state
fully_mixed
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
