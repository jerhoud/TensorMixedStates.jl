# Systems and States

## Systems

```@docs
System
length(::System)
⊗(::System, ::System)
siteind(::Type{Pure}, ::System, ::Int64)
```

## States

```@docs
TensorMixedStates.Representation
Pure
Mixed
Limits
State
TensorMixedStates.AbstractState
length(::TensorMixedStates.AbstractState)
maxlinkdim(::State)
mix
truncate(::State)
trace(::State{Pure})
tr(::TensorMixedStates.AbstractState)
purity(::State{Pure})
norm(::State)
normalize(::State{Pure})
dag(::State{Pure})
hermitianize
hermiticity
inner
fidelity
hs_fidelity
random_state
partial_trace
```

## Prepared states

Besides the product states the local states give, TMS builds states with an exact MPS of small
bond dimension, as the superpositions and the mixtures of product states, the GHZ, Dicke and W
states, the states of singlets on pairs of sites, the
fully mixed state of a sector, from which [`thermalize`](@ref) gives the canonical thermal
state, the Slater determinants of free fermions, and the Fermi sea
of a quadratic hamiltonian, on sites of `Fermion` or `Electron`, beside sites of other types if
need be, an impurity for instance. Any other state of which one knows the tensors of an MPS is
built by `mps_state`, and a state small enough to be written down, by its dense vector or
density matrix, with `dense_state`, whose vector `Vector` gives back.
The graph states of qubits are built by [`graph_state`](@ref), the AKLT state of spins one by
[`aklt_state`](@ref), and the thermal states by [`thermalize`](@ref).

```@docs
mps_state
dense_state
Vector(::State{Pure})
superposition
mixture
ghz_state
dicke_state
w_state
dimer_state
fully_mixed
slater_state
fermi_sea
```

## Products and permutations

Two states give the state of their two systems together, `a ⊗ b` on `a.system ⊗ b.system`, the
sites of `b` after those of `a`, and `permute_sites` moves the sites of a state, which
interleaves two states made apart: a system and its ancillas, two copies of a state, the two
chains of a ladder. An operator of the second state is placed on the product by
`map_sites(k -> k + length(a), op)`, and on the permuted state by `map_sites(k -> p[k], op)`.

```@docs
⊗(::State{R}, ::State{R}) where R
permute_sites
```

## Saving and loading

States can be written to disk in the HDF5 format and read back later. Several states may be
stored in the same file under different names. The same thing is available in the high level
interface with the `SaveState` and `LoadState` phases.

```@docs
save_state
load_state
```
