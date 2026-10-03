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

## Conserving a quantity

A site can be told that a quantity is conserved. The indices it draws then carry that charge,
the tensors become block sparse, and a contraction only pairs blocks whose charges agree,
which makes a large computation smaller and faster, see [What it saves](@ref). In exchange a
state is confined to the sector it was built in.

The quantity is named by one of the site's own operators, which has to be diagonal with
eigenvalues that are either all integers or all roots of unity:

```julia
Fermion(conserve = N)                  # the number of particles
Qubit(conserve = 2Sz)                  # twice the magnetisation, so that it is an integer
Electron(conserve = (Ntot, 2Sz))       # two quantities at once
Boson(4, conserve = parity(N))         # only the parity of the number, a charge modulo 2
```

Half integer quantities are written doubled, `2Sz` rather than `Sz`, so that the charges they
give are integers.

Sites that declare nothing may sit in a system beside sites that do. They take a trivial
charge and keep their whole vocabulary, so a `Qubit()` next to conserving fermions still
accepts `X` and the state `"+"`.

### What it forbids

Conserving is a promise about the whole computation, and what breaks it is refused with a
message naming the culprit rather than discovered in the middle of a run.

A **state** lives in one sector. `Fermion(conserve = N)` takes `"Occ"` and `"Emp"` and refuses
`[1, 1] / √2`, which superposes two numbers of particles and so has no number of its own.

An **operator** must carry a definite charge, its flux, which is the difference between the
charges of the states it connects. [`flux`](@ref) gives it:

```jldoctest
julia> using TensorMixedStates, .Fermions

julia> flux(N, Fermion(conserve = N))
QN("N",0)

julia> flux(dag(C), Fermion(conserve = N))
QN("N",1)

julia> flux(C, Fermion(conserve = N))
QN("N",-1)
```

`X` connects the two states of a qubit in both directions at once, so under
`Qubit(conserve = N)` it has no flux and is refused. Under `Qubit(conserve = parity(N))` it
has one, `QN("parity(N)",1,2)`, the two differences becoming the same one modulo 2.

The flux is asked of every factor of one site, and not only of the whole operator: each
factor becomes a tensor of its own. `Sx(1) * Sx(2) + Sy(1) * Sy(2)` conserves the
magnetisation as a whole, but under `Spin(1, conserve = 2Sz)` its factors have none, and it is
refused. Written with factors of a definite flux, as `(Sp(1) * Sm(2) + Sm(1) * Sp(2)) / 2`, it
is accepted.

### Weak and strong symmetries

For a **mixed** representation there are two ways of conserving a quantity, and they are not
the same promise.

By default the symmetry is **weak**: what is asked is that the density matrix commute with the
charge. Every jump operator of definite charge preserves that, particle loss and gain
included, and the state may spread over several sectors, as a thermal state does.

[`strong`](@ref) asks more: that every jump operator commute with the charge. The charge of
the ket and that of the bra are then conserved separately, which cuts the blocks finer, and in
exchange the state lives in a single sector, exactly as a pure one does.

|  | `conserve = N` | `conserve = strong(N)` |
|---|---|---|
| jump operators | any of definite charge, `C`, `dag(C)`, `N` | only those of zero flux, `N`, `dag(C) * C` |
| states | may mix sectors | one sector only |
| blocks of the mixed index of a fermion | 3 | 4 |
| `partial_trace` | yes | no, what is left spreads over sectors; [`weaken`](@ref) first |

Dephasing, whose jump operator is `N` itself, is the usual strong case; particle loss, whose
jump is `C`, is not. Declaring `strong` and then using a jump that moves the charge is refused
by a message saying which quantity the jump changes and pointing at the weak form.

A strong quantity takes two of the four charge components ITensors allows, where a weak one
takes a single one, so at most two quantities can be declared strong.

### Moving between the levels

A quantity is thus conserved at one of three levels: strongly, weakly, or not at all.
[`symmetries`](@ref) tells which, in the form a declaration takes, and [`weaken`](@ref) takes a
state, a system or a simulation down to a lower level, building the system it lands on:

```jldoctest
julia> using TensorMixedStates, .Electrons

julia> mysystem = System(4, Electron(conserve = (strong(Ntot), 2Sz)));

julia> mystate = State{Mixed}(mysystem, ["Up", "Dn", "Up", "Dn"]);

julia> symmetries(mysystem)
(2Sz, strong(Ntot))

julia> symmetries(weaken(mystate).system)          # one level down: every strong quantity made weak
(2Sz, Ntot)

julia> symmetries(weaken(weaken(mystate)).system)  # or, when none is strong, every quantity dropped
()

julia> symmetries(weaken(mystate, 2Sz).system)     # exactly these
2Sz

julia> weaken(mystate, symmetries(mysystem)) === mystate   # the identity
true
```

This is a step of a simulation in its own right. A phase may evolve under a strong symmetry,
which dephasing allows, and the next one continue under a weak one, where particle loss
becomes possible. There is no way back: the finer blocks of a strong symmetry cannot be
recovered from the coarser ones, and a target asking for more than the state has is refused.

Weakening a state costs one local tensor per site, and it is also how a strongly conserving
state gets a `partial_trace`.

### What it saves

Ten steps of `Tdvp` on the fermion chain with dephasing of the example
`examples/high_level/fermion_chain_conserved.jl`, measured on one machine:

|  | no charge | `conserve = N` | `conserve = strong(N)` |
|---|---|---|---|
| 8 sites, `maxdim = 32` | 21.7 s | 23.0 s | 40.5 s |
| 12 sites, `maxdim = 64` | 208.9 s | 80.8 s | 121.8 s |

At 8 sites the charges bring nothing; at 12 the weak form is 2.6 times faster than no charge
at all. The strong form is the slower of the two here: what it offers is physical rather than
speed, a state confined to one sector and a jump moving the charge refused.

```@docs
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
conserve_string
show(::IO, ::AbstractSite)
@def_states
@def_operators
@create_site_module
```
