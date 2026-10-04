# Extending TMS

TMS can be extended without touching its code, from a script or from a package of your
own:

- [site types](@ref "Site types of one's own"), with their states and operators;
- operators and measurements, see [`named`](@ref) and [`StateFunc`](@ref);
- [phases](@ref own-phases) of a simulation;
- [algorithms](@ref "Algorithms of one's own") of time evolution;
- [representations](@ref "Representations of one's own") of a state, with a state type of
  their own.

## Site types of one's own

To define a new site type, you need to define a new subtype of [`AbstractSite`](@ref) and define [`dim`](@ref) and possibly [`string_state`](@ref) on it (to overload do not forget to use the full name e.g. `TensorMixedStates.dim`). Then define its specific states and operators using [`@def_states`](@ref) and [`@def_operators`](@ref).
Don't forget to define the [`F`](@ref) operator for fermionic sites.

An atom with a ground state, an excited one and a Rydberg one, driven from one to the next and
decaying from the Rydberg state, with an interaction between Rydberg atoms side by side:

```@example extending
using TensorMixedStates

struct Atom <: AbstractSite end          # a ground state g, an excited one e, a Rydberg one r
TensorMixedStates.dim(::Atom) = 3        # the full name, to add a method to dim

@def_states(Atom(), [
    "G" => [1., 0., 0.],
    "E" => [0., 1., 0.],
    "R" => [0., 0., 1.],
])

@def_operators(Atom(), [
    plain_op => [
        Sge = [0. 1. 0. ; 0. 0. 0. ; 0. 0. 0.],     # |g><e|
        Ser = [0. 0. 0. ; 0. 0. 1. ; 0. 0. 0.],     # |e><r|
    ],
    selfadjoint_op => [
        Nr = [0. 0. 0. ; 0. 0. 0. ; 0. 0. 1.],      # the Rydberg population
    ],
])

hamiltonian = sum(Sge(i) + dag(Sge)(i) + Ser(i) + dag(Ser)(i) for i in 1:4) +
              sum(Nr(i) * Nr(i + 1) for i in 1:3)
myrho = State{Mixed}(System(4, Atom()), "G")
mystate = tdvp(-im * hamiltonian + sum(Dissipator(Ser)(i) for i in 1:4), 0.5, myrho;
               nsweeps = 5)
measure(mystate, [Nr, Trace])
```

The operators declared are constants of the module the declaration is in. Declared in a module
of your own, the site type and the operator names are exported for the programs that use it.

### Conserved quantities

A site type may carry a field named `conserve`, of type `String`, recording the quantities it
conserves. Declare it only if your site can have some: a site type that cannot simply leaves
it out and conserves nothing, which is why the field is optional rather than part of the
interface.

You do not fill it by hand. Give your site a keyword constructor in the manner of the ones of
this package, which reads the operators the user names and records the charges they give on
each basis state:

```julia
struct MySite <: AbstractSite
    conserve::String
end

MySite(; conserve = ()) = MySite(conserve_string(MySite(""), conserve))
```

A conserved quantity is an operator of your site, diagonal, whose eigenvalues are either all
integers or all roots of unity. `MySite(conserve = N)` and `MySite(conserve = (N, 2Sz))` are
then written the same way as for the site types of this package.

### The fields of a site

A site is rebuilt from the values of its fields, in their order, by `MySite(values...)`: when a
state file, which records them, is read, and when a `Weaken` phase changes what the site
conserves. Julia gives this constructor to every struct, unless an inner constructor replaces
it, in which case one taking all the fields in their order has to be kept. For the states on
your site to be saved, its fields must be numbers, booleans, symbols, strings or `nothing`.

### Free fermions

A site holding fermions that are free under a quadratic hamiltonian gives their annihilation
operators, one per species, by a method of `TensorMixedStates.fermion_species`, as `(C,)` for
a `Fermion` and `(Cup, Cdn)` for an `Electron`. `slater_state` and `fermi_sea` then build its
Slater determinants and Fermi seas, its local state `"0"` being the one with no fermion.

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

Which name to share is a question of how many there are to count. A site with a single
number of particles calls it `N`, as `Fermion`, `Boson`, `Qboson`, `Qudit`, `Qubit` and
`Spin` all do; a site holding several species names them apart, as `Electron` and `Tj` do
with `Nup`, `Ndn` and `Ntot`, where a bare `N` would leave the reader guessing which one it
meant.

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

## Operators and measurements of one's own

An operator of your own is defined with [`named`](@ref), from a matrix, a function of the
sites or an expression, or with [`Operator`](@ref) for one of several sites. Its type, which
`simplify` reasons with, see [`OpType`](@ref), is read off its matrix: here an involution.

```@example extending
using TensorMixedStates, .Qubits

myop = named([1 1 ; 1 -1] / √2, "MyOp")
myop.type
```

A type given to `named` with a matrix is checked against it at once, as far as it can be
without a site, and any type is checked against the matrix of the operator when that matrix is
laid on a site. `simplify` relies on the type before, squaring an involution to the identity
for instance, so a wrong type given to a function or an expression gives a wrong result.

A matrix of several sites given alone does not say how many sites it acts on, which is then
given in braces, with the type:

```@example extending
myswap = Operator{2}("MySwap", [1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1], involution_op)
myswap(4, 7)
```

Such an operator, or a function of an operator of several sites, as
`exp(-im * t * (X⊗X + Y⊗Y) / 4)`, can only be applied as a gate. To measure it or to put it
in a hamiltonian, give the sites it acts on when creating it, one per index or a single one
for identical sites, whose number is then read off the size of the matrix:

```@example extending
myswap2 = named([1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1], "MySwap2", Qubit())
```

It is then split into a sum of products of one site operators, which becomes its definition,
the way `Swap` is defined by an expression:

```@example extending
myswap2.expr
```

The factors carry a definite charge of what their sites conserve. On a fermionic site the
matrix has to commute with `F`, since it is taken as it is, with no Jordan-Wigner string:
an operator moving fermions between sites is written with `C` and `dag(C)` instead.

A measurement of your own is a [`StateFunc`](@ref), a function of the state, or a
[`TimeFunc`](@ref), a function of the simulation time.

## [Phases of one's own](@id own-phases)

A phase of your own is a subtype of [`AbstractPhase`](@ref) with the three fields every phase has,
`name`, `time_start` and `final_measures`, and a method of [`TensorMixedStates.run_phase`](@ref)
for it, which returns the simulation the phase leaves behind. The full name is needed,
`run_phase` not being exported. The phases of the library are written the same way, each in a
file of `src/phases`, and make as many examples.
Written as a loop with [`run_steps`](@ref), the phase is checkpointed, stopped and resumed between two
steps, as those of the library are between two sweeps:

```julia
using TensorMixedStates, .Qubits

Base.@kwdef struct Kicks <: AbstractPhase
    name::String = "Kicks"
    time_start = nothing
    final_measures = []
    nkicks::Int
    measures = []
end

TensorMixedStates.run_phase(sim::Simulation, p::Kicks) =
    run_steps(sim, p.nkicks) do sim, k
        sim = apply(exp(-0.3im * X)(1), sim)                # the kick
        sim = Simulation(sim, sim.state, sim.time + 0.1)    # and the time it takes
        output(sim, p.measures; sweep = k)                  # k for the :sweep measurement
        return sim
    end

runTMS(SimData(name = "kicks", phases = [
    CreateState{Pure}(2, Qubit(), "Up"),
    Kicks(nkicks = 10, measures = "data" => Z(1)),
]))
```

Within the method, `output` measures the simulation, `log_msg` writes to its log and
`get_sim_file` gives a file of the simulation to write anything else to, which is cut back on
a resume as the others are. A phase driving a solver is described with
[`TensorMixedStates.run_phase`](@ref).

A phase creating the state, the adapter of another library for instance, says so by a method
of [`TensorMixedStates.creates_state`](@ref), which lets it be the first phase of a simulation,
and gives the system it creates the state on by one of
[`TensorMixedStates.phase_system`](@ref), on which a resumed run puts the state of its
checkpoint back. Within a phase, [`stopped`](@ref) tells whether the run is stopping,
[`resume_time`](@ref) the time a resumed phase goes on from and [`committed_time`](@ref) the time
of the last step committed, which a phase stopped in its course has reached, and
`save_state(file, name, sim)` saves the state, refusing a file of the simulation, as its
checkpoint.

A value the steps carry from one to the next, a sum, a count or the state of a random number
generator, is given to `run_steps` as `carry`: each step receives it and returns it updated
with the simulation, and a resumed run gets it back, see [`run_steps`](@ref). A phase evolving
one step at a time prepares its evolver once with [`PreMPO`](@ref), which the solvers take in
place of the operator.

## Algorithms of one's own

An `Evolve` phase hands the evolution of its state to [`TensorMixedStates.evolve`](@ref), whose method
is chosen by the type of its algorithm and by that of the state. An algorithm of your own is a
subtype of [`Algo`](@ref) with a method for the states it evolves, which returns the simulation it
leaves behind. The method is called before the phase has read its resume point, so that a loop
written with `run_steps` is checkpointed, stopped and resumed between two steps:

```julia
using TensorMixedStates, .Qubits

struct Stepwise <: Algo end     # tdvp, one step at a time

function TensorMixedStates.evolve(::Stepwise, ::State, sim::Simulation, phase::Evolve;
                                  evolver, coefs, nsteps, kwargs...)
    dt = phase.duration / nsteps
    return run_steps(sim, nsteps) do sim, k
        sim = tdvp(evolver, dt, sim; coefs, phase.limits)
        if mod(k, phase.measures_period) == 0
            output(sim, phase.measures; sweep = k)
        end
        return sim
    end
end

runTMS(SimData(name = "stepwise", phases = [
    CreateState{Pure}(2, Qubit(), "Up"),
    Evolve(duration = 1., time_step = 0.1, algo = Stepwise(), evolver = -im * X(1),
           measures = "data" => Z(1)),
]))
```

`evolver` is the evolver of the phase, and `coefs` the functions of time of a time dependent
one, `nothing` otherwise, which `tdvp` takes as they are. The method ends its keywords with
`kwargs...`, see [Keywords of the methods you give](@ref extending-keywords).

## Representations of one's own

A representation of your own is a subtype of [`Representation`](@ref), which `CreateState`
takes as its `type`, and needs nothing more. Its states are of a type of your own, a subtype of
[`AbstractState`](@ref) with a field `system` holding the `System` of their physical sites.
TMS reads it to save a state, to choose the threading, to log what a `Weaken` phase changes,
and to put the reference of `Fidelity` or `Overlap` on the system of the measured state.

That system need not be the one the representation keeps its tensors on: a mixed state held as
a pure state on a doubled system, a purification, has the system of its physical sites in
`system`, and the doubled one in a field of its own.

The phases and the measurements reach these states through the methods you give:

- `TensorMixedStates.run_phase(sim, phase::CreateState{MyRepresentation})` creates the state
  from the fields of the phase;
- `expect(state, op::IndexedOp)`, `expect1(state, ops)` and `expect2(state, pairs)` compute
  what `measure`, and so `output`, asks for: the expectation value of an operator placed on
  sites, normalised by the trace, those of one site operators on every site, and the
  correlations of pairs of them. The operator is given as the user wrote it, on the physical
  sites, its factors in any order and without the Jordan-Wigner strings of fermionic
  operators, which [`map_sites`](@ref) places on another system. `measure` gives all the
  operators it measures at once to `expect(state, ops)`, which calls the method of one
  operator on each, unless you give one for `ops::Vector{<:IndexedOp}` that computes them
  together;
- `apply(gates, state; limits, kwargs...)` applies the gates of a `Gates` phase, whose channels
  [`kraus_operators`](@ref) gives as their Kraus operators;
- `TensorMixedStates.evolve(algo, state, sim, phase; evolver, coefs, nsteps, kwargs...)`
  evolves the state
  in an `Evolve` phase, for each algorithm it supports, the hamiltonian and the jump operators
  of the evolver being given by [`lindblad_terms`](@ref);
- `TensorMixedStates.write_state(group, state)` and
  `TensorMixedStates.read_state(::Type{MyState}, group, sites, system)` save the state and read
  it back, for `SaveState`, `LoadState` and the checkpoints: the first writes its tensors in the
  HDF5 group, where its type and its sites are already written, and the second rebuilds it from
  them, on `system` when that is not `nothing`. The type is written by its name alone: a type
  with parameters writes them in the group as well, and is read back by a method for the type
  without them, `read_state(::Type{<:MyState}, group, sites, system)`. A state holding a
  `State`, on the physical sites or on a system of its own, writes it in a subgroup,
  `write_state(create_group(group, "inner"), inner)`, `create_group` coming from HDF5, and
  reads it back with `read_state(State{Pure}, group["inner"], sites, nothing)`, given the sites
  it lies on. Whatever else the state needs to go on as it would have, the state of a random
  number generator for instance, is written in the group as numbers or arrays of them.

The other phases go through `truncate`, `mix`, `partial_trace`, `weaken`, `dmrg` and
`steady_state`, and the state functions of the measurements through `trace`, `norm`,
`entanglement_entropy` and the like. A state gets what its type has a method for, and a phase or
a measurement it does not support raises a `MethodError`. `LoadState` goes through `truncate`
only when it is given limits.

## [Keywords of the methods you give](@id extending-keywords)

A method you give for a function TMS calls with keywords, `evolve`, `apply` or
`truncate(state; limits)`, ends its keywords with `kwargs...`, as `Stepwise` above does, even
when it uses every keyword it is given. A later version of TMS may pass one more keyword, the
generator of random numbers of a trajectory for instance, and Julia refuses a keyword that a
method does not declare: a method declaring only those it uses would then fail, where one
taking `kwargs...` ignores what it does not need.
