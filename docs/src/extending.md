# Extending TMS

TMS can be extended without touching its code, from a script or from a package of your
own:

- [site types](@ref "Site types of one's own"), with their states and operators;
- operators and measurements, see [`named`](@ref), [`Operator`](@ref) and the functions of the
  context that [`measure`](@ref) takes;
- [phases](@ref own-phases) of a simulation;
- [methods](@ref "Methods of one's own") of evolution and of search;
- [destinations](@ref "Destinations of one's own") of the measurements;
- [states](@ref "States of one's own") of a representation of their own.

!!! warning "Experimental"
    The interfaces of the [phases](@ref own-phases), the
    [methods](@ref "Methods of one's own"), the
    [destinations](@ref "Destinations of one's own") and the
    [states](@ref "States of one's own") of one's own, and the functions they name, are
    experimental. They were written for extensions still to come, quantum trajectories,
    purifications or adapters of other libraries, and may change in a version 2.x when the
    first of these shows what they lack. The functions they name are not exported: a program
    adds a method to a function under its full name, as `TensorMixedStates.run_phase`, and
    imports a type it uses, as with `using TensorMixedStates: AbstractPhase`. The site types,
    the operators and the measurements of one's own are not concerned.

## Site types of one's own

A site type is a subtype of [`AbstractSite`](@ref) that defines [`dim`](@ref), declares its
states and its operators with [`@def_states`](@ref) and [`@def_operators`](@ref), and, if it is
fermionic, gives the parity of its basis states by [`fermion_parities`](@ref), from which its
`F` is derived. A method is added to a function of TMS under its full name,
`TensorMixedStates.dim`.

An atom with a ground state, an excited one and a Rydberg one, driven from one to the next and
decaying from the Rydberg state, with an interaction between Rydberg atoms side by side:

```@example extending
using TensorMixedStates

struct Atom <: AbstractSite end          # a ground state g, an excited one e, a Rydberg one r
TensorMixedStates.dim(::Atom) = 3        # the full name, to add a method to dim

@def_states Atom [
    "G" => [1., 0., 0.],
    "E" => [0., 1., 0.],
    "R" => [0., 0., 1.],
]

@def_operators Atom [
    Sge = [0. 1. 0. ; 0. 0. 0. ; 0. 0. 0.],     # |g><e|
    Ser = [0. 0. 0. ; 0. 0. 1. ; 0. 0. 0.],     # |e><r|
    Nr = [0. 0. 0. ; 0. 0. 0. ; 0. 0. 1.],      # the Rydberg population
]

hamiltonian = sum(Sge(i) + dag(Sge)(i) + Ser(i) + dag(Ser)(i) for i in 1:4) +
              sum(Nr(i) * Nr(i + 1) for i in 1:3)
myrho = State{Mixed}(System(4, Atom()), "G")
mystate = evolve(-im * hamiltonian + sum(Dissipator(Ser)(i) for i in 1:4), 0.5, myrho;
                 nsteps = 5)
measure(mystate, [Nr, Trace])
```

The names of the operators become constants of the module the declaration is in, the first time
they are declared, which a package exports with its site type. The help of a site type, its
docstring, lists its operators in a section `# Operators`, as the site types of this package do,
which the help of an operator, `?Nr`, gathers from every site type naming it:

```
# Operators

- `Sge, Ser` : the transitions from g to e and from e to r
- `Nr`       : the Rydberg population
```

A site type declaring `Sp` and `Sz` also has `Sm`, `Sx`, `Sy` and `S2`, one declaring `A` has
`Q` and `P`, and one declaring `Cup` and `Cdn` has `Nup`, `Ndn` and `Ntot`, as the help of
[`AbstractSite`](@ref) lists, unless it declares these names itself.

### Conserved quantities

A site type that can conserve quantities records them in a field `conserve` of type
[`Conservation`](@ref), which its keyword constructor fills with [`conserving`](@ref), reading
the charges the operators the user names give on each basis state:

```julia
struct MySite <: AbstractSite
    conserve::Conservation
end

MySite(; conserve = ()) = conserving(MySite; conserve)
```

A conserved quantity is an operator of your site, diagonal with integer eigenvalues, or such
an operator modulo an integer, `mod(Q, m)` or `parity(Q)`. `MySite(conserve = N)` and
`MySite(conserve = (strong(N), 2Sz))` are then written as for the site types of this package. A
site type that cannot conserve anything leaves the field out.

### The fields of a site

A site is rebuilt from the values of its fields, in their order, by `MySite(values...)`: when a
state file, which records them, is read, and when `weaken` changes what the site conserves.
Julia gives this constructor to every struct, unless an inner constructor replaces it, in which
case one taking all the fields in their order has to be kept.

### Free fermions

A site holding fermions that are free under a quadratic hamiltonian gives their annihilation
operators, one per species, by a method of `TensorMixedStates.fermion_species`, as `(C,)` for
a `Fermion` and `(Cup, Cdn)` for an `Electron`. `slater_state` and `fermi_sea` then build its
Slater determinants and Fermi seas, its local state `"0"` being the one with no fermion.

### Reusing an operator name

Site types are meant to share operator names: `N` means the same thing for a `Fermion`, a
`Boson`, a `Qboson` and a `Qudit`, and each of them gives it its own matrix. Your own site
can join in, and there is nothing to do for that: name your operator `N` and
`@def_operators` declares your matrix for your site under that name.

A name becomes a `const` of your module the first time it is declared, and only then. A name
that is already in scope, because you loaded a site module exporting it or because you
declared it for an earlier site of your own, is declared for the new site and left bound as it
is, so that it goes on standing for one single operator. A name already standing for something
else, an operator with a definition of its own or anything that is not an operator, is
refused. A declaration is a method of the site type: declaring the same name again for the same
site type replaces the first declaration, with the warning Julia gives when a method is
overwritten in a script, and an error when a package is precompiled.

Which name to share is a question of how many there are to count. A site with a single
number of particles calls it `N`, as `Fermion`, `Boson`, `Qboson`, `Qudit`, `Qubit` and
`Spin` all do; a site holding several species names them apart, as `Electron` and `Tj` do
with `Nup`, `Ndn` and `Ntot`, where a bare `N` would leave the reader guessing which one it
meant.

!!! note "Overloading `dim`"
    `dim` follows a different rule, being a function rather than a name you declare: it is
    exported by `TensorMixedStates`, so `dim(::MySite) = 2` would try to define a new
    function of your own instead of adding a method. This is why it has to be written
    `TensorMixedStates.dim(::MySite) = 2`, and the same goes for `string_state` and
    `fermion_parities`.

## Operators and measurements of one's own

An operator of your own is defined with [`named`](@ref), from a matrix, a function of the site
or an expression of operators, or with [`Operator`](@ref), which takes the number of its sites:

```@example extending
using TensorMixedStates, .Qubits

myop = named([1 1 ; 1 -1] / √2, "MyOp")
myswap = Operator{2}("MySwap", [1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1])
expect(State{Pure}(System(2, Qubit()), ["Up", "+"]), myswap(1, 2) * myop(1))
```

Placed on sites, it is measured, applied or put in an evolver as any other. A matrix of several
sites is read in the basis of the sites as it is: on fermionic sites, it carries no
Jordan-Wigner string, and an operator moving fermions between sites is written with `C` and
`dag(C)` instead.

A measurement of your own is a function of the context, `"name" => c -> f(c.state)`, see
[`measure`](@ref), which may read the time, `c.time`, and the values of the simulation.

## [Phases of one's own](@id own-phases)

A phase of your own is a subtype of [`TensorMixedStates.AbstractPhase`](@ref) with the three
fields every phase has, `name`, `time_start` and `final_measurements`. A phase made of steps
gives [`TensorMixedStates.prepare_phase`](@ref), what its steps need, `step_count`, their
number, and `run_step`, which makes one step and returns the simulation it leaves: the run
measures, checkpoints, stops and resumes it between two steps, as it does the phases of the
library. A field `measurements`, with `measurements_period`, is measured after the steps:

```julia
using TensorMixedStates, .Qubits
using TensorMixedStates: AbstractPhase

Base.@kwdef struct Kicks <: AbstractPhase
    name::String = "Kicks"
    time_start = nothing
    final_measurements = []
    nkicks::Int
    measurements = []
end

TensorMixedStates.prepare_phase(p::Kicks, _) = p
TensorMixedStates.step_count(p::Kicks) = p.nkicks
TensorMixedStates.run_step(::Kicks, sim, _) =
    Simulation(sim, apply(exp(-0.3im * X)(1), sim.state), sim.time + 0.1)

runTMS(SimData(name = "kicks", phases = [
    CreateState{Pure}(2, Qubit(), "Up"),
    Kicks(nkicks = 10, measurements = "data" => Z(1)),
]))
```

`Simulation(sim, state, time; values...)` gives the simulation another state, and possibly
another time and values, and `sim.rng` is the generator of random numbers of the run, which
its checkpoints keep. A phase of a single call gives instead
[`TensorMixedStates.run_phase`](@ref), called with the number of steps a resumed run has already
done, and calls [`TensorMixedStates.step_done`](@ref) at the end of each of its steps, which
measures, checkpoints and stops the run there.

A phase providing values to the measurements, as `:energy`, names them by a method of
[`TensorMixedStates.provided_values`](@ref). Within a phase, `output` measures the simulation
and `log_message` writes to its log. A function changing the simulation in a few steps, without
a type of its own, is an [`Update`](@ref) phase.

## Methods of one's own

An evolution goes through [`TensorMixedStates.evolve_step`](@ref), whose method is chosen by the
type of its `algo`: a method of evolution of your own is a subtype of
[`TensorMixedStates.EvolutionAlgorithm`](@ref) with a method of `evolve_step`, which makes one
step, and which `evolve`, `thermalize` and their phases then take as they take `Tdvp`. Two steps
of `Tdvp` of half the length make one:

```julia
using TensorMixedStates, .Qubits
using TensorMixedStates: EvolutionAlgorithm

struct HalfSteps <: EvolutionAlgorithm end

function TensorMixedStates.evolve_step(::HalfSteps, pre, coefs, state, time, dt, step, limits)
    half = TensorMixedStates.evolve_step(Tdvp(), pre, coefs, state, time, dt / 2, step, limits)
    return TensorMixedStates.evolve_step(Tdvp(), pre, coefs, half, time + dt / 2, dt / 2, step,
                                         limits)
end

evolve(-im * sum(X(i) for i in 1:4), 0.5, State{Pure}(System(4, Qubit()), "Up");
       algo = HalfSteps(), nsteps = 5)
```

A search goes through [`TensorMixedStates.search`](@ref) in the same way: a method of search of
your own is a subtype of [`TensorMixedStates.SearchAlgorithm`](@ref) with a method of `search`,
which `ground_state`, `steady_state` and their phases then take as they take `Dmrg`.

## Destinations of one's own

A destination of the measurements is a subtype of [`TensorMixedStates.Destination`](@ref),
compared by value, which `open` opens once and which is then given each row of values by
[`TensorMixedStates.write_row`](@ref), and closed at the end of the run. A destination printing
the time and the values on a line:

```julia
using TensorMixedStates, .Qubits
using TensorMixedStates: Destination

struct Line <: Destination end
struct OpenLine end

Base.open(::Line) = OpenLine()
Base.close(::OpenLine) = nothing
TensorMixedStates.write_row(::OpenLine, time, row) = println(time, ": ", join(last.(row), ", "))

runTMS(SimData(name = "line", phases = [
    CreateState{Pure}(2, Qubit(), "Up"),
    Gates(gates = X(1), final_measurements = Line() => [Z(1), Z(2)]),
]))
```

[`TensorMixedStates.check_list`](@ref) refuses a list of measurements it cannot take before
anything is measured, and [`TensorMixedStates.resume_mark`](@ref), with `open(destination,
mark)`, lets a resumed run go on where the checkpoint left it, rather than start anew.

## States of one's own

A state of your own, of a representation of its own, is a subtype of
[`TensorMixedStates.AbstractState`](@ref) with a field `system` holding the `System` of its
physical sites, which a phase of your own creates, see [Phases of one's own](@ref own-phases).
TMS reads that system to save the state, to choose the threading, to log what weakening
changes, and to put the reference of `Fidelity` or `Overlap` on the system of the measured
state. It need not be the one the state keeps its tensors on: a mixed state held as a pure state
on a doubled system, a purification, has the system of its physical sites in `system`, and the
doubled one in a field of its own.

The phases and the measurements reach the state through the methods you give:

- `expect(state, op::TensorMixedStates.IndexedOp)` computes the expectation value of an
  operator placed on sites, normalised by the trace, which `measure`, and so `output`, asks for,
  those of the operators of one site and of the pairs of them included. The operator is given
  as the user wrote it, on the physical sites, which [`map_sites`](@ref) places on another
  system. It gives the value as computed: `measure`, `expect1` and `expect2` make that of a
  Hermitian operator real;
- `apply(op::TensorMixedStates.IndexedOp, state; kwargs...)` applies the gates of a `Gates` phase;
- `TensorMixedStates.write_state(group, state)` and
  `TensorMixedStates.read_state(::Type{MyState}, group, sites)` save the state and read it back,
  for `save_state`, `load_state`, `SaveState`, `LoadState` and the checkpoints: the first writes
  in the HDF5 group, where its type and its sites are already written, and the second rebuilds
  it from them. The type is written by its name alone: a type with parameters writes them in
  the group as well, and is read back by a method for the type without them,
  `read_state(::Type{<:MyState}, group, sites)`. A state holding a `State` writes it in a
  subgroup, `write_state(create_group(group, "inner"), inner)`, `create_group` coming from
  HDF5, and reads it back with `read_state(State{Pure}, group["inner"], sites)`, given the
  sites it lies on.

The measurements of the library go through `trace`, `norm`, `entanglement_entropy` and the
like: a state gets what its type has a method for, and a measurement or a phase it does not
support raises a `MethodError`.
