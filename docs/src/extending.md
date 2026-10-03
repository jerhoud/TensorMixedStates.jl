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
sites or an expression, or with [`Operator`](@ref) for one of several sites. A measurement
of your own is a [`StateFunc`](@ref), a function of the state, or a [`TimeFunc`](@ref), a
function of the simulation time.

## [Phases of one's own](@id own-phases)

A phase of your own is a struct with the three fields every phase has, `name`, `time_start`
and `final_measures`, and a method of [`TensorMixedStates.run_phase`](@ref) for it, which returns the
simulation the phase leaves behind. The full name is needed, `run_phase` not being exported.
Written as a loop with [`run_steps`](@ref), the phase is checkpointed, stopped and resumed between two
steps, as those of the library are between two sweeps:

```julia
using TensorMixedStates, .Qubits

Base.@kwdef struct Kicks
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
                                  evolver, coefs, nsweeps)
    dt = phase.duration / nsweeps
    return run_steps(sim, nsweeps) do sim, k
        sim = tdvp(evolver, dt, sim; phase.limits)
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
one, `nothing` otherwise.

## Representations of one's own

A representation of your own is a subtype of [`Representation`](@ref), which `CreateState`
takes as its `type`, and needs nothing more. Its states are of a type of your own, a subtype
of [`AbstractState`](@ref) with a field `system`, the `System` of their physical sites, which
the package reads to save a state, its sites going in the file, to choose the threading, to
log what a `Weaken` phase changes, and to put the reference of `Fidelity` or `Overlap` on the
system of the measured state. It need not be the system the representation keeps its tensors
on: a mixed state held as a pure state on a doubled system, a purification, has the system of
its physical sites in `system`, and the doubled one in a field of its own. The phases and the
measurements reach these states through the methods you give:

- `TensorMixedStates.run_phase(sim, phase::CreateState{MyRepresentation})` creates the state
  from the fields of the phase;
- `TensorMixedStates.expect_norm(state, terms)`, `expect1(state, ops)` and
  `expect2(state, pairs)` compute what `measure`, and so `output`, asks for: the expectation
  values of the terms of the operators, already simplified, those of one site operators on
  every site, and the correlations of pairs of them;
- `apply(gates, state; limits)` applies the gates of a `Gates` phase;
- `TensorMixedStates.evolve(algo, state, sim, phase; evolver, coefs, nsweeps)` evolves the state
  in an `Evolve` phase, for each algorithm it supports;
- `TensorMixedStates.write_state(group, state)` and
  `TensorMixedStates.read_state(::Type{MyState}, group, sites, system)` save the state and read
  it back, for `SaveState`, `LoadState` and the checkpoints: the first writes its tensors in the
  HDF5 group, where its type and its sites are already written, and the second rebuilds it from
  them, on `system` when that is not `nothing`. The type is written by its name alone: a type
  with parameters writes them in the group as well, and is read back by a method for the type
  without them, `read_state(::Type{<:MyState}, group, sites, system)`.

The other phases go through `truncate`, `mix`, `partial_trace`, `weaken`, `dmrg` and
`steady_state`, and the state functions of the measurements through `trace`, `norm`,
`entanglement_entropy` and the like. A state gets what its type has a method for, and a phase or
a measurement it does not support raises a `MethodError`.
