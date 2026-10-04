# The phase types a simulation is made of, as runTMS takes them: CreateState, LoadState,
# SaveState, ToMixed, the evolutions Tdvp, ApproxW, Evolve and Gates, GroundState, SteadyState,
# PartialTrace and Weaken.

export Phases, Algo, CreateState, LoadState, SaveState, ToMixed, Tdvp, ApproxW, Evolve, Gates
export GroundState, Dmrg, PartialTrace, SteadyState, Weaken

"""
    CreateState(; type, system, state, randomize, seed, name, time_start, final_measures)
    CreateState{Pure|Mixed}(n, site, state; options...)
    CreateState{Pure|Mixed}(sites, state; options...)

a phase that creates the state of the simulation. The first phase of a simulation is this one
or `LoadState`.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `type`: the representation of the state, `Pure()` or `Mixed()`, or a `Representation` an
  extension defines, together with the method of `run_phase` creating its state
- `system`: the `System` of the state, unused when `state` is a `State`
- `state`: a description of the state, or a `State`, which is mixed if `type` asks for it (a
  mixed one cannot be made pure)
- `randomize`: the link dimension of a random state to create (default 0, none). With a
  `state`, a pure state is randomized from it and a mixed one drawn from a purification
  starting from it, see `RandomState`; a `State` can only be randomized into a pure state
- `seed`: the seed the global random generator is given at the start of the phase (default
  `nothing`, none). A checkpoint does not save the generator, so a run resumed after this
  phase does not set it again

# Examples

    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up")
    CreateState(type = Mixed(), system = System(3, Qubit()), state = ["Up", "Dn", "Up"])
    CreateState(type = Pure(), system = System(10, Qubit()), randomize = 50)
    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up", randomize = 50)
    CreateState{Mixed}(4, Fermion(conserve = N), ["Occ", "Emp", "Occ", "Emp"]; randomize = 16)
    CreateState{Pure}(10, Qubit(), "Up")                                      # simple form
    CreateState{Mixed}([Qubit(), Boson(4), Fermion()], ["Up", "2", "Occ"])    # other simple form
"""
@kwdef struct CreateState{R <: Representation}
    name::String = "Creating state"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    type::R
    system::Union{Nothing, System} = nothing
    state = nothing
    randomize::Int = 0
    seed::Union{Nothing, Int} = nothing
end

CreateState{R}(n, site, state; kwargs...) where R =
    CreateState(;type = R(), system = System(n, site), state, kwargs...)
CreateState{R}(sites, state; kwargs...) where R =
    CreateState(;type = R(), system = System(sites), state, kwargs...)

"""
    SaveState(; file, statename = "state", name, time_start, final_measures)

a phase that saves the state in an HDF5 file, see `save_state`. Several states can be saved in
one file under different `statename`; saving under a name already in the file replaces it.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `file`: the name of the HDF5 file to write to, taken in the simulation directory when it is
  relative, since `runTMS` runs the phases there
- `statename`: the name under which the state is stored in the file

# Examples

    SaveState(file = "myfile.h5")
    SaveState(file = "myfile.h5", statename = "after_evolution")
"""
@kwdef struct SaveState
    name::String = "Saving state"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    file::String
    statename::String = "state"
end

"""
    LoadState(; file, statename = "state", limits, name, time_start, final_measures)

a phase that loads the state from an HDF5 file written by `SaveState` or `save_state`, see
`load_state`, and truncates it to `limits`, unless they are the default `Limits()`.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `file`: the name of the HDF5 file to read from, taken in the simulation directory when it is
  relative: a state another simulation saved is found under `../othername/`
- `statename`: the name under which the state is stored in the file
- `limits`: the truncation applied to the state once loaded, see `Limits` (default `Limits()`,
  none)

# Examples

    LoadState(file = "myfile.h5")
    LoadState(file = "myfile.h5", statename = "after_evolution")
"""
@kwdef struct LoadState
    name::String = "Loading state"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    file::String
    statename::String = "state"
    limits::Limits = Limits()
end

"""
    ToMixed(; limits, name, time_start, final_measures)

a phase that switches the state to the mixed representation and truncates it to `limits`; a
state already mixed is only truncated.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `limits`: constraints on the mixed state, see `Limits` (default `Limits()`, none)

# Examples

    ToMixed()
    ToMixed(limits = Limits(cutoff = 1e-10, maxdim = 10))
"""
@kwdef struct ToMixed
    name::String = "Switching to mixed state representation"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    limits::Limits = Limits()
end

"""
    abstract type Algo

the supertype of the time evolution algorithms the `algo` field of `Evolve` takes: `Tdvp` and
`ApproxW`, and those an extension defines with a method of `TensorMixedStates.evolve`.
"""
abstract type Algo end

"""
    Tdvp(; n_expand = 0, n_hermitianize = 0, krylov = Krylov())

the tdvp algorithm, for the `algo` field of `Evolve`, see `tdvp`.

- `n_expand`: enlarge the bond dimension of the state by a global Krylov expansion before the
  first step and then every `n_expand` steps (default 0, never)
- `n_hermitianize`: make a mixed state hermitian every `n_hermitianize` steps (default 0,
  never)
- `krylov`: the parameters of the Krylov exponentiation of each local step, see `Krylov`
  (default `Krylov()`)

# Examples

    Tdvp()
    Tdvp(n_expand = 5)                   # tdvp with expansion steps every 5 steps
    Tdvp(n_hermitianize = 3)             # tdvp, make hermitian every 3 steps
    Tdvp(krylov = Krylov(tol = 1e-10))   # tdvp, local steps at a lower precision
"""
@kwdef struct Tdvp <: Algo
    n_expand::Int = 0
    n_hermitianize::Int = 0
    krylov::Krylov = Krylov()
end

"""
    ApproxW(; order, w = 2, n_hermitianize = 0, apply_algo = "densitymatrix")

time evolution by WI or WII approximations of the exponential, combined into an approximation
of the given order, for the `algo` field of `Evolve`, see `approx_W`.

- `order`: the order of the approximation, from 1 to 4, required
- `w`: 1 or 2 for WI or WII (default 2)
- `n_hermitianize`: make a mixed state hermitian every `n_hermitianize` steps (default 0,
  never)
- `apply_algo`: the algorithm of the product of the state by each MPO, `"densitymatrix"`
  (default) or `"naive"`, see `approx_W`

# Examples

    ApproxW(order = 2)                       # order 2, WII
    ApproxW(order = 4, w = 1)                # order 4, WI
    ApproxW(order = 4, n_hermitianize = 3)   # order 4, make hermitian every 3 steps
    ApproxW(order = 2, apply_algo = "naive") # order 2, naive products
"""
@kwdef struct ApproxW <: Algo
    order::Int
    w::Int = 2
    n_hermitianize::Int = 0
    apply_algo::String = "densitymatrix"
    # checked when the phase is written rather than when it runs: corrected then, it could not
    # resume its checkpoint, which belongs to a simulation of other phases
    function ApproxW(order, w, n_hermitianize, apply_algo)
        check_w_approx(order, w)
        check_apply_algo(apply_algo)
        return new(order, w, n_hermitianize, apply_algo)
    end
end

# an algorithm is printed on one line, field by field, read from its type as a phase is, so
# that a field added to it shows up in the log and in `prog.jl` with nothing else to change
show(io::IO, s::Algo) =
    print(io, nameof(typeof(s)), "(",
          join(("$f = $(repr(getfield(s, f)))" for f in fieldnames(typeof(s))), ", "), ")")

"""
    Evolve(; duration, time_step, algo, evolver, measures, measures_period, limits, options...)

a phase of time evolution.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `limits`: constraints on the state, see `Limits` (default `Limits()`, none)
- `duration`: the duration of the evolution
- `time_step`: the time step, adjusted to the nearest one that divides the duration into a
  whole number of steps, and taken with the sign of the duration (the phase is skipped when
  that number is zero)
- `algo`: the algorithm, `Tdvp(...)` or `ApproxW(...)`
- `evolver`: `-im * H` for a hamiltonian `H`, plus dissipators for a mixed state, or
  `evolvers => coefs` for a time dependent one, see the `coefs` option of `tdvp`
- `measures`: the measurements to make during the evolution, see `output` (default `[]`),
  after every `measures_period` time steps. The state the phase starts from is not measured
  here: the `final_measures` of the phase before measure it
- `measures_period`: the number of time steps between two measurements (default 1)

# Examples

    Evolve(duration = 2., time_step = 0.1, algo = Tdvp(),
           evolver = -im * (Z(1)Z(2) + Z(2)Z(3)), measures = "data" => [X, Y, Z])
"""
@kwdef struct Evolve
    name::String = "Time evolution"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    limits::Limits = Limits()
    duration::Number
    time_step::Number
    algo::Algo
    evolver::Union{IndexedOp, Pair}
    measures_period::Int = 1
    measures = []
    # refused here rather than when it runs, as the fields of ApproxW
    function Evolve(name, time_start, final_measures, limits, duration, time_step, algo, evolver,
                    measures_period, measures)
        if iszero(time_step)
            error("the time step of an Evolve cannot be zero")
        end
        # a single term and its function of time, written without the vectors, which failed on
        # a MethodError at the first step, and a function per term, which only the first step
        # checked
        if evolver isa Pair
            ops, fs = evolver
            if !(ops isa AbstractVector) && !(fs isa AbstractVector)
                evolver = [ops] => [fs]
            elseif !(ops isa AbstractVector && fs isa AbstractVector) || length(ops) ≠ length(fs)
                error("a time dependent evolver gives one function of time per term, as " *
                      "[A, B] => [f, g] or A => f")
            end
        end
        return new(name, time_start, final_measures, limits, duration, time_step, algo, evolver,
                   measures_period, measures)
    end
end

"""
    Gates(; gates, limits, name, time_start, final_measures)

a phase that applies gates to the state.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `gates`: the gates to apply
- `limits`: the truncations made while applying a gate of several sites, see `apply` (default
  `Limits()`, none)

# Examples

    Gates(gates = controlled(X)(1, 3)*controlled(Z)(2, 4), limits = Limits(cutoff=1e-10, maxdim = 20))
"""
@kwdef struct Gates
    name::String = "Applying gates"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    gates::IndexedOp
    limits::Limits = Limits()
end

"""
    GroundState(; hamiltonian, limits, nsweeps, noise, tolerance, measures, options...)

a phase that searches the ground state of a hamiltonian by dmrg, see `dmrg`, on a pure state.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `hamiltonian`: the hamiltonian whose ground state is searched
- `limits`: constraints on the state, see `Limits`, required
- `nsweeps`: the maximum number of sweeps, required
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `krylov`: the parameters of the Krylov search of each local step, see `Krylov` (default
  `Krylov()`)
- `measures`: the measurements to make during the search, see `output` (default `[]`)
- `measures_period`: the number of sweeps between two measurements (default 1)
- `tolerance`: the search stops when the energy changes by less than this from one sweep to
  the next (default 0, never)

# Examples

    GroundState(hamiltonian = X(1)X(2), nsweeps = 10,
                limits = Limits(cutoff = 1e-10, maxdim = [10, 20, 30]), tolerance = 1e-6)
"""
@kwdef struct GroundState
    name::String = "Ground state computation using Dmrg"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    hamiltonian::IndexedOp{Pure}
    limits::Limits
    nsweeps::Int
    noise::Union{Float64, Vector{Float64}} = 0.
    krylov::Krylov = Krylov()
    measures = []
    measures_period::Int = 1
    tolerance::Real = 0.
    # a noise is a real number, an integer one included, stored as the Float64 dmrg takes
    GroundState(name, time_start, final_measures, hamiltonian, limits, nsweeps,
                noise::Union{Real, AbstractVector{<:Real}}, krylov, measures, measures_period,
                tolerance) =
        new(name, time_start, final_measures, hamiltonian, limits, nsweeps,
            noise isa AbstractVector ? Vector{Float64}(noise) : Float64(noise),
            krylov, measures, measures_period, tolerance)
end

# the docstring goes through `@doc` rather than sitting above the call, because the macro
# expands to a toplevel block and a docstring cannot be attached to one. It comes before the
# deprecation: after it, attaching the docstring used the deprecated name and warned
@doc """
    Dmrg

deprecated, use [`GroundState`](@ref) instead.

`Dmrg` is an alias of `GroundState`, which works until it is removed in a future version. Its
deprecation warning only shows with `--depwarn=yes`, as when running the tests, so this line
is the notice.
""" Dmrg

Base.@deprecate_binding Dmrg GroundState false ", use GroundState instead."

"""
    PartialTrace(; trace_positions | keep_positions, name, time_start, final_measures)

a phase that traces out part of the sites, given by exactly one of `trace_positions` and
`keep_positions`.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `trace_positions`: the sites to trace out
- `keep_positions`: the sites to keep, all the others being traced out

The state must be mixed, see `ToMixed`, and conserve nothing strongly, see `Weaken`. The sites
kept make a new system, numbered from 1 in their order, which the operators of the phases
after this one refer to.

# Examples

    PartialTrace(trace_positions = [2, 3, 6])
    PartialTrace(keep_positions = [1, 4, 5])
"""
@kwdef struct PartialTrace
    name::String = "Computing partial trace"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    trace_positions::Union{Nothing, Vector{Int}} = nothing
    keep_positions::Union{Nothing, Vector{Int}} = nothing
end

"""
    Weaken(; target, name, time_start, final_measures)

a phase that takes the state down to a lower level of conservation, see `weaken`. A phase may
evolve under a strong symmetry, which every dissipator commuting with the charge allows, and
the next one continue under a weak one, where a jump that moves the charge becomes possible.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `target`: what the state must still conserve, as `weaken` takes it (default `nothing`, one
  level down: every strong quantity made weak or, when none is strong, every quantity
  dropped)

# Examples

    Weaken()
    Weaken(target = (strong(Ntot), 2Sz))
    Weaken(target = ())                     # no charges at all
"""
@kwdef struct Weaken
    name::String = "Weakening the symmetries"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    target = nothing
end

"""
    SteadyState(; lindbladian, limits, nsweeps, noise, tolerance, measures, options...)

a phase that searches the steady state of a Lindbladian, see `steady_state`, on a mixed state.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `Phases`
- `lindbladian`: the Lindbladian ``L`` whose steady state is searched, of the form
  `-im * hamiltonian + dissipators`
- `mpo_limits`: the truncation of the MPO of ``L^\\dagger L`` (default `Limits()`, none)
- `mpo_algo`: the algorithm computing ``L^\\dagger L``, `"naive"` (default) or `"zipup"`
- `limits`: constraints on the state, see `Limits`, required
- `nsweeps`: the maximum number of sweeps, required
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `krylov`: the parameters of the Krylov search of each local step, see `Krylov` (default
  `Krylov(dim = 8, maxiter = 3)`, see `steady_state`)
- `measures`: the measurements to make during the search, see `output` (default `[]`)
- `measures_period`: the number of sweeps between two measurements (default 1)
- `tolerance`: the search stops when the dmrg energy changes by less than this from one sweep
  to the next (default 0, never)

# Examples

    SteadyState(
        lindbladian = -im * hamiltonian + dissipators,
        limits = Limits(cutoff = 1e-20, maxdim = [10, 20, 50]),
        nsweeps = 10,
        tolerance = 1e-5,
    )
"""
@kwdef struct SteadyState
    name::String = "Steady state optimization"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    lindbladian::IndexedOp{Mixed}
    mpo_limits::Limits = Limits()
    mpo_algo::String = "naive"
    limits::Limits
    nsweeps::Int
    noise::Union{Float64, Vector{Float64}} = 0.
    krylov::Krylov = Krylov(dim = 8, maxiter = 3)
    measures = []
    measures_period::Int = 1
    tolerance::Real = 0.
    # a noise is stored as the Float64 dmrg takes, as in GroundState, and mpo_algo is checked
    # when the phase is written, as the fields of ApproxW
    function SteadyState(name, time_start, final_measures, lindbladian, mpo_limits, mpo_algo,
                         limits, nsweeps, noise::Union{Real, AbstractVector{<:Real}}, krylov,
                         measures, measures_period, tolerance)
        check_mpo_algo(mpo_algo)
        return new(name, time_start, final_measures, lindbladian, mpo_limits, mpo_algo, limits,
                   nsweeps, noise isa AbstractVector ? Vector{Float64}(noise) : Float64(noise),
                   krylov, measures, measures_period, tolerance)
    end
end

"""
    Phases

the union of the phase types of the library: `CreateState`, `SaveState`, `LoadState`,
`ToMixed`, `Evolve`, `Gates`, `GroundState`, `PartialTrace`, `SteadyState` and `Weaken`.
Every phase has at least these three fields:

- `name`: the name of the phase, written in the log
- `time_start`: the simulation time the clock is set to when the phase starts (default
  `nothing`, keeping the current time)
- `final_measures`: the measurements to make at the end of the phase, see `output` (default
  `[]`)

A phase of your own can be defined, see `TensorMixedStates.run_phase`.
"""
const Phases = Union{CreateState, SaveState, LoadState, ToMixed, Evolve, Gates, GroundState, PartialTrace, SteadyState, Weaken}

# The phases are printed field by field, read from the type rather than written out one by
# one, so that a field added to a phase shows up in the log and in `prog.jl` without
# anything else to change.
function show(io::IO, s::Phases)
    t = typeof(s)
    print(io, "\n", nameof(t))
    if !isempty(t.parameters)
        print(io, "{", join(nameof.(t.parameters), ", "), "}")
    end
    print(io, "(")
    fs = fieldnames(t)
    for (i, f) in enumerate(fs)
        print(io, "\n    ", f, " = ", repr(getfield(s, f)), i < length(fs) ? "," : ")")
    end
end