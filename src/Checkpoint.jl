# 2: the complex numbers and matrices of `Data` destinations are marked, see
# `checkpoint_value`, and a measurement becoming complex would continue a file of version 1
# in another layout
# 3: each value of a dictionary destination records the call of `output` it came from, see
# `next_event`; the state is in a file the metadata names, see `write_checkpoint`
const checkpoint_file_version = 3

############### the fingerprint of the phases ###############

# FNV-1a, whose definition is fixed: `Base.hash` is not, and changes between versions of
# Julia, which refused a checkpoint after an upgrade as belonging to another simulation
const fnv_offset = 0xcbf29ce484222325
const fnv_prime = 0x00000100000001b3

fnv(h::UInt64, bytes) = foldl((h, b) -> (h ⊻ b) * fnv_prime, bytes; init = h)

# the bytes a value is mixed in by, none of which depends on the version of Julia
mix(h::UInt64, x::UInt64) = fnv(h, reinterpret(UInt8, [x]))
mix(h::UInt64, x::Integer) =
    typemin(Int64) ≤ x ≤ typemax(Int64) ? fnv(h, reinterpret(UInt8, [Int64(x)])) : mix(h, string(x))
# a signed zero mixes in as the zero it is equal to, `==` holding the two phases equal
mix(h::UInt64, x::AbstractFloat) = fnv(h, reinterpret(UInt8, [no_signed_zero(Float64(x))]))
mix(h::UInt64, x::Rational) = mix(mix(h, numerator(x)), denominator(x))
mix(h::UInt64, x::Complex) = mix(mix(h, real(x)), imag(x))
mix(h::UInt64, x::AbstractString) = fnv(mix(h, ncodeunits(x)), codeunits(x))
mix(h::UInt64, x::Union{Number, Symbol, Char}) = mix(h, string(x))
mix(h::UInt64, ::Nothing) = mix(h, "nothing")

"""
    type_key(T)

a type written by the full path of its module, its name and its parameters, which does not
depend on the names visible where it is written: `string(typeof(Qubit()))` gives `Qubit` or
`TensorMixedStates.Qubit` depending on the `using` of the program, so the same simulation took
two fingerprints and its checkpoint was refused as belonging to another one.
"""
function type_key(T::DataType)
    name = join((fullname(parentmodule(T))..., nameof(T)), ".")
    if isempty(T.parameters)
        return name
    end
    return name * "{" * join(map(type_parameter_key, T.parameters), ",") * "}"
end
type_key(T::UnionAll) = type_key(Base.unwrap_unionall(T))
type_key(T::Union) = "Union{" * join(sort(map(type_key, Base.uniontypes(T))), ",") * "}"
type_key(T::TypeVar) = string(T.name)
# `Union{}`, the one type of none of the kinds above
type_key(T::Type) = string(T)

# a parameter is a type, a type variable, or a value such as the dimension of an array or the
# names of a named tuple, which is written as it reads
type_parameter_key(p::Union{Type, TypeVar}) = type_key(p)
type_parameter_key(p) = repr(p)

"""
    phases_id(phases)

a fingerprint of the phases of a simulation, so that a checkpoint can tell whether it
belongs to the simulation being run. It walks the phases field by field, reading the
structure from the types themselves, so a field added to a phase counts without anything
else to change and nothing depends on how phases are printed.

Functions are a blind spot: coefficients of a time dependent evolver, or the body of a
`StateFunc`, cannot be told apart. Only their type is hashed, which for an anonymous
function reflects where it sits in the source.

ITensor indices are left out, since they carry an identity drawn afresh in every session
and say nothing about the simulation: a `System` is what its sites are.
"""
phases_id(phases) = string(phase_hash(fnv_offset, phases))

phase_hash(h::UInt64, x::Union{Number, AbstractString, Symbol, Char, Nothing}) = mix(h, x)
phase_hash(h::UInt64, x::Type) = mix(h, type_key(x))
# an enumeration value has no field to walk into, its name is what it is
phase_hash(h::UInt64, x::Enum) = mix(h, type_key(typeof(x)) * "." * string(x))
phase_hash(h::UInt64, ::Index) = h
phase_hash(h::UInt64, x::System) = phase_hash(mix(h, "System"), x.sites)
# a State is its system and its tensors: `preobs` is a cache filled as measurements are
# made, so the same state would hash differently once it has been measured
phase_hash(h::UInt64, x::State{R}) where R =
    phase_hash(phase_hash(mix(h, "State{$R}"), x.system), x.state)
phase_hash(h::UInt64, x::Function) = mix(h, type_key(typeof(x)))
phase_hash(h::UInt64, x::Union{Tuple, Pair}) = foldl(phase_hash, (x...,); init = mix(h, "()"))
phase_hash(h::UInt64, x::AbstractArray) = foldl(phase_hash, x; init = foldl(mix, size(x); init = h))

# what a dictionary or a set holds, in no order: its storage order is not part of it, and its
# fields are the internals of a hash table, some of them undefined
function phase_hash(h::UInt64, x::Union{AbstractDict, AbstractSet})
    s = zero(UInt64)
    for y in x
        s += phase_hash(zero(UInt64), y)
    end
    return mix(mix(h, type_key(typeof(x))), s)
end

function phase_hash(h::UInt64, x)
    h = mix(h, type_key(typeof(x)))
    for f in fieldnames(typeof(x))
        if isdefined(x, f)
            h = phase_hash(h, getfield(x, f))
        end
    end
    return h
end

############### commits ###############

"""
    Commit

a point a simulation can be resumed from: the phase it is in, the sweeps of it done, the time
the phase started from, the time reached, the state, what dmrg compares its next sweep with,
and how far every destination had got. It is taken at once, in `commit!`, at the end of a
sweep or at a phase boundary, once everything the uninterrupted run writes before that point
is written.

A checkpoint is a commit written down, and nothing else: whatever happens between a commit
and the writing of it, an interrupt included, the checkpoint holds the state, the counts and
the outputs of one and the same moment.

- `phase`:      index of the phase to resume, one past the last when the run is over
- `sweep`:      sweeps of that phase done
- `phase_time`: simulation time the phase started from, which its solver counts sweeps from
- `time`:       simulation time reached
- `state`:      the state reached, `nothing` before the first phase has made one
- `energy`:     the energy of the last dmrg sweep, which a resumed search compares its first
                sweep with
- `reached`:    how far every destination had got, see `output_marks`
"""
struct Commit
    phase::Int
    sweep::Int
    phase_time::Number
    time::Number
    state::Union{Nothing, State}
    energy::Union{Nothing, Float64}
    reached::NamedTuple
end

"""
    struct Checkpointer

holds the checkpointing machinery of a simulation. One is created by `runTMS` and carried
by the `Simulation`, so that the solvers can reach it through their observers.

A simulation stops cleanly when its deadline is past, when the file `<simulation>/stop`
appears, or on an interrupt. With a directory, a checkpoint is written first.

# Fields

- `dir`:        the simulation directory, where the checkpoint is written, empty for none
- `id`:         a fingerprint of the phases, a checkpoint of another simulation is refused
- `interval`:   seconds between two checkpoints, zero or less disables checkpointing
- `deadline`:   time after which the simulation stops cleanly, `Inf` for no limit
- `next`:       time of the next checkpoint
- `last`:       the last commit, which a checkpoint writes
- `resume`:     the commit a resumed run starts from, until the phase it belongs to has run
- `sweeps`:     whether the phase being run has read its resume point, which is what lets its
                sweeps be committed, see `resume_sweeps!`
- `written`:    the commit the checkpoint on the disk holds, which is not written again
- `generation`: which of the two state files the checkpoint on the disk names
- `stopping`:   set once a stop has been requested, so that every loop unwinds
"""
mutable struct Checkpointer
    dir::String
    id::String
    interval::Float64
    deadline::Float64
    next::Float64
    last::Union{Nothing, Commit}
    resume::Union{Nothing, Commit}
    sweeps::Bool
    written::Union{Nothing, Commit}
    generation::Int
    stopping::Bool
end

Checkpointer(dir::String = "", id::String = ""; interval::Real = 0, max_time::Real = Inf) =
    Checkpointer(dir, id, interval,
                 max_time == Inf ? Inf : time() + max_time,
                 interval ≤ 0 ? Inf : time() + interval,
                 nothing, nothing, false, nothing, 0, false)

"""
    stop_file(::Checkpointer)

the file whose presence asks the simulation to stop cleanly
"""
stop_file(c::Checkpointer) = joinpath(c.dir, "stop")

"""
    stop_requested(::Checkpointer)

whether the simulation should stop now
"""
stop_requested(c::Checkpointer) =
    c.stopping || time() ≥ c.deadline || (!isempty(c.dir) && isfile(stop_file(c)))

"""
    checkpoint_due(::Checkpointer)

whether enough time has passed since the last checkpoint

An interval of zero or less means no checkpointing, the same rule `sweep_due` applies to
the sweep counters. A negative interval used to put the next checkpoint in the past and
keep it there, writing the whole state to disk on every sweep.
"""
checkpoint_due(c::Checkpointer) = c.interval > 0 && time() ≥ c.next

"""
    commit!(::Checkpointer, ::Outputs, phase, sweep, phase_time, time, state; energy)

record a point the simulation can be resumed from, see `Commit`. The destinations are read
here, at the same moment as the rest.
"""
commit!(c::Checkpointer, o::Outputs, phase::Int, sweep::Int, phase_time::Number, t::Number,
        state::Union{Nothing, State}; energy = nothing) =
    c.last = Commit(phase, sweep, phase_time, t, state, energy, output_marks(o))

checkpoint_json(dir::String) = joinpath(dir, "checkpoint.json")

"""
    write_checkpoint(::Checkpointer, ::Outputs)

write the last commit down. The state goes to one of two files, `checkpoint-1.h5` and
`checkpoint-2.h5`, the one the checkpoint on the disk does not name, and the metadata, which
names it, is renamed into place last. That rename is the only step that changes which
checkpoint is on the disk, so a crash at any point leaves either the previous checkpoint or
this one, whole: renaming the state and then the metadata used to pair, after a kill between
the two, the new state with the previous counts, and the resume ran again the sweeps the state
already held.

Nothing is written without a directory, nor before the first phase has made a state, and a
commit already on the disk is not written again: a phase that commits none of its sweeps is
checkpointed at its start however long it runs, see `resume_sweeps!`.
"""
function write_checkpoint(c::Checkpointer, o::Outputs)
    k = c.last
    if isempty(c.dir) || isnothing(k) || !(k.state isa State)
        return nothing
    end
    if k === c.written
        c.next = time() + c.interval
        return nothing
    end
    g = c.generation == 1 ? 2 : 1
    file = "checkpoint-$g.h5"
    h5 = joinpath(c.dir, file)
    rm(h5; force = true)
    save_state(h5, "checkpoint", k.state)
    js = checkpoint_json(c.dir)
    tjs = js * ".tmp"
    open(tjs, "w") do io
        JSON.print(io, Dict(
            "version" => checkpoint_file_version,
            "id" => c.id,
            "phase" => k.phase,
            "sweep" => k.sweep,
            "energy" => checkpoint_value(k.energy),
            # through `checkpoint_value`, as the times of the destinations, so that a complex
            # time with no imaginary part stays complex and keeps its two columns
            "phase_time" => checkpoint_value(k.phase_time),
            "time" => checkpoint_value(k.time),
            "state" => file,
            "outputs" => persist_outputs(o, k.reached),
        ))
    end
    mv(tjs, js; force = true)
    rm(joinpath(c.dir, "checkpoint-$(3 - g).h5"); force = true)
    c.generation = g
    c.written = k
    c.next = time() + c.interval
    return nothing
end

"""
    has_checkpoint(dir)

whether a checkpoint is present in the given directory
"""
has_checkpoint(dir::String) = isfile(checkpoint_json(dir))

"""
    load_checkpoint(dir)

read the checkpoint of the given directory, as a named tuple of the fields of the commit it
records (`phase`, `sweep`, `phase_time`, `time`, `state`, `energy`), the fingerprint `id` of
its phases, the `outputs` to put back with `restore_outputs!`, and the `generation` of its
state file.
"""
function load_checkpoint(dir::String)
    meta = JSON.parsefile(checkpoint_json(dir))
    if meta["version"] ≠ checkpoint_file_version
        error("checkpoint of $dir has version $(meta["version"]), expected $checkpoint_file_version")
    end
    file = meta["state"]
    if !isfile(joinpath(dir, file))
        error("the checkpoint of $dir names the state file $file, which is missing")
    end
    return (phase = Int(meta["phase"]), sweep = Int(meta["sweep"]),
            phase_time = restored_value(meta["phase_time"]), time = restored_value(meta["time"]),
            state = load_state(joinpath(dir, file), "checkpoint"),
            energy = restored_value(meta["energy"]), id = meta["id"], outputs = meta["outputs"],
            generation = file == "checkpoint-2.h5" ? 2 : 1)
end

"""
    resume_sweeps!(::Checkpointer)

the sweeps the phase being run has done already and the energy dmrg had reached at the last
of them, `(0, nothing)` unless it is the phase a resumed run starts from. The resume point is
consumed, so that it is read once.

Reading it is also what lets the sweeps of the phase be committed, see `sweep_commit!`. A
resume hands the phase the state reached at the sweep committed, and only a phase that reads
the sweeps done and starts its solver after them continues correctly from there: any other
would run all its sweeps again on that state. So the sweeps of a phase that never calls this
are not committed, and a checkpoint written while it runs resumes it from its start.
"""
function resume_sweeps!(c::Checkpointer)
    c.sweeps = true
    r = c.resume
    if isnothing(r) || isnothing(c.last) || r.phase ≠ c.last.phase
        return (0, nothing)
    end
    c.resume = nothing
    return (r.sweep, r.energy)
end

"""
    resume_schedule(x, done)

`maxdim`, `mindim`, `cutoff` and `noise` may be given one value per sweep. A phase
resuming at sweep `done + 1` has to be handed the tail of those schedules, otherwise it
would start them over and run the remaining sweeps with the wrong ones. A schedule that
runs out is continued with its last value, which is what ITensor does with a schedule
shorter than its sweeps.

This is for `dmrg` alone, which is handed the whole schedule and counts its sweeps from 1
on every call. The evolution solvers keep the sweep numbers of the phase across a resume
and pick their value out sweep by sweep with `sweep_limits`, so handing them a tail would
skip part of the schedule twice.
"""
resume_schedule(x, ::Int) = x
resume_schedule(x::Vector, done::Int) =
    isempty(x) || done < length(x) ? x[done + 1:end] : x[end:end]
resume_schedule(l::Limits, done::Int) =
    Limits(resume_schedule(l.cutoff, done), resume_schedule(l.maxdim, done),
           resume_schedule(l.mindim, done))
