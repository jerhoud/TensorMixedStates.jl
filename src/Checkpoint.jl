# Checkpointing a running simulation: the fingerprint that identifies its phases, when a
# checkpoint is due or a stop requested, and how a checkpoint is written, committed and loaded
# back.

"""
    checkpoint_file_version

the layout version of a checkpoint; `load_checkpoint` refuses a checkpoint of another one.
"""
const checkpoint_file_version = 3

############### the fingerprint of the phases ###############

"""
    fnv_offset

the offset basis of the 64 bit FNV-1a hash, by which the phases are fingerprinted. Its
definition is fixed, unlike that of `Base.hash`, which changes between versions of Julia and
would make a checkpoint look like another simulation's after an upgrade.
"""
const fnv_offset = 0xcbf29ce484222325

"""
    fnv_prime

the prime of the FNV-1a hash, see `fnv_offset`.
"""
const fnv_prime = 0x00000100000001b3

"""
    fnv(h, bytes)

the FNV-1a hash `h` with the given bytes mixed in.
"""
fnv(h::UInt64, bytes) = foldl((h, b) -> (h ⊻ b) * fnv_prime, bytes; init = h)

"""
    fnv_mix(h, x)

the FNV-1a hash `h` with the value `x` mixed in, by bytes none of which depends on the version
of Julia: an integer by its 64 bits, or by its digits beyond them, a float by the bits of its
`Float64`, a rational or a complex number by its two parts, a string by its length and its code
units, anything else by the string it prints as.
"""
fnv_mix(h::UInt64, x::UInt64) = fnv(h, reinterpret(UInt8, [x]))
fnv_mix(h::UInt64, x::Integer) =
    typemin(Int64) ≤ x ≤ typemax(Int64) ? fnv(h, reinterpret(UInt8, [Int64(x)])) : fnv_mix(h, string(x))
# a signed zero mixes in as the zero it is equal to, `==` holding the two phases equal
fnv_mix(h::UInt64, x::AbstractFloat) = fnv(h, reinterpret(UInt8, [no_signed_zero(Float64(x))]))
fnv_mix(h::UInt64, x::Rational) = fnv_mix(fnv_mix(h, numerator(x)), denominator(x))
fnv_mix(h::UInt64, x::Complex) = fnv_mix(fnv_mix(h, real(x)), imag(x))
fnv_mix(h::UInt64, x::AbstractString) = fnv(fnv_mix(h, ncodeunits(x)), codeunits(x))
fnv_mix(h::UInt64, x::Union{Number, Symbol, Char}) = fnv_mix(h, string(x))
fnv_mix(h::UInt64, ::Nothing) = fnv_mix(h, "nothing")

"""
    type_key(T)

a type written by the full path of its module, its name and its parameters, which does not
depend on the names visible where it is written. `string(typeof(Qubit()))` gives `Qubit` or
`TensorMixedStates.Qubit` depending on the `using` of the program, which would give one
simulation two fingerprints.
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

"""
    type_parameter_key(p)

a parameter of a type as `type_key` writes it: a type or a type variable by its `type_key`, a
value, such as the dimension of an array or the names of a named tuple, as `repr` writes it.
"""
type_parameter_key(p::Union{Type, TypeVar}) = type_key(p)
type_parameter_key(p) = repr(p)

"""
    phase_hash(h, x)

the hash `h` with `x` mixed in, walking into its fields, see `phases_id`.
"""
phase_hash(h::UInt64, x::Union{Number, AbstractString, Symbol, Char, Nothing}) = fnv_mix(h, x)
phase_hash(h::UInt64, x::Type) = fnv_mix(h, type_key(x))
# an enumeration value has no field to walk into, its name is what it is
phase_hash(h::UInt64, x::Enum) = fnv_mix(h, type_key(typeof(x)) * "." * string(x))
phase_hash(h::UInt64, ::Index) = h
phase_hash(h::UInt64, x::System) = phase_hash(fnv_mix(h, "System"), x.sites)
# a State is its system and its tensors: `preobs` is a cache filled as measurements are
# made, so the same state would hash differently once it has been measured
phase_hash(h::UInt64, x::State{R}) where R =
    phase_hash(phase_hash(fnv_mix(h, "State{$R}"), x.system), x.state)
phase_hash(h::UInt64, x::Function) = fnv_mix(h, type_key(typeof(x)))
phase_hash(h::UInt64, x::Union{Tuple, Pair}) = foldl(phase_hash, (x...,); init = fnv_mix(h, "()"))
phase_hash(h::UInt64, x::AbstractArray) = foldl(phase_hash, x; init = foldl(fnv_mix, size(x); init = h))

# what a dictionary or a set holds, in no order: its storage order is not part of it, and its
# fields are the internals of a hash table, some of them undefined
function phase_hash(h::UInt64, x::Union{AbstractDict, AbstractSet})
    s = zero(UInt64)
    for y in x
        s += phase_hash(zero(UInt64), y)
    end
    return fnv_mix(fnv_mix(h, type_key(typeof(x))), s)
end

function phase_hash(h::UInt64, x)
    h = fnv_mix(h, type_key(typeof(x)))
    for f in fieldnames(typeof(x))
        if isdefined(x, f)
            h = phase_hash(h, getfield(x, f))
        end
    end
    return h
end

"""
    phases_id(phases)

a fingerprint of the phases of a simulation, by which a checkpoint tells whether it belongs
to the simulation being run. The phases are walked field by field, the structure being read
from the types, so a field added to a phase counts without anything else to change, and
nothing depends on how phases are printed.

Functions are a blind spot: only their type is hashed, so the coefficients of a time
dependent evolver or the bodies of two `StateFunc` cannot be told apart. The type of an
anonymous function reflects where it sits in the source.

ITensor indices are left out, since they carry an identity drawn afresh in every session: a
`System` is what its sites are.
"""
phases_id(phases) = string(phase_hash(fnv_offset, phases))

############### commits ###############

"""
    Commit

a point a simulation can be resumed from, taken at once by `commit!` at the end of a sweep or
at a phase boundary, once everything the uninterrupted run writes before that point is
written. A checkpoint is a commit written down and nothing else, so whatever happens between
a commit and its writing, an interrupt included, the checkpoint holds the state, the counts
and the outputs of one and the same moment.

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
    Checkpointer(dir = "", id = ""; interval = 0, max_time = Inf)

the checkpointing machinery of a simulation, created by `runTMS` and carried by the
`Simulation`, so that the solvers reach it through their observers. A simulation stops
cleanly when `max_time` seconds have passed and, with a directory, when the file `stop`
appears in it or on an interrupt; with a directory, a checkpoint is written first.

# Fields

- `dir`:        the simulation directory, where the checkpoint is written, empty for none
- `id`:         the fingerprint of the phases, see `phases_id`: a checkpoint with another one
                is refused
- `interval`:   seconds between two checkpoints, zero or less for no periodic checkpoint
- `deadline`:   the `time()` after which the simulation stops cleanly, `Inf` for none
- `next`:       the `time()` of the next checkpoint
- `last`:       the last commit, which a checkpoint writes
- `resume`:     the commit a resumed run starts from, until the phase it belongs to has run
- `sweeps`:     whether the phase being run has read its resume point, which lets its sweeps
                be committed, see `resume_sweeps!`
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

whether the simulation should stop now: a stop is under way, the deadline is past, or the
stop file is present.
"""
stop_requested(c::Checkpointer) =
    c.stopping || time() ≥ c.deadline || (!isempty(c.dir) && isfile(stop_file(c)))

"""
    checkpoint_due(::Checkpointer)

whether a periodic checkpoint is due. An interval of zero or less means none, the rule
`sweep_due` applies to the sweep counters: a negative one would otherwise keep the next
checkpoint in the past and write the whole state on every sweep.
"""
checkpoint_due(c::Checkpointer) = c.interval > 0 && time() ≥ c.next

"""
    commit!(::Checkpointer, ::Outputs, phase, sweep, phase_time, time, state; energy)

record a point the simulation can be resumed from, see `Commit`. How far the destinations
have got is read here, at the same moment as the rest.
"""
commit!(c::Checkpointer, o::Outputs, phase::Int, sweep::Int, phase_time::Number, t::Number,
        state::Union{Nothing, State}; energy = nothing) =
    c.last = Commit(phase, sweep, phase_time, t, state, energy, output_marks(o))

"""
    checkpoint_json(dir)

the metadata file of the checkpoint of a directory.
"""
checkpoint_json(dir::String) = joinpath(dir, "checkpoint.json")

"""
    state_file(generation)

the name of the state file of the given generation, 1 or 2, see `write_checkpoint`
"""
state_file(g::Int) = "checkpoint-$g.h5"

"""
    write_checkpoint(::Checkpointer, ::Outputs)

write the last commit down. The state goes to whichever of `checkpoint-1.h5` and
`checkpoint-2.h5` the checkpoint on the disk does not name, and the metadata, which names it,
is renamed into place last. That rename alone changes which checkpoint is on the disk, so a
crash at any point leaves the previous checkpoint or this one, whole, and never pairs a new
state with the previous counts.

Nothing is written without a directory, nor before the first phase has made a state, and a
commit already on the disk is not written again: a phase that commits none of its sweeps is
checkpointed at its start however long it runs, see `resume_sweeps!`.
"""
function write_checkpoint(c::Checkpointer, o::Outputs)
    k = c.last
    if isempty(c.dir) || isnothing(k) || !(k.state isa State)
        return nothing
    end
    if k !== c.written
        g = c.generation == 1 ? 2 : 1
        file = state_file(g)
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
                # through `checkpoint_value`, as the times of the destinations, so that a
                # complex time with no imaginary part stays complex and keeps its two columns
                "phase_time" => checkpoint_value(k.phase_time),
                "time" => checkpoint_value(k.time),
                "state" => file,
                "outputs" => persist_outputs(o, k.reached),
            ))
        end
        mv(tjs, js; force = true)
        rm(joinpath(c.dir, state_file(3 - g)); force = true)
        c.generation = g
        c.written = k
    end
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
state file. A checkpoint of another version, or whose state file is missing, is refused.
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
            generation = file == state_file(2) ? 2 : 1)
end

"""
    resume_sweeps!(::Checkpointer)

the sweeps the phase being run has already done and the energy dmrg had reached at the last
of them: `(0, nothing)` unless it is the phase a resumed run starts from. The resume point is
then consumed, so that it is read once.

Calling it is also what lets the sweeps of the phase be committed, see `sweep_commit!`. A
resume hands the phase the state of its last committed sweep, and only a phase that reads the
sweeps done and starts its solver after them continues correctly from there; any other would
run all its sweeps again on that state. So the sweeps of a phase that never calls this are
not committed, and a checkpoint written while it runs resumes it from its start.
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

the tail of a per sweep schedule, of `noise` or of the fields of a `Limits`, for a dmrg
resuming at sweep `done + 1`, which would otherwise start the schedule over. A schedule that
runs out is continued with its last value, as ITensor does; a plain value is returned as it
is.

This is for `dmrg` alone, which is handed the whole schedule and counts its sweeps from 1 on
every call. The evolution solvers keep the sweep numbers of the phase across a resume and
pick their value sweep by sweep with `sweep_limits`, so a tail would skip part of the
schedule twice.
"""
resume_schedule(x, ::Int) = x
resume_schedule(x::Vector, done::Int) =
    isempty(x) || done < length(x) ? x[done + 1:end] : x[end:end]
resume_schedule(l::Limits, done::Int) =
    Limits(resume_schedule(l.cutoff, done), resume_schedule(l.maxdim, done),
           resume_schedule(l.mindim, done))
