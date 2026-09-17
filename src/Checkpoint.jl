const checkpoint_file_version = 1

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
phases_id(phases) = string(phase_hash(zero(UInt), phases))

phase_hash(h::UInt, x::Union{Number, AbstractString, Symbol, Char, Nothing}) = hash(x, h)
phase_hash(h::UInt, x::Type) = hash(string(x), h)
# an enumeration value has no field to walk into, its name is what it is
phase_hash(h::UInt, x::Enum) = hash(string(typeof(x), ".", x), h)
phase_hash(h::UInt, ::Index) = h
phase_hash(h::UInt, x::System) = phase_hash(hash("System", h), x.sites)
# a State is its system and its tensors: `preobs` is a cache filled as measurements are
# made, so the same state would hash differently once it has been measured
phase_hash(h::UInt, x::State{R}) where R =
    phase_hash(phase_hash(hash("State{$R}", h), x.system), x.state)
phase_hash(h::UInt, x::Function) = hash(string(typeof(x)), h)
phase_hash(h::UInt, x::Union{Tuple, Pair}) = foldl(phase_hash, (x...,); init = hash("()", h))
phase_hash(h::UInt, x::AbstractArray) = foldl(phase_hash, x; init = hash(size(x), h))

function phase_hash(h::UInt, x)
    h = hash(string(typeof(x)), h)
    for f in fieldnames(typeof(x))
        h = phase_hash(h, getfield(x, f))
    end
    return h
end

"""
    struct Checkpointer

holds the checkpointing machinery of a simulation. One is created by `runTMS` and carried
by the `Simulation`, so that the solvers can reach it through their observers.

A simulation stops cleanly when its deadline is past, when the file `<simulation>/stop`
appears, or on an interrupt. In all three cases a checkpoint is written first.

# Fields

- `dir`:       the simulation directory, where the checkpoint is written
- `id`:        a fingerprint of the phases, a checkpoint of another simulation is refused
- `interval`:  seconds between two checkpoints, `0` disables checkpointing
- `deadline`:  time after which the simulation stops cleanly, `Inf` for no limit
- `next`:      time of the next checkpoint
- `phase`:     index of the phase being run
- `sweep`:     sweeps done so far in that phase
- `phase_time`: simulation time at the start of that phase, what a resume restores
- `simtime`:   simulation time reached, used when a phase is cut short
- `state`:     the state as of the last sweep, so that an interrupt can still save it
- `resuming`:  set while the first phase of a resumed run has not started yet
- `skip`:      sweeps to skip when resuming the current phase, `0` when not resuming
- `appending`: whether output files are being continued rather than created
- `stopping`:  set once a stop has been requested, so that every loop unwinds
"""
mutable struct Checkpointer
    dir::String
    id::String
    interval::Float64
    deadline::Float64
    next::Float64
    phase::Int
    sweep::Int
    phase_time::Number
    simtime::Number
    state::Union{Nothing, State}
    resuming::Bool
    skip::Int
    appending::Bool
    stopping::Bool
end

Checkpointer(dir::String = "", id::String = ""; interval::Real = 0, max_time::Real = Inf) =
    Checkpointer(dir, id, interval,
                 max_time == Inf ? Inf : time() + max_time,
                 interval == 0 ? Inf : time() + interval,
                 1, 0, 0., 0., nothing, false, 0, false, false)

"""
    stop_file(::Checkpointer)

the file whose presence asks the simulation to stop cleanly
"""
stop_file(c::Checkpointer) = joinpath(c.dir, "stop")

checkpoint_h5(dir::String) = joinpath(dir, "checkpoint.h5")
checkpoint_json(dir::String) = joinpath(dir, "checkpoint.json")

"""
    has_checkpoint(dir)

whether a complete checkpoint is present in the given directory
"""
has_checkpoint(dir::String) = isfile(checkpoint_h5(dir)) && isfile(checkpoint_json(dir))

"""
    stop_requested(::Checkpointer)

whether the simulation should stop now
"""
stop_requested(c::Checkpointer) =
    c.stopping || time() ≥ c.deadline || (!isempty(c.dir) && isfile(stop_file(c)))

"""
    checkpoint_due(::Checkpointer)

whether enough time has passed since the last checkpoint
"""
checkpoint_due(c::Checkpointer) = c.interval ≠ 0 && time() ≥ c.next

"""
    save_checkpoint(::Checkpointer, ::Simulation, state, time, sweep)

write a checkpoint recording the given state, simulation time and number of sweeps done in
the current phase. The state and the metadata are written to temporary files and moved into
place afterwards, so that a crash during the write leaves the previous checkpoint intact.
"""
function save_checkpoint(c::Checkpointer, sim, state::State, t::Number, sweep::Int)
    isempty(c.dir) && return nothing
    positions = Dict{String, Int}()
    for (name, f) in sim.files
        if f isa IO && f ∉ (stdout, stderr, devnull)
            flush(f)
            positions[name] = position(f)
        end
    end
    h5, js = checkpoint_h5(c.dir), checkpoint_json(c.dir)
    th5, tjs = h5 * ".tmp", js * ".tmp"
    rm(th5; force = true)
    save_state(th5, "checkpoint", state)
    open(tjs, "w") do io
        JSON.print(io, Dict(
            "version" => checkpoint_file_version,
            "id" => c.id,
            "phase" => c.phase,
            "sweep" => sweep,
            "time" => [real(c.phase_time), imag(c.phase_time)],
            "positions" => positions,
            "data" => sim.data,
        ))
    end
    mv(th5, h5; force = true)
    mv(tjs, js; force = true)
    c.next = time() + c.interval
    return nothing
end

"""
    load_checkpoint(dir)

read the checkpoint of the given directory and return `(state, phase_time, phase, sweep,
positions, data, id)`. `phase_time` is the simulation time at the start of the interrupted
phase, which is what the solvers count their sweeps from.
"""
function load_checkpoint(dir::String)
    meta = JSON.parsefile(checkpoint_json(dir))
    if meta["version"] ≠ checkpoint_file_version
        error("checkpoint of $dir has version $(meta["version"]), expected $checkpoint_file_version")
    end
    id = get(meta, "id", "")
    state = load_state(checkpoint_h5(dir), "checkpoint")
    re, im = meta["time"]
    t = im == 0 ? re : complex(re, im)
    positions = Dict{String, Int}(k => Int(v) for (k, v) in meta["positions"])
    data = Dict{String, Dict}(k => Dict(v) for (k, v) in meta["data"])
    return (state, t, Int(meta["phase"]), Int(meta["sweep"]), positions, data, id)
end

"""
    truncate_outputs(dir, positions)

cut the output files back to the length they had when the checkpoint was written, so that
the lines produced after it are not duplicated when the simulation resumes
"""
function truncate_outputs(dir::String, positions::Dict{String, Int})
    for (name, pos) in positions
        path = joinpath(dir, name)
        isfile(path) || continue
        if filesize(path) > pos
            open(path, "a") do io
                Base.truncate(io, pos)
            end
        end
    end
    return nothing
end

"""
    checkpoint_step!(::Checkpointer, sim, state, time, sweep)

record the progress of a solver, write a checkpoint if one is due or if the simulation is
about to stop, and return whether the loop should break out
"""
function checkpoint_step!(c::Checkpointer, sim, state::State, t::Number, sweep::Int)
    c.sweep = sweep
    c.simtime = t
    # kept so that an interrupt, which unwinds before the phase returns its state, still
    # has something to checkpoint
    c.state = state
    stop = stop_requested(c)
    if stop || checkpoint_due(c)
        save_checkpoint(c, sim, state, t, sweep)
    end
    c.stopping = stop
    return stop
end

"""
    first_sweep!(::Checkpointer)

the sweep a solver must start from, consuming the resume point so that the next phase
starts from the beginning
"""
function first_sweep!(c::Checkpointer)
    s = c.skip
    c.skip = 0
    return s + 1
end
