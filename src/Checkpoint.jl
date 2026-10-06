# Checkpointing a running simulation: when a checkpoint is due or a stop requested, and how a
# checkpoint is written, committed and loaded back.

"""
    checkpoint_file_version

the layout version of a checkpoint; `load_checkpoint` refuses a checkpoint of another one.
"""
const checkpoint_file_version = 4

############### commits ###############

"""
    Commit

a point a simulation can be resumed from, taken by `commit!` at the end of a sweep or at a
phase boundary, once everything the uninterrupted run writes before that point is written. A
checkpoint writes a commit and nothing else, so that it holds the state, the counts and the
outputs of one and the same moment.

- `phase`:      index of the phase to resume, one past the last when the run is over
- `sweep`:      sweeps of that phase done
- `phase_time`: simulation time the phase started from, which its solver counts sweeps from
- `time`:       simulation time reached
- `state`:      the state reached, `nothing` before the first phase has made one
- `carried`:    the value the phase carries from one step to the next, as the energy of the
                last dmrg sweep or the logarithm of the trace `thermal_state` has reached,
                see `run_steps`
- `reached`:    how far every destination had got, see `output_marks`
"""
struct Commit
    phase::Int
    sweep::Int
    phase_time::Number
    time::Number
    state::Union{Nothing, AbstractState}
    carried::Any
    reached::Dict{Union{String, Data}, Any}
end

"""
    Checkpointer(dir = ""; interval = 0, max_time = Inf)

the checkpointing machinery of a simulation, created by `runTMS` and carried by the
`Simulation`. A simulation stops cleanly when `max_time` seconds have passed and, with a
directory, when the file `stop` appears in it or on an interrupt; with a directory, a
checkpoint is written first.

# Fields

- `dir`:        the simulation directory, where the checkpoint is written, empty for none
- `interval`:   seconds between two checkpoints, zero or less for no periodic checkpoint
- `deadline`:   the `time()` after which the simulation stops cleanly, `Inf` for none
- `next`:       the `time()` of the next checkpoint
- `last`:       the last commit, which a checkpoint writes
- `resume`:     the commit a resumed run starts from, until the phase it belongs to has run
- `sweeps`:     whether the phase being run has read its resume point, which lets its sweeps
                be committed, see `resume_step`
- `written`:    the commit the checkpoint on the disk holds, which is not written again
- `generation`: which of the two state files the checkpoint on the disk names, 0 for none
- `stopping`:   set once a stop has been requested, so that every loop unwinds
"""
mutable struct Checkpointer
    dir::String
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

Checkpointer(dir::String = ""; interval::Real = 0, max_time::Real = Inf) =
    Checkpointer(dir, interval,
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

whether a periodic checkpoint is due, an interval of zero or less meaning none, as for
`sweep_due`
"""
checkpoint_due(c::Checkpointer) = c.interval > 0 && time() ≥ c.next

"""
    commit!(::Checkpointer, ::Outputs, phase, sweep, phase_time, time, state; carried)

record a point the simulation can be resumed from, see `Commit`, with how far the
destinations have got at that same moment
"""
commit!(c::Checkpointer, o::Outputs, phase::Int, sweep::Int, phase_time::Number, t::Number,
        state::Union{Nothing, AbstractState}; carried = nothing) =
    c.last = Commit(phase, sweep, phase_time, t, state, carried, output_marks(o))

"""
    check_carried(x, step)

refuse the value `x` carried by step `step` of `run_steps` when a checkpoint would not give it
back equal and of the same type: a symbol, a tuple, a float of another width, a dictionary
whose keys are not strings, an object
"""
function check_carried(x, step::Int)
    back = try
        Some(restored_value(JSON.parse(JSON.json(checkpoint_value(x)))))
    catch e
        if e isa InterruptException
            rethrow()
        end
        nothing
    end
    shown(v) = sprint(show, v; context = :limit => true)
    if isnothing(back)
        error("step $step of run_steps carries $(shown(x)), which a checkpoint cannot write: " *
              "carry numbers, strings, and vectors and dictionaries with string keys of them")
    elseif typeof(something(back)) ≠ typeof(x) || !isequal(something(back), x)
        error("step $step of run_steps carries $(shown(x)), which a checkpoint gives back as " *
              "$(shown(something(back))): carry numbers, strings, and vectors and dictionaries " *
              "with string keys of them")
    end
    return nothing
end

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
is renamed into place last, so that a crash at any point leaves the previous checkpoint or
this one, whole.

Nothing is written without a directory, before the first phase has made a state, or for a
commit already on the disk.
"""
function write_checkpoint(c::Checkpointer, o::Outputs)
    k = c.last
    if isempty(c.dir) || isnothing(k) || !(k.state isa AbstractState)
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
                "phase" => k.phase,
                "sweep" => k.sweep,
                "carried" => checkpoint_value(k.carried),
                # so that a complex time with no imaginary part stays complex
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
    load_checkpoint(dir[, system])

read the checkpoint of `dir`, as a named tuple of the fields of its commit (`phase`, `sweep`,
`phase_time`, `time`, `state`, `carried`), the `outputs` to put back with `restore_outputs!`
and the `generation` of its state file. The state comes back on the system
`system(phase, sites)` gives, or on a system of its own when that is `nothing`. A checkpoint
of another version, or whose state file is missing, is refused.
"""
function load_checkpoint(dir::String, system = (_, _) -> nothing)
    meta = JSON.parsefile(checkpoint_json(dir))
    if meta["version"] ≠ checkpoint_file_version
        error("checkpoint of $dir has version $(meta["version"]), expected $checkpoint_file_version")
    end
    file = meta["state"]
    path = joinpath(dir, file)
    if !isfile(path)
        error("the checkpoint of $dir names the state file $file, which is missing")
    end
    phase = Int(meta["phase"])
    return (phase, sweep = Int(meta["sweep"]),
            phase_time = restored_value(meta["phase_time"]), time = restored_value(meta["time"]),
            state = load_state(path, "checkpoint";
                               system = system(phase, saved_sites(path, "checkpoint"))),
            carried = restored_value(meta["carried"]),
            outputs = meta["outputs"],
            generation = file == state_file(2) ? 2 : 1)
end

"""
    resume_schedule(x, done)

the tail of a per sweep schedule, of `noise` or of the fields of a `Limits`, for a dmrg
resuming at sweep `done + 1`, a schedule that runs out being continued with its last value,
as ITensor does; a plain value is returned as it is. For `dmrg` alone: the evolution solvers
keep the sweep numbers across a resume, see `sweep_limits`.
"""
resume_schedule(x, ::Int) = x
resume_schedule(x::Vector, done::Int) =
    isempty(x) || done < length(x) ? x[done + 1:end] : x[end:end]
resume_schedule(l::Limits, done::Int) =
    Limits(resume_schedule(l.cutoff, done), resume_schedule(l.maxdim, done),
           resume_schedule(l.mindim, done))
