# Checkpointing a running simulation: the commit, a point it can be resumed from, when a
# checkpoint is due or a stop requested, how a commit is written, in two generations so that a
# crash always leaves one whole, and how it is read back.

"""
    checkpoint_file_version

the layout version of a checkpoint; `load_checkpoint` refuses a checkpoint of another one
"""
const checkpoint_file_version = 5

"""
    Commit(phase, step, last, phase_time, time, state, values, marks)

a point a simulation can be resumed from, taken once everything the uninterrupted run writes
before that point is written: a checkpoint writes a commit and nothing else, so that it holds the
state, the counts and the outputs of one moment.

- `phase`: the index of the phase to run or to go on with, one past the last once they are done
- `step`: the steps of that phase done
- `last`: whether the phase ended with that step, the commit holding its result
- `phase_time`: the time that phase started from, from which its steps count
- `time`, `state`: the time and the state reached, `nothing` before a phase has made one
- `values`: the values of the simulation, its named states and its generator, see `Simulation`
- `marks`: how far each destination had got, see `resume_marks`
"""
struct Commit
    phase::Int
    step::Int
    last::Bool
    phase_time::Number
    time::Number
    state::Union{Nothing, AbstractState}
    values::Dict{Symbol, Any}
    marks::Dict{Destination, Any}
end

"""
    committed_values(values)

the values of a simulation as a commit keeps them: a dictionary of its own, and a copy of the
generator, which the steps after the commit advance
"""
committed_values(values::AbstractDict) =
    Dict{Symbol, Any}(k => v isa AbstractRNG ? copy(v) : v for (k, v) in values)

"""
    Checkpointer(dir, seed; interval = 0, max_time = Inf)

the checkpointing of a run: its directory, empty for a run without one, which writes no
checkpoint, the seed of its generator, which a checkpoint keeps, the seconds between two
checkpoints, zero for none but those of a stop, and the deadline after which the run stops; the
last commit, the commit the checkpoint on the disk holds and which of the two generations of
state files it names, 0 for none
"""
mutable struct Checkpointer
    dir::String
    seed::Integer
    interval::Float64
    deadline::Float64
    next::Float64
    last::Union{Nothing, Commit}
    written::Union{Nothing, Commit}
    generation::Int
end

Checkpointer(dir::String, seed::Integer; interval::Real = 0, max_time::Real = Inf) =
    Checkpointer(dir, seed, interval, max_time == Inf ? Inf : time() + max_time,
                 interval ≤ 0 ? Inf : time() + interval, nothing, nothing, 0)

"""
    stop_requested(c)

whether the run should stop now: its deadline is past, or the file `stop` is in its directory
"""
stop_requested(c::Checkpointer) =
    time() ≥ c.deadline || (!isempty(c.dir) && isfile(joinpath(c.dir, "stop")))

"""
    checkpoint_due(c)

whether a periodic checkpoint is due
"""
checkpoint_due(c::Checkpointer) = c.interval > 0 && time() ≥ c.next

"""
    checkpoint_json(dir)
    state_file(generation)

the file of a checkpoint that names its generation and the state file of a generation, 1 or 2
"""
checkpoint_json(dir::String) = joinpath(dir, "checkpoint.json")

state_file(g::Int) = "checkpoint-$g.h5"

"""
    checkpoint_files

the files a checkpoint writes in the directory of a simulation
"""
const checkpoint_files = ("checkpoint.json", "checkpoint.json.tmp", state_file(1), state_file(2))

"""
    relative(d, dir)
    rooted(d, dir)

the destination `d` with the path of a file of the directory `dir` made relative to it, and
back, so that a checkpoint still resumes once the directory has moved
"""
relative(d::TextFile, dir::String) = TextFile(relative(d.path, dir), d.format)
relative(d::JsonFile, dir::String) = JsonFile(relative(d.path, dir))
relative(path::String, dir::String) =
    startswith(path, joinpath(dir, "")) ? relpath(path, dir) : path
relative(d, _) = d

rooted(d::TextFile, dir::String) = TextFile(rooted(d.path, dir), d.format)
rooted(d::JsonFile, dir::String) = JsonFile(rooted(d.path, dir))
rooted(path::String, dir::String) = isabspath(path) ? path : normpath(joinpath(dir, path))
rooted(d, _) = d

"""
    state_name(name)

the name in a state file of the value `name` that is a state
"""
state_name(name::Symbol) = "value $name"

"""
    sync_file(path)

wait until what was written in the file `path` is on the disk, which the system may otherwise
keep in memory for a while after the file is closed, for a power cut to lose. Julia has no
function for it, `fsync` is called; nothing is done on Windows.
"""
function sync_file(path::String)
    if !Sys.iswindows()
        open(path) do io
            systemerror("fsync of $path", ccall(:fsync, Cint, (Cint,), fd(io)) ≠ 0)
        end
    end
    return nothing
end

"""
    write_checkpoint(c)

write the last commit down. Its state and the values that are states go to whichever of
`checkpoint-1.h5` and `checkpoint-2.h5` the checkpoint on the disk does not name, with the rest
of the commit serialized, and `checkpoint.json`, which names it, is renamed into place last, so
that a crash at any point leaves the previous checkpoint or this one, whole. Both are synced to
the disk before the rename, so that a power cut does too, but on Windows. Nothing is written
without a directory or for a commit already on the disk.
"""
function write_checkpoint(c::Checkpointer)
    k = c.last
    if isempty(c.dir) || isnothing(k) || k === c.written
        c.next = time() + c.interval
        return nothing
    end
    g = c.generation == 1 ? 2 : 1
    file = joinpath(c.dir, state_file(g))
    rm(file; force = true)
    # none before the first state, or after a CopyState that moved it away
    if !isnothing(k.state)
        save_state(file, "state", k.state)
    end
    rest = Dict{Symbol, Any}()
    for (name, v) in k.values
        if v isa AbstractState
            save_state(file, state_name(name), v)
        else
            rest[name] = v
        end
    end
    io = IOBuffer()
    marks = Dict{Destination, Any}(relative(d, c.dir) => kept(m) for (d, m) in k.marks)
    serialize(io, (; c.seed, k.phase, k.step, k.last, k.phase_time, k.time, values = rest,
                   states = [ name for (name, v) in k.values if v isa AbstractState ], marks))
    # the file made here when no state was saved there
    h5open(file, "cw") do f
        f["commit"] = take!(io)
    end
    tmp = checkpoint_json(c.dir) * ".tmp"
    open(tmp, "w") do f
        JSON.print(f, Dict("version" => checkpoint_file_version, "julia" => string(VERSION),
                           "seed" => string(c.seed), "phase" => k.phase, "step" => k.step,
                           "time" => string(k.time),
                           "state" => state_file(g)))
    end
    # on the disk before checkpoint.json names them: a power cut could otherwise leave it naming
    # a state file the system had not written yet, or empty, the other generation removed
    sync_file(file)
    sync_file(tmp)
    # at once, where the mv of Julia 1.10 removes checkpoint.json before renaming
    Base.Filesystem.rename(tmp, checkpoint_json(c.dir))
    # before the other generation goes: a checkpoint written again after an interrupt there goes
    # over that one, not over the one checkpoint.json names
    c.generation = g
    c.written = k
    rm(joinpath(c.dir, state_file(3 - g)); force = true)
    c.next = time() + c.interval
    return nothing
end

"""
    has_checkpoint(dir)

whether the directory holds a checkpoint
"""
has_checkpoint(dir::String) = isfile(checkpoint_json(dir))

"""
    load_checkpoint(dir)

the commit of the checkpoint of `dir`, the generation of its state file and the seed of the
generator of the run. A checkpoint of
another layout, or written by another minor version of Julia, whose serialization may differ,
is refused.
"""
function load_checkpoint(dir::String)
    meta = JSON.parsefile(checkpoint_json(dir))
    if meta["version"] ≠ checkpoint_file_version
        error("the checkpoint of $dir has the layout $(meta["version"]) of another version of " *
              "TensorMixedStates: start the simulation over")
    end
    v = VersionNumber(meta["julia"])
    if (v.major, v.minor) ≠ (VERSION.major, VERSION.minor)
        error("the checkpoint of $dir was written by Julia $v: resume it with Julia " *
              "$(v.major).$(v.minor), or start the simulation over")
    end
    file = joinpath(dir, meta["state"])
    if !isfile(file)
        error("the checkpoint of $dir names the state file $(meta["state"]), which is missing")
    end
    k = h5open(f -> deserialize(IOBuffer(read(f["commit"]))), file, "r")
    values = Dict{Symbol, Any}(k.values)
    for name in k.states
        values[name] = load_state(file, state_name(name))
    end
    marks = Dict{Destination, Any}(rooted(d, dir) => m for (d, m) in k.marks)
    state = h5open(f -> haskey(f, "state"), file, "r") ? load_state(file, "state") : nothing
    commit = Commit(k.phase, k.step, k.last, k.phase_time, k.time, state, values, marks)
    return (commit, meta["state"] == state_file(2) ? 2 : 1, k.seed)
end
