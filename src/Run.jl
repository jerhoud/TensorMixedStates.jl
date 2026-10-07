# runTMS, which runs a simulation described by a SimData phase after phase, in a directory of its
# own, with its log, its checkpoints and the resumption of an interrupted run.

export runTMS, SimData, stopped

"""
    check_threading(threading)

refuse a `threading` of `SimData` other than `nothing`, `:auto`, `:dense` and `:blocks`
"""
function check_threading(threading)
    if !(threading in (nothing, :auto, :dense, :blocks))
        error("unknown threading $(repr(threading)), expected :auto, :dense, :blocks or nothing")
    end
    return threading
end

"""
    SimData(; name = "simulation", phases, options...)

the description of a simulation, which `runTMS` runs.

# Fields

- `name`: the name of the simulation, and of the directory its results are written to
- `phases`: the phases of the simulation, see `TensorMixedStates.AbstractPhase`, as a vector
  which may contain vectors to any depth and is flattened. The first phase must create the
  state, as `CreateState` and `LoadState` do, see `TensorMixedStates.creates_state`
- `description`: the text of the `description` file of the simulation (default `""`)
- `time_start`: the initial simulation time (default 0.)
- `final_measurements`: the measurements to make at the end of the simulation, see `output`
  (default `[]`)
- `time_format`: the C like format of the simulation times written (default
  `$default_time_format`)
- `data_format`: the C like format of the measured values written (default
  `$default_data_format`)
- `checkpoint_interval`: seconds between two checkpoints (default 0, no periodic checkpoint;
  a stop or an interrupt still writes one)
- `max_time`: seconds after which the simulation stops cleanly (default `Inf`)
- `threading`: how the tensor contractions are threaded, see `set_threading`: `:dense`
  (default), `:blocks`, `:auto`, chosen before each phase from the system of the state,
  `:dense` before there is one, or `nothing`, leaving the settings of the process as they are.
  The settings before the run are put back when `runTMS` returns

The simulation stops cleanly, writing a checkpoint in its directory that `runTMS` resumes
from, when `max_time` is past, when the file `<name>/stop` appears, or on an interrupt. Run
with the `output` of `runTMS`, it has no directory: only `max_time` stops it, and no
checkpoint is written.

# Examples

    SimData(
        name = "ising",
        phases = [
            CreateState{Pure}(10, Qubit(), "Up"),
            GroundState(
                hamiltonian = -sum(Z(i)Z(i + 1) for i in 1:9) - sum(X(i) for i in 1:10),
                limits = Limits(maxdim = [10, 20, 50]), nsweeps = 10),
        ],
        final_measurements = "data" => [X, Z],
    )
"""
@kwdef struct SimData
    description::String = ""
    name::String = "simulation"
    time_start::Number = 0.
    final_measurements = []
    time_format::String = default_time_format
    data_format::String = default_data_format
    checkpoint_interval::Real = 0
    max_time::Real = Inf
    threading::Union{Nothing, Symbol} = :dense
    phases
    # flattened once, here, so that the phase loop and the position a checkpoint records agree
    SimData(description, name, time_start, final_measurements, time_format, data_format,
            checkpoint_interval, max_time, threading, phases) =
        new(description, name, time_start, final_measurements, time_format, data_format,
            checkpoint_interval, max_time, check_threading(threading),
            check_first_phase(check_phases(flatten_phases(phases))))
end

"""
    flatten_phases(phases)

the phases given as nested vectors, flattened into a single list, where a phase has the
position a checkpoint records
"""
flatten_phases(p::Vector) = reduce(vcat, map(flatten_phases, p); init = [])
flatten_phases(p) = [p]

"""
    check_phases(phases)

refuse, when the simulation is written, an object among the phases that is not one, see
`check_is_phase`: a program corrected later could not resume its checkpoint
"""
function check_phases(phases::Vector)
    foreach(check_is_phase, phases)
    return phases
end

"""
    check_first_phase(phases)

refuse a first phase that does not create the state, see `creates_state`
"""
function check_first_phase(phases::Vector)
    if isempty(phases) || !creates_state(first(phases))
        error("the first phase must be CreateState, LoadState or a phase of one's own creating " *
              "the state, see TensorMixedStates.creates_state")
    end
    return phases
end

# A `SimData` has the fields of a phase, `runTMS` running it through `log_phase`, but within
# `phases` its loop would share the phase counter of the loop around it and skip phases.
flatten_phases(sd::SimData) =
    error("the SimData \"$(sd.name)\" cannot be a phase of another simulation, " *
          "nest plain vectors instead")

show(io::IO, s::SimData) =
    print(io,
    """
    SimData(
        description = $(repr(s.description)),
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measurements = $(s.final_measurements),
        time_format = $(repr(s.time_format)),
        data_format = $(repr(s.data_format)),
        checkpoint_interval = $(s.checkpoint_interval),
        max_time = $(s.max_time),
        threading = $(repr(s.threading)),
        phases =
    $(s.phases))"""
    )

"""
    same_program(dir, src_path)

whether the checkpoint of the directory `dir` may be resumed by the program `src_path`
given the arguments `ARGS`: the same, byte for byte, as `prog.jl`, given the arguments
`prog_args.json` holds, none when it is absent. A program with no file, run from the REPL, is
always accepted.
"""
function same_program(dir, src_path)
    if isnothing(src_path) || src_path == ""
        return true
    end
    prog, args_file = joinpath(dir, "prog.jl"), joinpath(dir, "prog_args.json")
    args = isfile(args_file) ? JSON.parsefile(args_file) : []
    return isfile(prog) && read(src_path) == read(prog) && args == ARGS
end

"""
    live_marker(text)

whether the file `running` holding `text`, the machine and the process of a run, marks a run
that may still be going on: not when it names a process of this machine that no longer exists,
or this very process, whose number the dead one had. A file naming another machine, or none, is
taken as live.
"""
function live_marker(text::AbstractString)
    parts = split(text)
    if length(parts) ≠ 2 || parts[1] ≠ gethostname() || isnothing(tryparse(Int32, parts[2]))
        return true
    end
    pid = parse(Int32, parts[2])
    # libuv answers signal 0 on every system: 0 for a process that exists, a refusal of
    # permission for one of another user, and ESRCH for none
    return pid ≠ getpid() && ccall(:uv_kill, Cint, (Cint, Cint), pid, 0) ≠ Base.UV_ESRCH
end

"""
    threading_stamp(mode)

the lines of the `stamp` file saying how the run is threaded when it starts: the `threading`
of its `SimData`, and the settings of the process that its running time depends on, see
`threading_settings`.
"""
function threading_stamp(mode)
    s = threading_settings()
    return """
        Threading $(repr(mode))
        BLAS $(s.blas_library), $(s.blas) threads
        Julia threads $(s.julia)
        GC threads $(s.gc)
        Strided threads $(s.strided)
        Block sparse multithreading $(s.blocksparse ? "on" : "off")
        CPU threads $(Sys.CPU_THREADS)
        """
end

"""
    resume_system(phases, phase, sites)

the system the state of a checkpoint taken at the start of the phase `phase`, or in it, is put
back on: that of the last phase creating the state before it, see `phase_system`, when it has
the `sites` of the state, `nothing` otherwise, the state then coming back on a system of its
own
"""
function resume_system(phases::Vector, phase::Int, sites)
    i = findlast(creates_state, phases[1:min(phase - 1, end)])
    if isnothing(i)
        return nothing
    end
    system = phase_system(phases[i])
    return !isnothing(system) && system.sites == sites ? system : nothing
end

"""
    stopped(sim)

whether the simulation `runTMS` returned stopped before the end of its phases, at `max_time`,
on the file `stop` or on an interrupt: it is resumed by running it again, see
[High Level Interface](@ref). Within a phase of one's own, whether the run is stopping, its
solver or its steps having stopped for a checkpoint: the phase then leaves out what it writes
at its end, which its resume writes.

# Examples

    sim = runTMS(sim_data)
    if stopped(sim)
        println("stopped at time ", sim.time, ", to be resumed")
    end
"""
stopped(sim::Simulation) = sim.checkpoint.stopping

"""
    runTMS(::SimData)
    runTMS(::SimData; clean = true)
    runTMS(::SimData; restart = true)
    runTMS(::SimData; output = myoutput)

run the given simulation, see `SimData`, in a directory named after it, and return the
`Simulation` it ends with. A checkpoint found in the directory is resumed from, and refused
if it was written by another program, or one given other arguments, see [Resuming](@ref).

- `restart` (default `false`): remove the simulation directory first
- `clean` (default `false`): remove the simulation directory and return without running
- `output`: a stream to redirect everything to, `stdout` or `devnull` for instance, instead of
  writing a directory

With a directory, an interrupt is a clean stop: a checkpoint is written and `runTMS` returns.
For that, Ctrl-C raises `InterruptException` during the run, a process wide setting put back
to the default of Julia when `runTMS` returns.

`runTMS` is meant to run one simulation at a time in a process: it sets for the whole process
the threading of the contractions (unless `threading = nothing`), the global random generator
(for a `CreateState` with a `seed`) and, with a directory, the behaviour of Ctrl-C. Run
simultaneous simulations, in parallel or one calling the other, in separate processes. The
working directory of the process is left as it is.
"""
function runTMS(sim_data::SimData; restart::Bool=false, clean::Bool=false, output::Union{Nothing, IO} = nothing)
    live = isnothing(output)
    if live && (restart || clean)
        name = sim_data.name
        # rm would empty the current directory, the program included, before failing; a link
        # is removed without its target
        if ispath(name) && !islink(name) &&
           startswith(joinpath(realpath(pwd()), ""), joinpath(realpath(name), ""))
            error("cannot remove \"$name\", which contains the current directory")
        end
        rm(name; recursive = true, force = true)
    end
    if clean 
        return
    end
    saved_threading = isnothing(sim_data.threading) ? nothing : ThreadingState()
    # whether the run has taken its directory: only then does the way out write `error` and
    # remove `running`
    started = false
    # the empty name would be the working directory itself
    if live && isempty(sim_data.name)
        error("a simulation run in a directory needs a name, that of its directory")
    end
    dir = live ? abspath(sim_data.name) : ""
    file(name) = joinpath(dir, name)
    try
        # set at once, for the stamp: `:auto` starts dense, with no state to choose from
        if !isnothing(sim_data.threading)
            set_threading(sim_data.threading == :auto ? :dense : sim_data.threading)
        end
        c = Checkpointer(dir;
                         interval = sim_data.checkpoint_interval, max_time = sim_data.max_time)
        # the machine and the process named in the file running of a run that ended without
        # removing it
        stale = nothing
        if live
            mkpath(dir)
            # refused before anything is written or loaded, the directory holding another run
            if isfile(file("running"))
                marker = read(file("running"), String)
                if live_marker(marker)
                    error("\"$(sim_data.name)\" holds the file running of another run, which " *
                          "may still be going on: if none is, remove $(sim_data.name)/running")
                end
                stale = split(marker)
            end
            src_path = Base.source_path()
            if has_checkpoint(dir) && !same_program(dir, src_path)
                error("the checkpoint of \"$(sim_data.name)\" was written by another program, " *
                      "or one given other arguments: to resume with this one, copy it onto " *
                      "prog.jl and its arguments into prog_args.json, or start over with " *
                      "restart = true")
            end
            started = true
            write(file("running"), "$(gethostname()) $(getpid())\n")
            # left over from the previous run
            rm(file("stop"); force = true)
            rm(file("error"); force = true)
            # an exception rather than an exit, so that the run writes its checkpoint
            Base.exit_on_sigint(false)
            if sim_data.description ≠ ""
                write(file("description"), sim_data.description)
            end
            write(file("stamp"), """
                    Julia $VERSION
                    TensorMixedStates $(pkgversion(TensorMixedStates))
                    Date $(now())
                    """ * threading_stamp(sim_data.threading))
            # cp refuses to copy the program onto itself, when it is the copy that is run
            if !isnothing(src_path) && src_path ≠ ""
                if !(isfile(file("prog.jl")) && samefile(src_path, file("prog.jl")))
                    cp(src_path, file("prog.jl"); force = true)
                end
                # no arguments, no file, as `same_program` reads it
                if isempty(ARGS)
                    rm(file("prog_args.json"); force = true)
                else
                    write(file("prog_args.json"), JSON.json(ARGS))
                end
            end
        end
        sim = Simulation(nothing; output, sim_data.time_format, sim_data.data_format, checkpoint = c)
        if !isnothing(stale)
            log_message(sim, "Warning: the run before, process $(stale[2]) on $(stale[1]), " *
                             "ended without removing running")
        end
        try
            if live && has_checkpoint(dir)
                system(phase, sites) = resume_system(sim_data.phases, phase, sites)
                k = load_checkpoint(dir, system)
                # before anything is written: a destination created on first use would empty a
                # file the checkpoint continues
                shortened = restore_outputs!(sim.outputs, k.outputs)
                c.generation = k.generation
                log_message(sim, "Resuming from checkpoint: phase $(k.phase), sweep $(k.sweep), simulation time $(k.time)")
                for name in shortened
                    log_message(sim, "Warning: $name is shorter than at the checkpoint, the lines " *
                                     "it lost are not written again")
                end
                # also the last commit, so that an interrupt before its phase has begun writes
                # it back as it was
                c.resume = Commit(k.phase, k.sweep, k.phase_time, k.time, k.state, k.carried,
                                  output_marks(sim.outputs))
                c.last = c.resume
            end
            try
                sim = log_phase(sim, sim_data)
            catch e
                # an interrupt is a clean stop, unless there is no directory to save in
                if !(e isa InterruptException) || isempty(c.dir)
                    rethrow()
                end
                log_message(sim, "\n***** Interrupted, writing a checkpoint *****")
                # the last commit, whatever was written since, so that the checkpoint and the
                # outputs it resumes are those of one moment
                write_checkpoint(c, sim.outputs)
                k = c.last
                if !isnothing(k) && k.state isa AbstractState
                    # what was reached, not what the phase was handed
                    sim = Simulation(sim, k.state, k.time)
                end
                c.stopping = true
            end
        finally
            # on a failure too: the json destinations are only written when the files are
            # closed
            close_sim_files(sim)
        end
        return sim
    catch
        if started
            touch(file("error"))
        end
        rethrow()
    finally
        # what the run set up is undone here only, whichever way it is left
        if started
            rm(file("running"); force = true)
            # process wide and impossible to read back: the default of Julia goes back, on in
            # a script, off in the REPL and with `-i`
            Base.exit_on_sigint(!isinteractive())
        end
        # process wide as well, and read back before the run
        if !isnothing(saved_threading)
            set_threading(saved_threading)
        end
    end
end

"""
    adapt_threading(sim, threading)

set the `threading` of a `SimData` before a phase: `:dense` or `:blocks` as asked, or for
`:auto` the mode suited to the system of the state, once there is one. A change is logged.
"""
function adapt_threading(sim::Simulation, threading)
    mode = threading == :auto ?
        (sim.state isa AbstractState ? threading_mode(sim.state.system) : nothing) : threading
    if isnothing(mode)
        return
    end
    before = threading_settings()
    set_threading(mode)
    after = threading_settings()
    if after != before
        log_message(sim, "Threading set to :$mode: BLAS threads $(after.blas), Strided threads " *
                         "$(after.strided), block sparse multithreading " *
                         (after.blocksparse ? "on" : "off"))
    end
end

"""
    log_phase(sim, phases::Vector; threading)
    log_phase(sim, phase)

run a list of phases, from the one a resumed run starts at, committing each boundary and
writing a checkpoint when one is due or a stop is asked for, which ends the loop, and setting
the `threading` of the `SimData` before each phase. A single phase is logged, given its
`time_start`, run by `run_phase` and measured by its `final_measurements`, unless it stopped for
a checkpoint.
"""
function log_phase(sim::Simulation, phases::Vector; threading = nothing)
    c = sim.checkpoint
    r = c.resume
    # the state and time of a checkpoint are put back here, after the `time_start` of the
    # simulation is applied: a checkpoint written after the last phase resumes none, and the
    # final measurements take the time it reached
    if isnothing(r)
        commit!(c, sim.outputs, 1, 0, sim.time, sim.time, sim.state)
    else
        sim = Simulation(sim, r.state, r.phase_time)
        c.last = r
    end
    for i in c.last.phase:length(phases)
        # a phase commits its sweeps only once it has read its resume point
        c.sweeps = false
        adapt_threading(sim, threading)
        sim = log_phase(sim, phases[i])
        # consumed by the phase it belongs to, whether it read it or not
        c.resume = nothing
        if c.stopping
            # a stopped run hands back what it resumes from, its last commit
            sim = Simulation(sim, c.last.state, c.last.time)
            log_stop(sim, i; within = true)
            break
        end
        # committed after the last phase too, for an interrupt in the final measurements
        commit!(c, sim.outputs, i + 1, 0, sim.time, sim.time, sim.state)
        stop = stop_requested(c)
        if stop || checkpoint_due(c) || (i == length(phases) && (c.interval > 0 || c.generation ≠ 0))
            # after the last phase, so that a new run does nothing; without periodic
            # checkpoints, only when one is on the disk, which a new run would resume
            write_checkpoint(c, sim.outputs)
        end
        if stop
            c.stopping = true
            log_stop(sim, i)
            break
        end
    end
    return sim
end

"""
    log_stop(sim, i; within = false)

log that the simulation stops after phase `i`, or in it when `within`, and whether it can be
resumed.
"""
function log_stop(sim::Simulation, i::Int; within::Bool = false)
    at = within ? "in" : "after"
    log_message(sim, isempty(sim.checkpoint.dir) ?
        "***** Stopping $at phase $i, with no directory to save it in: it cannot be resumed *****" :
        "***** Stopping $at phase $i, the simulation can be resumed *****")
end

"""
    check_is_phase(phase)

refuse an object that is not a subtype of `AbstractPhase`, or lacks one of the fields `name`,
`time_start` and `final_measurements`; a missing method is reported by the fallback of
`run_phase`
"""
function check_is_phase(phase)
    if !(phase isa AbstractPhase)
        error("$(typeof(phase)) is not a phase: a phase is a subtype of " *
              "TensorMixedStates.AbstractPhase, with the fields name, time_start and " *
              "final_measurements and a method of TensorMixedStates.run_phase")
    end
    for f in (:name, :time_start, :final_measurements)
        if hasfield(typeof(phase), f)
            continue
        end
        error("$(typeof(phase)) is not a phase, it has no `$f`. A phase needs name, " *
              "time_start, final_measurements and a run_phase method")
    end
end

function log_phase(sim::Simulation, phase)
    log_message(sim, "\n***** Starting phase \"$(phase.name)\" *****")
    if !isnothing(phase.time_start)
        sim = Simulation(sim, sim.state, phase.time_start)
    end
    td = @timed begin
        sim = run_phase(sim, phase)
        # a phase stopped for a checkpoint takes its final measurements once resumed
        if !sim.checkpoint.stopping
            output(sim, phase.final_measurements)
        end
    end
    elapsed = round(td.time; digits=3)
    comp =
        if haskey(td, :compile_time)
            round(td.compile_time + td.recompile_time; digits=3)
        else
            nothing
        end
    # a phase stopped for a checkpoint has not ended, and is resumed
    ends = sim.checkpoint.stopping ? "Stopping" : "Ending"
    log_message(sim, "***** $ends phase \"$(phase.name)\" after $elapsed seconds, $(Base.format_bytes(td.bytes)) allocated *****")
    if !isnothing(comp)
        log_message(sim, "compilation time was $comp seconds")
    end
    return sim
end

run_phase(sim::Simulation, sd::SimData) =
    log_phase(sim, sd.phases; sd.threading)
