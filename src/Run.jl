# runTMS, which runs a simulation described by a SimData phase after phase, in a directory of its
# own, with its log, its checkpoints and the resumption of an interrupted run.

export runTMS, SimData

"""
    check_threading(threading)

refuse a `threading` of `SimData` other than `nothing`, `:auto`, `:dense` and `:blocks` when
the simulation is described, rather than when its phases start.
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
- `phases`: the phases of the simulation, see `Phases`, as a vector which may contain vectors
  to any depth and is flattened. The first phase must be `CreateState` or `LoadState`, the
  simulation having no state before it
- `description`: the text of the `description` file of the simulation (default `""`)
- `time_start`: the initial simulation time (default 0.)
- `final_measures`: the measurements to make at the end of the simulation, see `output`
  (default `[]`)
- `time_format`: the C like format of the simulation times written (default
  `$default_time_format`)
- `data_format`: the C like format of the measured values written (default
  `$default_data_format`)
- `checkpoint_interval`: seconds between two checkpoints (default 0, no periodic checkpoint;
  a stop or an interrupt still writes one, so that the simulation can be resumed)
- `max_time`: seconds after which the simulation stops cleanly (default `Inf`)
- `threading`: how the tensor contractions are threaded, see `set_threading`: `:dense`
  (default), `:blocks`, `:auto`, which chooses the mode before each phase from the system of
  the state, `:dense` until there is one, or `nothing`, which leaves the settings of the
  process as they are. The settings in force before the run are put back when `runTMS`
  returns

A checkpoint is written in the directory of the simulation, and `runTMS` resumes from it on
its own when it finds one. The simulation stops cleanly, writing a checkpoint, when
`max_time` is past, when the file `<name>/stop` appears, or on an interrupt. Run with the
`output` of `runTMS`, it has no directory: only `max_time` stops it, and no checkpoint is
written.

# Examples

    SimData(
        name = "ising",
        phases = [
            CreateState{Pure}(10, Qubit(), "Up"),
            GroundState(
                hamiltonian = -sum(Z(i)Z(i + 1) for i in 1:9) - sum(X(i) for i in 1:10),
                limits = Limits(maxdim = [10, 20, 50]), nsweeps = 10),
        ],
        final_measures = "data" => [X, Z],
    )
"""
@kwdef struct SimData
    description::String = ""
    name::String = "simulation"
    time_start::Number = 0.
    final_measures = []
    time_format::String = default_time_format
    data_format::String = default_data_format
    checkpoint_interval::Real = 0
    max_time::Real = Inf
    threading::Union{Nothing, Symbol} = :dense
    phases
    # the phases are flattened once, here, so that everything downstream works on a single
    # list: the phase loop, the position a checkpoint records, the fingerprint that tells
    # one simulation from another. None of them has to remember to do it, and none of them
    # can disagree on what the phases of a simulation are.
    SimData(description, name, time_start, final_measures, time_format, data_format,
            checkpoint_interval, max_time, threading, phases) =
        new(description, name, time_start, final_measures, time_format, data_format,
            checkpoint_interval, max_time, check_threading(threading),
            check_first_phase(flatten_phases(phases)))
end

"""
    flatten_phases(phases)

the phases given as nested vectors, as is convenient when a program builds them in pieces,
flattened into a single list by `SimData`. A phase then has one well defined position, which
is what a checkpoint records.
"""
flatten_phases(p::Vector) = reduce(vcat, map(flatten_phases, p); init = [])
flatten_phases(p) = [p]

"""
    check_first_phase(phases)

refuse a first phase other than `CreateState` or `LoadState`. A simulation starts without a
state, which every other phase transforms: it would fail on `nothing` deep inside its solver.
"""
function check_first_phase(phases::Vector)
    if isempty(phases) || !(first(phases) isa Union{CreateState, LoadState})
        error("the first phase must be CreateState or LoadState, which give the simulation its state")
    end
    return phases
end

# A `SimData` has the shape of a phase — `name`, `time_start`, `final_measures` — because
# `runTMS` runs the top level one through `log_phase` like any other phase, which is where
# the first line of the log comes from. That makes `run_phase(::Simulation, ::SimData)`
# reachable for a `SimData` sitting inside `phases`, and there it silently misbehaves: the
# loop it opens shares the phase counter of the loop around it, so it skips every phase
# whose index is below the one the outer loop had reached. Phases are grouped with plain
# vectors, which flatten properly, so this is refused rather than half supported.
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
        final_measures = $(s.final_measures),
        time_format = $(repr(s.time_format)),
        data_format = $(repr(s.data_format)),
        checkpoint_interval = $(s.checkpoint_interval),
        max_time = $(s.max_time),
        threading = $(repr(s.threading)),
        phases =
    $(s.phases))"""
    )

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
    runTMS(::SimData)
    runTMS(::SimData; clean = true)
    runTMS(::SimData; restart = true)
    runTMS(::SimData; output = myoutput)

run the given simulation, see `SimData`, in a directory named after it, and return the
`Simulation` it ends with. A checkpoint found in the directory is resumed from, and refused
if it belongs to a simulation with other phases.

- `restart` (default `false`): remove the simulation directory first
- `clean` (default `false`): remove the simulation directory and return without running
- `output`: a stream to redirect everything to, `stdout` or `devnull` for instance, instead of
  writing a directory

With a directory, an interrupt is a clean stop: a checkpoint is written and `runTMS` returns
instead of killing the program. For that, Ctrl-C raises `InterruptException` during the run,
a process wide setting that is given back the default Julia applies when `runTMS` returns.

`runTMS` is meant to be called once at a time in a process. It sets the threading of the
contractions for the whole process, unless its `SimData` has `threading = nothing`, and a
`CreateState` with a `seed` reseeds the global random generator. With a directory, it also
changes the working directory of the process and the Ctrl-C behaviour for the duration of the
run. Two simulations run at once in the same process, in parallel or one calling the other,
would fight over all of these: run them in separate processes. Passing `output` only spares
the working directory and the Ctrl-C behaviour.
"""
function runTMS(sim_data::SimData; restart::Bool=false, clean::Bool=false, output::Union{Nothing, IO} = nothing)
    live = isnothing(output)
    if live && (restart || clean)
        rm(sim_data.name; recursive = true, force = true)
    end
    if clean 
        return
    end
    start_dir = pwd()
    saved_threading = isnothing(sim_data.threading) ? nothing : ThreadingState()
    try
        # set at once, so that the stamp records the settings the run starts with: `:auto`
        # starts dense, having no state to choose from yet
        if !isnothing(sim_data.threading)
            set_threading(sim_data.threading == :auto ? :dense : sim_data.threading)
        end
        if live
            mkpath(sim_data.name);
            cd(sim_data.name);
            touch("running")
            # a stop left over from the previous run would stop this one immediately, and the
            # marker of a failed run would go on describing this one once it has succeeded
            rm("stop"; force = true)
            rm("error"; force = true)
            # scripts exit straight away on an interrupt, which would lose the state.
            # asking for an exception instead lets the simulation checkpoint and quit.
            Base.exit_on_sigint(false)
            if sim_data.description ≠ ""
                write("description", sim_data.description)
            end
            write("stamp", """
                    Julia $VERSION
                    TensorMixedStates $(pkgversion(TensorMixedStates))
                    Date $(now())
                    """ * threading_stamp(sim_data.threading))
            src_path = Base.source_path()
            if !isnothing(src_path) && src_path ≠ ""
                cp(src_path, "prog.jl"; force = true)
            end
        end
        c = Checkpointer(live ? "." : "", phases_id(sim_data.phases);
                         interval = sim_data.checkpoint_interval, max_time = sim_data.max_time)
        sim = Simulation(nothing; output, sim_data.time_format, sim_data.data_format, checkpoint = c)
        try
            if live && has_checkpoint(".")
                k = load_checkpoint(".")
                if k.id ≠ c.id
                    error("the checkpoint of \"$(sim_data.name)\" belongs to another simulation, " *
                          "its phases are not the ones being run. Use restart = true to start over " *
                          "and erase it, or choose another name.")
                end
                # put back before anything is written: a destination is created on first use,
                # which would empty a file the checkpoint continues
                restore_outputs!(sim.outputs, k.outputs)
                c.generation = k.generation
                log_msg(sim, "Resuming from checkpoint: phase $(k.phase), sweep $(k.sweep), simulation time $(k.time)")
                # the resume point is the last commit from the start, so that an interrupt
                # before the phase it belongs to has begun writes it back as it was, rather
                # than the state it holds as the start of that phase
                c.resume = Commit(k.phase, k.sweep, k.phase_time, k.time, k.state, k.energy,
                                  output_marks(sim.outputs))
                c.last = c.resume
            end
            try
                sim = log_phase(sim, sim_data)
            catch e
                # an interrupt is a request to stop cleanly, anything else is a real failure.
                # Without a directory nothing can be saved and nothing resumed, so the
                # interrupt goes on to the caller, as it would outside `runTMS`
                if !(e isa InterruptException) || isempty(c.dir)
                    rethrow()
                end
                log_msg(sim, "\n***** Interrupted, writing a checkpoint *****")
                # the last commit, whatever was written since: the checkpoint and the outputs
                # it resumes are those of one moment
                write_checkpoint(c, sim.outputs)
                k = c.last
                if !isnothing(k) && k.state isa AbstractState
                    # the returned simulation must carry what was reached, not what the phase
                    # was handed when it started
                    sim = Simulation(sim, k.state, k.time)
                end
                c.stopping = true
            end
        finally
            # a failing phase must not take away what was collected before it: the json
            # destinations are only written when the files are closed, so that has to
            # happen on the way out of an exception too
            close_sim_files(sim)
        end
        return sim
    catch
        # all the catch has of its own: the marker is written while the working directory
        # is still the simulation's, the finally below leaving it just afterwards
        if live
            touch("error")
        end
        rethrow()
    finally
        # leaving the simulation, by whichever way, is described here and nowhere else, so
        # that a step added later cannot be put on one path and forgotten on the other
        if live
            rm("running"; force = true)
            cd(start_dir)
            # the flag is process wide and would otherwise change how Ctrl-C behaves for
            # everything the caller runs afterwards. There is no way to read it back, so
            # what goes back is the default Julia itself applies: on in a script, off in
            # the REPL and in a session started with `-i`
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
        log_msg(sim, "Threading set to :$mode: BLAS threads $(after.blas), Strided threads " *
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
`time_start`, run by `run_phase` and measured by its `final_measures`, unless it stopped for
a checkpoint.
"""
function log_phase(sim::Simulation, phases::Vector; threading = nothing)
    c = sim.checkpoint
    r = c.resume
    # a resumed run starts again from the state and the time of its checkpoint, past the phases
    # it had completed. Put back here rather than where the checkpoint is read, since the
    # `time_start` of the simulation is applied in between, and a checkpoint written after the
    # last phase resumes none, the final measurements then still having to be taken at the
    # time the simulation reached
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
            # what a stopped run hands back is what it resumes from: the state and the time
            # of its last commit, the start of the phase when that one commits no sweep
            sim = Simulation(sim, c.last.state, c.last.time)
            log_stop(sim, i)
            break
        end
        # the next phase is the one to resume from, committed at a clean boundary. Marked here
        # also so that an interrupt falling after the last phase, in the final measurements,
        # still checkpoints what the simulation reached
        commit!(c, sim.outputs, i + 1, 0, sim.time, sim.time, sim.state)
        stop = stop_requested(c)
        if stop || checkpoint_due(c) || (i == length(phases) && (c.interval > 0 || c.generation ≠ 0))
            # the last phase done, a checkpoint records it whether one is due or not, so that
            # running the simulation again resumes past every phase and does nothing. Without
            # periodic checkpoints, only when one is on the disk, left by a stop: the next run
            # would resume from it and compute the end of the simulation again
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
    log_stop(sim, i)

log that the simulation stops after phase `i`, and whether it can be resumed.
"""
log_stop(sim::Simulation, i::Int) =
    log_msg(sim, isempty(sim.checkpoint.dir) ?
        "***** Stopping after phase $i, with no directory to save it in: it cannot be resumed *****" :
        "***** Stopping after phase $i, the simulation can be resumed *****")

"""
    check_is_phase(phase)

refuse an object without the three fields every phase is read through, `name`, `time_start`
and `final_measures`, so that it says so instead of failing with a `FieldError` in the middle
of a run. A missing `run_phase` method is then reported by the fallback of `run_phase`.
"""
function check_is_phase(phase)
    for f in (:name, :time_start, :final_measures)
        if hasfield(typeof(phase), f)
            continue
        end
        error("$(typeof(phase)) is not a phase, it has no `$f`. A phase needs name, " *
              "time_start, final_measures and a run_phase method")
    end
end

function log_phase(sim::Simulation, phase)
    check_is_phase(phase)
    log_msg(sim, "\n***** Starting phase \"$(phase.name)\" *****")
    if !isnothing(phase.time_start)
        sim = Simulation(sim, sim.state, phase.time_start)
    end
    td = @timed begin
        sim = run_phase(sim, phase)
        # a phase stopped for a checkpoint takes its final measurements when it is resumed and
        # finished, from the state an uninterrupted run takes them from
        if !sim.checkpoint.stopping
            output(sim, phase.final_measures)
        end
    end
    elapsed = round(td.time; digits=3)
    comp =
        if haskey(td, :compile_time)
            round(td.compile_time + td.recompile_time; digits=3)
        else
            nothing
        end
    log_msg(sim, "***** Ending phase \"$(phase.name)\" after $elapsed seconds, $(Base.format_bytes(td.bytes)) allocated *****")
    if !isnothing(comp)
        log_msg(sim, "compilation time was $comp seconds")
    end
    return sim
end

run_phase(sim::Simulation, sd::SimData) =
    log_phase(sim, sd.phases; sd.threading)
