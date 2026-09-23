export runTMS, SimData

"""
    SimData(name = "my_simulation", phases::Vector{Phases} = [phase1, phase2...])

A type for describing a simulation to use with `runTMS`

# Fields

- `name`:            the name of the simulation used as the name of the directory to store the results
- `phases`:          the list of phases of the simulation (see Phases for a list of possible values),
  which may itself contain lists, to any depth, and is flattened on construction
- `description`:     text put in the description file of the simulation (default "")
- `time_start`:      initial simulation time (default 0.)
- `final_measures`:  measures to make at the end of simulation (default []) see `measure` and `output`
- `time_format`:     C like format for output of simulation time (default `$default_time_format`)
- `data_format`:     C like format for output of simulation data (default `$default_data_format`)
- `checkpoint_interval`: seconds between two checkpoints (default 0, no checkpointing)
- `max_time`:        seconds after which the simulation stops cleanly (default `Inf`)

A simulation with a checkpoint interval writes its state to `<name>/checkpoint.h5` and
`runTMS` resumes from it on its own if it finds one. It stops cleanly, after writing a
checkpoint, when `max_time` is past, when the file `<name>/stop` appears, or on an
interrupt.
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
    phases
    # the phases are flattened once, here, so that everything downstream works on a single
    # list: the phase loop, the position a checkpoint records, the fingerprint that tells
    # one simulation from another. None of them has to remember to do it, and none of them
    # can disagree on what the phases of a simulation are.
    SimData(description, name, time_start, final_measures, time_format, data_format,
            checkpoint_interval, max_time, phases) =
        new(description, name, time_start, final_measures, time_format, data_format,
            checkpoint_interval, max_time, flatten_phases(phases))
end

"""
    flatten_phases(phases)

phases may be given as nested vectors, for convenience when a program builds its phases in
pieces, and `SimData` flattens them into a single list. A phase then has one well defined
position, which is what a checkpoint records.
"""
flatten_phases(p::Vector) = reduce(vcat, map(flatten_phases, p); init = [])
flatten_phases(p) = [p]

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
        phases =
    $(s.phases))"""
    )


"""
    runTMS(::SimData)
    runTMS(::SimData; clean = true)
    runTMS(::SimData; restart = true)
    runTMS(::SimData; output = myoutput)

run the given simulation (see SimData for details), write the output to file and return a Simulation object containing the result.
`clean` (default `false`) remove the simulation directory and exit,
`restart` (default `false`) remove the simulation directory and run the simulation,
`output` redirect all output to the given IO channel (no output directory created), useful values are stdout or devnull (to suppress all output).

A simulation writing to a directory turns an interrupt into a clean stop: it writes a
checkpoint and returns, instead of killing the program. See `SimData` for the checkpointing
options. This asks the runtime to raise `InterruptException` on Ctrl-C, a process wide
setting that is put back when `runTMS` returns.

`runTMS` drives a whole process and is meant to be called once at a time. Writing to a
directory, it changes the working directory of the process for the duration of the run, and
it sets the Ctrl-C behaviour; a `CreateState` phase given a `seed` also reseeds the global
random generator. Two simulations running at once in the same process, whether in parallel
or through one calling the other, would fight over all three. Run them in separate
processes, or pass `output` so that nothing touches the working directory. Parallelism
inside a single simulation is a different matter and works as usual: it comes from the
threads ITensor uses for its contractions.

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
    try
        if live
            mkpath(sim_data.name);
            cd(sim_data.name);
            touch("running")
            # a stop left over from the previous run would stop this one immediately
            rm("stop"; force = true)
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
                    """)
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
                state, phase_time, phase, sweep, positions, data, json, id = load_checkpoint(".")
                if id ≠ c.id
                    error("the checkpoint of \"$(sim_data.name)\" belongs to another simulation, " *
                          "its phases are not the ones being run. Use restart = true to start over " *
                          "and erase it, or choose another name.")
                end
                truncate_outputs(".", positions)
                merge!(sim.data, data)
                # put back before any measurement asks for them: `get_sim_file` creates a
                # json destination on first use and would otherwise start an empty one
                merge!(sim.files, json)
                c.phase, c.skip, c.phase_time = phase, sweep, phase_time
                c.appending = c.resuming = true
                sim = Simulation(sim, state, phase_time)
                log_msg(sim, "Resuming from checkpoint: phase $phase, sweep $sweep, simulation time $phase_time")
            end
            try
                sim = log_phase(sim, sim_data)
            catch e
                # an interrupt is a request to stop cleanly, anything else is a real failure
                if !(e isa InterruptException)
                    rethrow()
                end
                log_msg(sim, "\n***** Interrupted, writing a checkpoint *****")
                # the state of the interrupted sweep, the one the phase never got to return
                st = c.state isa State ? c.state : sim.state
                if st isa State
                    save_checkpoint(c, sim, st, c.sweep)
                    # the returned simulation must carry what was reached, not what the phase
                    # was handed when it started
                    sim = Simulation(sim, st, c.simtime)
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
    end
end

function log_phase(sim::Simulation, phases::Vector)
    c = sim.checkpoint
    for (i, phase) in enumerate(phases)
        # phases already completed before the checkpoint are not replayed
        if i < c.phase
            continue
        end
        c.phase = i
        if c.resuming
            # the interrupted phase restarts from the time it began with, its solver
            # counts sweeps from there and skips the ones already done
            sim = Simulation(sim, sim.state, c.phase_time)
            c.resuming = false
        else
            c.phase_time = sim.time
        end
        phase_start!(c, sim)
        sim = log_phase(sim, phase)
        if c.stopping
            log_msg(sim, "***** Stopping after phase $i, the simulation can be resumed *****")
            break
        end
        # the next phase is the one to resume from, record it at a clean boundary
        c.phase = i + 1
        c.phase_time = sim.time
        # a resume point belongs to the phase it was written for: a phase with no sweeps of
        # its own must not inherit the ones the previous phase was told to skip
        c.skip = 0
        # marked again here so that an interrupt falling after the last phase, in the final
        # measurements, still checkpoints what the simulation reached
        phase_start!(c, sim)
        if sim.state isa State && (checkpoint_due(c) || stop_requested(c))
            if checkpoint_step!(c, sim, sim.state, sim.time, 0) 
                log_msg(sim, "***** Stopping after phase $i, the simulation can be resumed *****")
                break
            end
        end
    end
    return sim
end

# The three fields every phase is read through, here rather than at the first `phase.name`
# so that something which is not a phase says so instead of surfacing as a `FieldError` from
# the middle of a run. What is missing afterwards is a `run_phase` method, and its fallback
# in `Phases.jl` says that in its turn.
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
        output(sim, phase.final_measures)
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
    log_phase(sim, sd.phases)
