# The phases of a simulation, those of the library (runner/phases) and those of one's own: their
# interface, a phase made of steps or of a single call, and the end of each step, where the run
# measures, logs, commits, writes a checkpoint and stops.

"""
    abstract type AbstractPhase

the supertype of the phases of a simulation. Every phase has these three fields:

- `name`: the name of the phase, written in the log
- `time_start`: the time the clock is set to when the phase starts (default `nothing`, keeping
  the time)
- `final_measurements`: the measurements to make at the end of the phase, see `output`
  (default `[]`)

A phase of one's own is a subtype with these fields and either of two interfaces, `runTMS`
logging it, applying its `time_start` and making its final measurements:

- a phase made of steps has methods of `TensorMixedStates.prepare_phase`,
  `TensorMixedStates.step_count` and `TensorMixedStates.run_step`; the run measures its field
  `measurements`, if it has one, every `measurements_period` steps (1 without that field), logs,
  commits and stops between two steps, and resumes it after the last step committed
- a phase of a single call has a method of `TensorMixedStates.run_phase`, which calls
  `TensorMixedStates.step_done` wherever it can stop

A phase providing values to the measurements, as `:energy`, declares them with a method of
`TensorMixedStates.provided_values`.
"""
abstract type AbstractPhase end

# the fields are read from the type, so that those of a phase of one's own are printed too
function show(io::IO, s::AbstractPhase)
    t = typeof(s)
    print(io, "\n", nameof(t))
    if !isempty(t.parameters)
        # a parameter may be a value, which has no `nameof`
        print(io, "{", join((p isa Union{DataType, UnionAll} ? nameof(p) : repr(p)
                             for p in t.parameters), ", "), "}")
    end
    print(io, "(")
    fs = fieldnames(t)
    for (i, f) in enumerate(fs)
        print(io, "\n    ", f, " = ", repr(getfield(s, f)), i < length(fs) ? "," : ")")
    end
end

"""
    TensorMixedStates.provided_values(phase)

the names of the values `phase` gives the simulation, `(:energy,)` for `GroundState` for
instance, none by default. A measurement of a value that no phase provides is refused when the
`SimData` is written, and one that is not there yet is `missing`.
"""
provided_values(_) = ()

"""
    phase_files(phase)

the files `phase` writes in the directory of the simulation besides its measurements, the file of
a `SaveState`, none by default, checked before the run, see `check_destinations`
"""
phase_files(_) = ()

"""
    TensorMixedStates.prepare_phase(phase, sim)
    TensorMixedStates.step_count(prepared)
    TensorMixedStates.run_step(prepared, sim, k)

the interface of a phase made of steps: `prepare_phase` returns what the steps need, from the
simulation as the phase starts, `step_count` the number of steps, and `run_step` makes step `k`,
returning the simulation it leaves, see `Simulation`. On a resume, `prepare_phase` is given
the state of the last step committed, with the time the phase started from, and the steps go on
after it: what it prepares may depend on the system of the state, not on the state.

# Examples

    struct Kicks <: AbstractPhase
        name::String
        time_start
        final_measurements
        nkicks::Int
        measurements
    end

    TensorMixedStates.prepare_phase(p::Kicks, _) = p
    TensorMixedStates.step_count(p::Kicks) = p.nkicks
    TensorMixedStates.run_step(::Kicks, sim, _) =
        Simulation(sim, apply(exp(-0.3im * X)(1), sim.state), sim.time + 0.1)
"""
function prepare_phase end
function step_count end
function run_step end

"""
    step_due(phase, k)

whether the measurements of `phase` during its steps, its field `measurements`, are due at step
`k`, every `measurements_period` steps, 1 when the phase does not have that field
"""
function step_due(phase, k::Int)
    period = hasfield(typeof(phase), :measurements_period) ? phase.measurements_period : 1
    return hasfield(typeof(phase), :measurements) && due(period, k)
end

"""
    StopRun()

the exception that unwinds a run asked to stop, once its checkpoint is written
"""
struct StopRun <: Exception end

"""
    commit!(sim, phase, step; last = false)

record the simulation as the point to resume from, at step `step` of the phase `phase`, the
last one when `last`, see `Commit`
"""
function commit!(sim::Simulation, phase::Int, step::Int; last::Bool = false)
    run = getfield(sim, :run)
    run.checkpointer.last = Commit(phase, step, last, run.phase_time, sim.time, sim.state,
                                   committed_values(getfield(sim, :values)),
                                   resume_marks(run.outputs))
    return nothing
end

"""
    checkpoint_or_stop(run)

write the last commit when a checkpoint is due or a stop requested, and stop then, see `StopRun`
"""
function checkpoint_or_stop(run::RunData)
    c = run.checkpointer
    stop = stop_requested(c)
    if stop || checkpoint_due(c)
        write_checkpoint(c)
    end
    if stop
        throw(StopRun())
    end
    return nothing
end

"""
    step_line(sim, phase, k)

the line the log takes at the end of step `k` of `phase`: its number, the time, and the values
the phase provides, but for states
"""
function step_line(sim::Simulation, phase, k::Int)
    values = getfield(sim, :values)
    return join([ "step $k, time $(sim.time)" ;
                  [ "$name $(values[name])" for name in provided_values(phase)
                    if haskey(values, name) && !(values[name] isa AbstractState) ] ], ", ")
end

"""
    TensorMixedStates.step_done(sim, phase, k; last = false)

the end of step `k` of `phase`, the simulation `sim` holding what it reached: its measurements,
when due, and its line in the log, then the commit of the step, a checkpoint when one is due,
and the stop of the run when one is asked for. It returns the simulation, its value `step` set
to `k`. A phase of a single call calls it wherever it can stop, a sweep of a search for
instance: its `run_phase(phase, sim, done)` is then called with the steps already done. With
`last`, the phase ends there, `sim` holding its result, which a resumed run takes without
running the phase again: a search ending on its tolerance would otherwise go on.
"""
function step_done(sim::Simulation, phase, k::Int; last::Bool = false)
    sim = with_values(sim; step = k)
    if step_due(phase, k)
        output(sim, phase.measurements)
    end
    log_message(sim, step_line(sim, phase, k))
    run = getfield(sim, :run)
    commit!(sim, run.phase, k; last)
    checkpoint_or_stop(run)
    return sim
end

"""
    TensorMixedStates.run_phase(phase, sim, done)

run `phase` on the simulation `sim`, its first `done` steps already done, and return the
simulation it leaves. For a phase made of steps, it prepares the phase and runs the steps after
`done`, ending each with `step_done`; a phase of a single call has a method of its own.
"""
function run_phase(phase::AbstractPhase, sim::Simulation, done::Int)
    if !hasmethod(prepare_phase, Tuple{typeof(phase), Simulation})
        error("$(typeof(phase)) has neither a method of TensorMixedStates.run_phase nor one of " *
              "TensorMixedStates.prepare_phase, so runTMS does not know how to run it")
    end
    # called for the phase and the state as they are, see output
    p = invokelatest(prepare_phase, phase, Simulation(sim, sim.state, getfield(sim, :run).phase_time))
    for k in done + 1:invokelatest(step_count, p)
        sim = step_done(invokelatest(run_step, p, sim, k), phase, k)
    end
    return sim
end
