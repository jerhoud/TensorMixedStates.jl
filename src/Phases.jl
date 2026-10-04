# The interface of the phases of a simulation, through which those of the library, each in a file
# of src/phases, and those of one's own are written: their supertype AbstractPhase, run_phase,
# resume_step, resume_time, committed_time and run_steps to resume a phase where it stopped, and
# creates_state and phase_system for a phase creating the state.

export AbstractPhase, resume_step, resume_time, committed_time, run_steps

"""
    abstract type AbstractPhase

the supertype of the phases of a simulation, those of the library, `CreateState`, `LoadState`,
`SaveState`, `ToMixed`, `Evolve`, `Gates`, `GroundState`, `SteadyState`, `Thermalize`,
`PartialTrace` and `Weaken`, and those of one's own, see `TensorMixedStates.run_phase`. Every
phase has at least these three fields:

- `name`: the name of the phase, written in the log
- `time_start`: the simulation time the clock is set to when the phase starts (default
  `nothing`, keeping the current time)
- `final_measurements`: the measurements to make at the end of the phase, see `output` (default
  `[]`)
"""
abstract type AbstractPhase end

# A phase is printed field by field, read from its type rather than written out, so that a field
# added to a phase shows up in the log and in `prog.jl` without anything else to change, those
# of one's own as those of the library.
function show(io::IO, s::AbstractPhase)
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

"""
    TensorMixedStates.creates_state(phase)

whether `phase` gives the simulation a state of its own, as `CreateState` and `LoadState` do,
false by default. A simulation has to start with such a phase, having no state before it. A
phase of one's own creating the state, the adapter of another library for instance, gets a
method of it: `TensorMixedStates.creates_state(::MyPhase) = true`.
"""
creates_state(_) = false

"""
    TensorMixedStates.phase_system(phase)

the system on which a phase creating the state creates it, when the phase knows it
beforehand, and `nothing` otherwise, the default. A resumed run puts the state of its checkpoint
back on the system of the last phase creating the state before the point it resumes from, so
that the state can still be compared with states built on that system. `CreateState` gives its
`system`, and a phase of one's own creating the state on a system it is given gets a method of
it.
"""
phase_system(_) = nothing

"""
    resume_step(sim)

the steps, or sweeps, the phase being run has already done, and the value it carried from the
last of them, as `(done, carried)`: the energy dmrg had reached, the logarithm of the trace
`thermal_state` had reached, or the value `run_steps` carries, `(0, nothing)` unless it is the
phase a resumed run starts from. The resume point is then consumed, so that it is read once.

Calling it is also what lets the steps of the phase be committed. A resume hands the phase the
state of its last committed step, and only a phase that reads the steps done and continues
after them goes on correctly from there; any other would run all its steps again on that
state. So the steps of a phase that never calls this are not committed, and a checkpoint
written while it runs resumes it from its start.

A phase of your own driving a solver starts it at `first_sweep = done + 1`, and hands `done`,
`energy` and its `nsweeps` to a `DmrgObserver`, which records a search stopped by its tolerance
as having done them all. The phase has to skip the solver when `done` has reached `nsweeps`, as
a search stopped by its tolerance or checkpointed on its last sweep has: run again for no
sweep, `dmrg` would give the energy of the state given rather than `energy`. A loop of your own is best
written with `run_steps`, which calls this for it and must then not have it called again
within its steps.
"""
function resume_step(sim::Simulation)
    c = sim.checkpoint
    c.sweeps = true
    r = c.resume
    if isnothing(r) || isnothing(c.last) || r.phase ≠ c.last.phase
        return (0, nothing)
    end
    c.resume = nothing
    return (r.sweep, r.carried)
end

"""
    resume_time(sim)

the simulation time the phase being run goes on from when a resumed run starts in it, that of
the step it had last committed, and otherwise the time of the simulation, without consuming the
resume point, which `resume_step` reads: a phase that drives a solver logs with it where it
starts from.
"""
function resume_time(sim::Simulation)
    c = sim.checkpoint
    r = c.resume
    if isnothing(r) || isnothing(c.last) || r.phase ≠ c.last.phase
        return sim.time
    end
    return r.time
end

"""
    committed_time(sim)

the simulation time of the last step committed, which a phase stopped for a checkpoint has
reached, see `stopped`: a phase whose solver stopped before its end gives the simulation that
time, which the solver, given the whole duration, does not.
"""
function committed_time(sim::Simulation)
    k = sim.checkpoint.last
    return isnothing(k) ? sim.time : k.time
end

"""
    run_steps(f, sim, nsteps)
    run_steps(f, sim, nsteps; carry)

run the steps `1:nsteps` of a phase of your own, `f(sim, k)` doing step `k` and returning the
simulation it leaves behind, its measurements written, with `output` for instance, which takes
`sweep = k` for the `:sweep` measurement. Between two steps, a checkpoint is written when one is
due and the phase stops when the simulation is asked to, and a resumed run continues after the
last step done, from the state and the simulation time it had reached. Outside `runTMS` the
steps simply run one after the other.

Given `carry`, a value other than `nothing`, the steps carry a value from one to the next, a
sum or a count for instance: `f(sim, k, carried)` returns `(sim, carried)`, the first step
receives `carry`, and `run_steps` returns `(sim, carried)` after the last. The value is
committed with each step and given back to a resumed run, so it is one a checkpoint can write:
numbers, strings, and arrays and dictionaries of them.

A step may run a solver with an observer of the package, `TdvpObserver` for instance: its
sweeps are not steps of the phase and are not committed, and a stop it honours ends the phase
with the last step done, the unfinished one being run again whole on a resume. A step does not
call `resume_step`, which `run_steps` has called already: calling it again would have the
sweeps of the solver committed as steps of the phase.

# Examples

    TensorMixedStates.run_phase(sim::Simulation, p::Kicks) =
        run_steps(sim, p.nkicks) do sim, k
            sim = apply(exp(-0.3im * X)(1), sim)
            sim = Simulation(sim, sim.state, sim.time + 0.1)
            output(sim, p.measurements; sweep = k)
            return sim
        end

    # the number of kicks the first qubit was found up after, carried through a resume
    function TensorMixedStates.run_phase(sim::Simulation, p::CountedKicks)
        sim, ups = run_steps(sim, p.nkicks; carry = 0) do sim, k, ups
            sim = apply(exp(-0.3im * X)(1), sim)
            return sim, ups + (real(expect(sim.state, Z(1))) > 0)
        end
        log_msg(sim, "up after \$ups kicks")
        return sim
    end
"""
function run_steps(f, sim::Simulation, nsteps::Int; carry = nothing)
    # read before resume_step consumes it
    r = sim.checkpoint.resume
    done, carried = resume_step(sim)
    if done == 0
        carried = carry
    else
        # a resume hands the phase the state of its last committed step, but the time the phase
        # started from, which the solvers count their steps from: a loop of one's own goes on
        # from the time it had reached
        sim = Simulation(sim, sim.state, r.time)
    end
    c = sim.checkpoint
    for k in done + 1:nsteps
        # a solver run within a step with an observer of the package would commit its sweeps as
        # steps of the phase: nothing is committed while the step runs. A stop the solver
        # honours still writes the checkpoint of the last step done, and leaves this one
        # unfinished, so that it is not committed either and is run again whole on a resume
        c.sweeps = false
        out = isnothing(carry) ? f(sim, k) : f(sim, k, carried)
        c.sweeps = true
        if isnothing(carry)
            sim = out
        elseif out isa Tuple && length(out) == 2
            sim, carried = out
        else
            error("step $k of run_steps returned a $(typeof(out)), where with carry it has to " *
                  "return the simulation it leaves behind and the value it carries")
        end
        if !(sim isa Simulation)
            error("step $k of run_steps returned a $(typeof(sim)), where it has to return " *
                  "the simulation it leaves behind")
        end
        if c.stopping
            break
        end
        if sweep_commit!(sim, sim.state, sim.time, k; carried)
            break
        end
    end
    return isnothing(carry) ? sim : (sim, carried)
end

"""
    run_phase(sim::Simulation, phase)

run one phase on a simulation and return the simulation it leaves behind.

This is where a phase of your own plugs in: define a subtype of `AbstractPhase` with the three
fields every phase has, `name`, `time_start` and `final_measurements`, and a method of
`TensorMixedStates.run_phase` for it. The full name is needed: `run_phase` is not exported, and a function of your own
of that name would shadow it rather than extend it. `runTMS` then logs the phase, applies its
`time_start`, calls the method and takes the final measurements, as for a phase of the
library. Within the method, `output` measures the simulation, `log_msg` writes to its log and
`get_sim_file` gives a file of the simulation to write anything else to.

The fields of a phase of your own are part of the fingerprint by which a checkpoint tells its
simulation, see `TensorMixedStates.phases_id`: keep in them what describes the phase, not what
changes from one run to the next or as it runs, which would have the checkpoint refused as
another simulation's.

A phase written as a loop of steps with `run_steps` is stopped, checkpointed and resumed
between two steps. One that drives a solver with `TdvpObserver`, `ApproxWObserver` or
`DmrgObserver` is stopped and checkpointed between two sweeps as those of the library are,
and resumed at the sweep it had reached if it reads `done, energy = resume_step(sim)` before
starting the solver, see `resume_step`. Any other phase has no point to be resumed from
inside: the `stop` file and `max_time` are only seen once it ends, an interrupt stops it at
once, and a resumed run runs it again from its start.

The method documented here is the fallback: it makes an object with no method of its own
say so, rather than fail with a bare `MethodError` inside a run. The phases of the library are
written each in a file of `src/phases` through the same interface, and are as many examples.
"""
run_phase(sim::Simulation, phase) =
    error("there is no run_phase method for $(typeof(phase)), so runTMS does not know " *
          "how to run it. It has the fields of a phase, so what is missing is the method " *
          "itself: define TensorMixedStates.run_phase(::Simulation, ::$(typeof(phase))), " *
          "returning the simulation the phase leaves behind")
