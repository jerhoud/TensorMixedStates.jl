# The interface of the phases of a simulation, those of the library (src/phases) and those of
# one's own: AbstractPhase, run_phase, the resumption of a phase where it stopped, and the
# phases creating the state.

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
    TensorMixedStates.creates_state(phase)

whether `phase` gives the simulation a state of its own, as `CreateState` and `LoadState` do,
false by default. A simulation starts with such a phase. A phase of one's own creating the
state gets a method of it: `TensorMixedStates.creates_state(::MyPhase) = true`.
"""
creates_state(_) = false

"""
    TensorMixedStates.phase_system(phase)

the system on which a phase creating the state creates it, when the phase knows it
beforehand, `nothing` otherwise (the default). A resumed run puts the state of its checkpoint
back on that system, so that it can still be compared with states built on it. A phase of
one's own creating the state on a system it is given gets a method of it.
"""
phase_system(_) = nothing

"""
    resume_step(sim)

the steps, or sweeps, the phase being run has already done, and the value carried from the last
of them (the energy of `dmrg`, the logarithm of the trace of `thermal_state`, the value of
`run_steps`), as `(done, carried)`; `(0, nothing)` unless it is the phase a resumed run starts
from. The resume point is consumed: it is read once.

Only a phase that calls it has its steps committed: a resume hands it the state of its last
committed step, which only a phase continuing after the steps done can go on from. A phase that
never calls it is resumed from its start.

A phase of one's own driving a solver starts it at `first_sweep = done + 1`, hands `done`, the
energy carried and its `nsweeps` to a `DmrgObserver`, and skips the solver when `done` has
reached `nsweeps`, as after a search stopped by its tolerance or checkpointed on its last
sweep: run for no sweep, `dmrg` would give the energy of the state rather than the one carried.
A loop of one's own is best written with `run_steps`, which calls this itself and must not have
it called again within its steps.
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

the simulation time the phase being run goes on from: that of its last committed step when a
resumed run starts in it, the time of the simulation otherwise. Unlike `resume_step`, it does
not consume the resume point.
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
reached, see `stopped`. A phase whose solver stopped early gives the simulation this time,
where the solver would count the whole duration.
"""
function committed_time(sim::Simulation)
    k = sim.checkpoint.last
    return isnothing(k) ? sim.time : k.time
end

"""
    run_steps(f, sim, nsteps)
    run_steps(f, sim, nsteps; carry)

run the steps `1:nsteps` of a phase of one's own, `f(sim, k)` doing step `k`, writing its
measurements (with `output(sim, measurements; sweep = k)` for instance) and returning the
simulation it leaves behind. Between two steps, a checkpoint is written when one is due and the
phase stops when asked to; a resumed run continues after the last step done, from the state and
the simulation time it had reached. Outside `runTMS` the steps simply run in turn.

Given `carry`, other than `nothing`, the steps carry a value, a sum or a count for instance:
`f(sim, k, carried)` returns `(sim, carried)`, the first step receives `carry`, and `run_steps`
returns `(sim, carried)`. The value is given back to a resumed run through the json of a
checkpoint, and has to come back as it was, of the same type, which is checked after each step
of a run with a directory: numbers, strings, and vectors, matrices and dictionaries with string
keys of them, typed by what they hold (a `Vector{Any}` of floats comes back a
`Vector{Float64}`). A random number generator is carried by its `UInt64` words,
`Xoshiro(words...)` giving it back.

A step may run a solver with an observer of the package, `TdvpObserver` for instance: its
sweeps are not committed, and a stop it honours ends the phase after the last step done, the
unfinished one being run again whole on a resume. A step must not call `resume_step`, which
`run_steps` has called.

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
        log_message(sim, "up after \$ups kicks")
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
        # a resume hands the state of the last committed step but the time the phase started
        # from, which the solvers count from
        sim = Simulation(sim, sim.state, r.time)
    end
    c = sim.checkpoint
    for k in done + 1:nsteps
        # nothing is committed while the step runs, or the sweeps of a solver with an observer
        # would be committed as steps of the phase
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
        # only a run with a directory writes checkpoints
        if !isnothing(carry) && !isempty(c.dir)
            check_carried(carried, k)
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

A phase of one's own is a subtype of `AbstractPhase` with the fields `name`, `time_start` and
`final_measurements`, and a method of `TensorMixedStates.run_phase`, named in full: `run_phase`
is not exported, and a function of one's own of that name would shadow it. `runTMS` logs the
phase, applies its `time_start`, calls the method and takes the final measurements. Within the
method, `output` measures the simulation, `log_message` writes to its log and `get_sim_file`
gives a file of the simulation to write anything else to.

A phase written with `run_steps` is stopped, checkpointed and resumed between two steps. One
driving a solver with `TdvpObserver`, `ApproxWObserver` or `DmrgObserver` is stopped and
checkpointed between two sweeps, and resumed at the sweep it had reached if it reads
`done, energy = resume_step(sim)` before starting the solver, see `resume_step`. Any other
phase is not resumed from inside: the `stop` file and `max_time` are only seen once it ends, an
interrupt stops it at once, and a resumed run runs it again from its start.

The method documented here is the fallback, saying that an object has no method. The phases of
the library, in `src/phases`, are as many examples.
"""
run_phase(sim::Simulation, phase) =
    error("there is no run_phase method for $(typeof(phase)), so runTMS does not know " *
          "how to run it. It has the fields of a phase, so what is missing is the method " *
          "itself: define TensorMixedStates.run_phase(::Simulation, ::$(typeof(phase))), " *
          "returning the simulation the phase leaves behind")
