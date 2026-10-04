# run_phase, which runs each type of phase on a simulation and returns the simulation it leaves
# behind; a phase type of one's own gets a method of it, with resume_step and run_steps to
# resume it where it stopped.

export resume_step, run_steps

"""
    resume_step(sim)

the steps, or sweeps, the phase being run has already done, and the energy dmrg had reached at
the last of them, or the logarithm of the trace `thermal_state` had reached, as
`(done, energy)`: `(0, nothing)` unless it is the phase a resumed run starts from. The resume point is then consumed, so that it is read once.

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
    return (r.sweep, r.energy)
end

"""
    run_steps(f, sim, nsteps)

run the steps `1:nsteps` of a phase of your own, `f(sim, k)` doing step `k` and returning the
simulation it leaves behind, its measurements written, with `output` for instance, which takes
`sweep = k` for the `:sweep` measurement. Between two steps, a checkpoint is written when one is
due and the phase stops when the simulation is asked to, and a resumed run continues after the
last step done, from the state and the simulation time it had reached. Outside `runTMS` the
steps simply run one after the other.

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
            output(sim, p.measures; sweep = k)
            return sim
        end
"""
function run_steps(f, sim::Simulation, nsteps::Int)
    # read before resume_step consumes it
    r = sim.checkpoint.resume
    done, _ = resume_step(sim)
    if done > 0
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
        sim = f(sim, k)
        c.sweeps = true
        if !(sim isa Simulation)
            error("step $k of run_steps returned a $(typeof(sim)), where it has to return " *
                  "the simulation it leaves behind")
        end
        if c.stopping
            break
        end
        if sweep_commit!(sim, sim.state, sim.time, k)
            break
        end
    end
    return sim
end

"""
    run_phase(sim::Simulation, phase)

run one phase on a simulation and return the simulation it leaves behind.

This is where a phase of your own plugs in: define a struct with the three fields every phase
has, `name`, `time_start` and `final_measures`, and a method of `TensorMixedStates.run_phase`
for it. The full name is needed: `run_phase` is not exported, and a function of your own
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
say so, rather than fail with a bare `MethodError` inside a run.
"""
run_phase(sim::Simulation, phase) =
    error("there is no run_phase method for $(typeof(phase)), so runTMS does not know " *
          "how to run it. It has the fields of a phase, so what is missing is the method " *
          "itself: define TensorMixedStates.run_phase(::Simulation, ::$(typeof(phase))), " *
          "returning the simulation the phase leaves behind")

"""
    as_representation(sim, R, ::State)

the given state in representation `R`, for a `CreateState` handed a `State` rather than a
description. A pure state is mixed if `R` is `Mixed`; a mixed state is refused if `R` is
`Pure`, since it holds no purification to go back to.
"""
as_representation(::Simulation, ::Type{R}, state::State{R}) where R = state
function as_representation(sim::Simulation, ::Type{Mixed}, state::State{Pure})
    log_msg(sim, "Creating mixed representation with $(length(state)) sites")
    return mix(state)
end
as_representation(::Simulation, ::Type{Pure}, ::State{Mixed}) =
    error("CreateState was asked for a pure state but given a mixed one, which cannot be " *
          "turned back into a pure state")

"""
    run_search(solve, sim, phase, what, final_line)

run a phase that searches a state by dmrg, `GroundState` or `SteadyState`, `solve(sim;
options...)` calling its solver: from the sweep a resumed run had reached, with a
`DmrgObserver`, the log saying `what` is being done and, unless the run stops for a
checkpoint, the line `final_line(e)` with the value reached. A search whose checkpoint fell on
its last sweep has only that line left to write, with the value the checkpoint recorded.
"""
function run_search(solve, sim::Simulation, phase, what::String, final_line)
    done, e = resume_step(sim)
    if done < phase.nsweeps
        log_msg(sim, "$what with $(phase.nsweeps - done) sweeps of Dmrg")
        e, sim = solve(sim; phase.nsweeps, first_sweep = done + 1, phase.limits,
            observer! = DmrgObserver(sim, phase.measures, phase.measures_period, phase.tolerance,
                                     done; phase.nsweeps, energy = e))
    end
    # a search stopped for a checkpoint is not done, and its resume writes the line
    if !sim.checkpoint.stopping
        log_msg(sim, final_line(e))
    end
    return sim
end

"""
    evolve(algo, state, sim, phase; evolver, coefs, nsteps)

the simulation `sim` once the `Evolve` phase `phase` has evolved its state `state` with the
algorithm `algo`, in `nsteps` steps covering `phase.duration`: `evolver` is the evolver of the
phase, and `coefs` the functions of time of a time dependent one, or `nothing`. The state is
given apart from the simulation so that a method is chosen by its type as well as by that of
the algorithm.

An algorithm of an extension, or an algorithm for the state of an extension, comes with a
method of its own. It is called before the phase has read its resume point, so that a method
can read it with `resume_step`, or run its steps with `run_steps`, which resumes, stops and
checkpoints them; `output(sim, phase.measures; sweep)` writes the measurements of the phase,
every `phase.measures_period` steps.
"""
function evolve(algo::Tdvp, state::State, sim::Simulation, phase::Evolve; evolver, coefs,
                nsteps)
    done, _ = resume_step(sim)
    st = tdvp(evolver, phase.duration, state; coefs, algo.n_hermitianize, nsweeps = nsteps,
              time_start = sim.time, phase.limits, first_sweep = done + 1, algo.n_expand,
              algo.krylov, observer! = TdvpObserver(sim, phase.measures, phase.measures_period))
    return Simulation(sim, st)
end

function evolve(algo::ApproxW, state::State, sim::Simulation, phase::Evolve; evolver, coefs,
                nsteps)
    done, _ = resume_step(sim)
    st = approx_W(evolver, phase.duration, state; coefs, algo.n_hermitianize, nsweeps = nsteps,
                  time_start = sim.time, phase.limits, first_sweep = done + 1, algo.order,
                  algo.w, algo.apply_algo,
                  observer! = ApproxWObserver(sim, phase.measures, phase.measures_period))
    return Simulation(sim, st)
end

evolve(algo, state, ::Simulation, ::Evolve; kwargs...) =
    error("TensorMixedStates.evolve has no method for $(typeof(algo)) on a $(typeof(state))")

function run_phase(sim::Simulation, phase::CreateState{R}) where {R <: PM}
    if !isnothing(phase.seed)
        Random.seed!(phase.seed)
    end
    if isnothing(phase.state)
        if phase.randomize == 0
            error("CreateState without state nor randomize: no state created !")
        elseif isnothing(phase.system)
            error("CreateState needs a system to create a random state")
        else
            state = RandomState{R}(phase.system, phase.randomize)
        end
    elseif phase.state isa State
        # a random mixed state is drawn from a purification, which a State does not give
        if R === Mixed && phase.randomize ≠ 0
            error("CreateState cannot randomize a State into a mixed state: give a description " *
                  "of the state, whose purification it is drawn from")
        end
        # `type` is what the phase was asked for, so a State given in the other
        # representation is converted rather than silently kept as it is
        state = as_representation(sim, R, phase.state)
        if phase.randomize ≠ 0
            state = RandomState(state, phase.randomize)
        end
    elseif isnothing(phase.system)
        error("CreateState needs a system or a State object")
    elseif phase.randomize == 0
        state = State{R}(phase.system, phase.state)
    elseif R === Mixed
        # a mixed state is drawn from the states its purification starts from, the only way
        # there is on a system that conserves something
        state = RandomState{Mixed}(phase.system, phase.state, phase.randomize)
    else
        state = RandomState(State{Pure}(phase.system, phase.state), phase.randomize)
    end
    return Simulation(sim, state)
end

function run_phase(sim::Simulation, phase::ToMixed)
    if sim.state isa State{Mixed}
        log_msg(sim, "State is already in mixed representation")
        sim = truncate(sim; phase.limits)
    else
        log_msg(sim, "Creating mixed representation with $(length(sim)) sites")
        sim = truncate(mix(sim); phase.limits)
        log_msg(sim, "State is now in mixed representation")
    end
    return sim
end

function run_phase(sim::Simulation, phase::Evolve)
    # the duration gives the direction, and a step of the other sign is adjusted as one that
    # does not divide it: of the opposite sign, it made no step while the time went on
    nsteps = round(Int, abs(phase.duration / phase.time_step))
    if nsteps == 0
        log_msg(sim, "Skipping an evolution of $(phase.duration), shorter than half a time step")
        return sim
    end
    # the step is adjusted rather than the duration, so that the phase ends where it was asked
    # to: a duration of 1 in steps of 0.3 stopped at 0.9
    duration = phase.duration
    if !(duration / nsteps ≈ phase.time_step)
        log_msg(sim, "Taking a time step of $(duration / nsteps) rather than $(phase.time_step) " *
                     "to cover the duration $duration in $nsteps steps")
    end
    time_stop = sim.time + duration
    # a phase resumed in its course goes on from the time of its checkpoint, not from its start
    r = sim.checkpoint.resume
    log_msg(sim, "Evolving state from simulation time $(isnothing(r) ? sim.time : r.time) to $(time_stop)")
    evolver, coefs = phase.evolver isa Pair ? (first(phase.evolver), last(phase.evolver)) :
                                              (phase.evolver, nothing)
    sim = evolve(phase.algo, sim.state, sim, phase; evolver, coefs, nsteps)
    # a phase cut short by a checkpoint stops at the time it actually reached
    c = sim.checkpoint
    return Simulation(sim, sim.state, c.stopping ? c.last.time : time_stop)
end

function run_phase(sim::Simulation, phase::Gates)
    log_msg(sim, "Applying $(length(prodsubs(phase.gates))) gates")
    return apply(phase.gates, sim; phase.limits)
end

run_phase(sim::Simulation, phase::GroundState) =
    run_search((sim; kwargs...) -> dmrg(phase.hamiltonian, sim; phase.noise, phase.krylov,
                                        kwargs...),
               sim, phase, "Optimizing state", e -> "Done, dmrg final energy is $e")

function run_phase(sim::Simulation, phase::SaveState)
    # a file of the simulation, as a checkpoint, would be overwritten by it, or overwrite it
    check_destination(sim, phase.file)
    save_state(phase.file, phase.statename, sim.state)
    return sim
end

function run_phase(sim::Simulation, phase::LoadState)
    st = load_state(phase.file, phase.statename)
    # Limits() truncates nothing, which a representation of one's own need not support
    return Simulation(sim, phase.limits == Limits() ? st : truncate(st; phase.limits))
end

function run_phase(sim::Simulation, phase::PartialTrace)
    pos = phase.trace_positions
    keep = phase.keep_positions
    if isnothing(pos) == isnothing(keep)
        error("PartialTrace requires one and only one of trace_positions and keep_positions")
    end
    return isnothing(pos) ? partial_trace(sim, keep; keepers = true) : partial_trace(sim, pos)
end

function run_phase(sim::Simulation, phase::Weaken)
    from = symmetries(sim.state.system)
    sim = isnothing(phase.target) ? weaken(sim) : weaken(sim, phase.target)
    log_msg(sim, "Symmetries weakened from $from to $(symmetries(sim.state.system))")
    return sim
end

function run_phase(sim::Simulation, phase::SteadyState)
    if sim.state isa State{Pure}
        error("state must be in mixed representation for computing steady state")
    end
    return run_search(
        (sim; kwargs...) -> steady_state(phase.lindbladian, sim; phase.mpo_limits,
                                         phase.mpo_algo, phase.noise, phase.krylov, kwargs...),
        sim, phase, "Searching for steady state",
        e -> "Done, dmrg final value is $e (0 for steady state)")
end

function run_phase(sim::Simulation, phase::Thermalize)
    # adjusted as the time step of Evolve, the step taking the sign of beta
    nsteps = round(Int, abs(phase.beta / phase.beta_step))
    if nsteps == 0
        log_msg(sim, "Skipping a thermalization to beta $(phase.beta), shorter than half a step")
        return sim
    end
    beta = phase.beta
    if !(beta / nsteps ≈ phase.beta_step)
        log_msg(sim, "Taking a step of $(beta / nsteps) rather than $(phase.beta_step) to reach " *
                     "beta $beta in $nsteps steps")
    end
    done, log_trace = resume_step(sim)
    log_msg(sim, "Thermalizing state to beta $beta")
    algo = phase.algo
    l, sim = thermal_state(phase.hamiltonian, beta, sim; nsteps, first_step = done + 1,
                           log_trace = something(log_trace, 0.), algo.n_expand,
                           algo.n_hermitianize, phase.limits, algo.krylov,
                           observer! = ThermalObserver(sim, phase.measures, phase.measures_period))
    log_msg(sim, "Done, log_trace is $l")
    return sim
end
