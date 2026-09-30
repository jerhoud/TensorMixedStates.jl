# run_phase, which runs each type of phase on a simulation and returns the simulation it leaves
# behind; a phase type of one's own gets a method of it.

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

A phase that does not drive a solver has no point to be resumed from inside: the `stop` file
and `max_time` are only seen once it ends, an interrupt stops it at once, and a resumed run
runs it again from its start.

A phase of your own that drives a solver with `TdvpObserver`, `ApproxWObserver` or
`DmrgObserver` is stopped and checkpointed as those of the library are, but a checkpoint
written while it runs resumes it from its start. To resume it at the sweep it had reached,
read `done, energy = TensorMixedStates.resume_sweeps!(sim.checkpoint)` before starting the
solver, start it at `first_sweep = done + 1`, and hand `done` and `energy` to a
`DmrgObserver`.

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
    done, e = resume_sweeps!(sim.checkpoint)
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

function run_phase(sim::Simulation, phase::CreateState{R}) where R
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
    nsweeps = round(Int, phase.duration / phase.time_step)
    if nsweeps == 0
        log_msg(sim, "Skipping an evolution of $(phase.duration), shorter than half a time step")
        return sim
    end
    # the step is adjusted rather than the duration, so that the phase ends where it was asked
    # to: a duration of 1 in steps of 0.3 stopped at 0.9
    duration = phase.duration
    if !(duration / nsweeps ≈ phase.time_step)
        log_msg(sim, "Taking a time step of $(duration / nsweeps) rather than $(phase.time_step) " *
                     "to cover the duration $duration in $nsweeps steps")
    end
    time_stop = sim.time + duration
    log_msg(sim, "Evolving state from simulation time $(sim.time) to $(time_stop)")
    evolver, coefs = phase.evolver isa Pair ? (first(phase.evolver), last(phase.evolver)) :
                                              (phase.evolver, nothing)
    state = sim.state
    # PreMPO adapts the evolver to the representation of the state, and handles the vector
    # form of a time dependent evolver
    pre = PreMPO(state, evolver)
    done, _ = resume_sweeps!(sim.checkpoint)
    algo = phase.algo
    common = (; coefs, algo.n_hermitianize, nsweeps, time_start = sim.time, phase.limits,
              first_sweep = done + 1)
    if algo isa ApproxW
        state = approx_W(pre, duration, state; common..., algo.order, algo.w,
            observer! = ApproxWObserver(sim, phase.measures, phase.measures_period))
    else
        state = tdvp(pre, duration, state; common..., algo.n_expand,
            observer! = TdvpObserver(sim, phase.measures, phase.measures_period))
    end
    # a phase cut short by a checkpoint stops at the time it actually reached
    c = sim.checkpoint
    return Simulation(sim, state, c.stopping ? c.last.time : time_stop)
end

function run_phase(sim::Simulation, phase::Gates)
    log_msg(sim, "Applying $(length(prodsubs(phase.gates))) gates")
    return apply(phase.gates, sim; phase.limits)
end

run_phase(sim::Simulation, phase::GroundState) =
    run_search((sim; kwargs...) -> dmrg(phase.hamiltonian, sim; phase.noise, kwargs...),
               sim, phase, "Optimizing state", e -> "Done, dmrg final energy is $e")

function run_phase(sim::Simulation, phase::SaveState)
    save_state(phase.file, phase.statename, sim.state)
    return sim
end

run_phase(sim::Simulation, phase::LoadState) =
    Simulation(sim, truncate(load_state(phase.file, phase.statename); phase.limits))

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
                                         phase.mpo_algo, kwargs...),
        sim, phase, "Searching for steady state",
        e -> "Done, dmrg final value is $e (0 for steady state)")
end