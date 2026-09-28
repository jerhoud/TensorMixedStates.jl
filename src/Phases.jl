"""
    run_phase(sim::Simulation, phase)

run one phase on the given simulation and return the simulation it leaves behind.

This is where a phase of your own plugs in. Define a struct carrying the three fields the
machinery around a phase reads — `name`, `time_start` and `final_measures` — and a method of
`TensorMixedStates.run_phase` for it. `runTMS` then logs it, applies its `time_start`, calls
your method and takes its final measurements, exactly as for a phase of the library. Note
the full name: `run_phase` is not exported, so it has to be written out to add a method to
it rather than shadowed by one of your own.

The fallback method below exists so that an object that is not a phase says so, instead of
surfacing as a bare `MethodError` from somewhere inside a run.
"""
run_phase(sim::Simulation, phase) =
    error("there is no run_phase method for $(typeof(phase)), so runTMS does not know " *
          "how to run it. It has the fields of a phase, so what is missing is the method " *
          "itself: define TensorMixedStates.run_phase(::Simulation, ::$(typeof(phase))), " *
          "returning the simulation the phase leaves behind")

"""
    as_representation(sim, R, ::State)

the given state in representation `R`, for a `CreateState` handed a `State` object rather
than a description of one. A pure state is mixed on the way in, which is what `type` asked
for; the other direction does not exist, a mixed state holds no purification to go back to.
"""
as_representation(::Simulation, ::Type{R}, state::State{R}) where R = state
as_representation(sim::Simulation, ::Type{Mixed}, state::State{Pure}) = begin
    log_msg(sim, "Creating mixed representation with $(length(state)) sites")
    mix(state)
end
as_representation(::Simulation, ::Type{Pure}, ::State{Mixed}) =
    error("CreateState was asked for a pure state but given a mixed one, which cannot be " *
          "turned back into a pure state")

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
    nsweeps = Int(round(phase.duration / phase.time_step))
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
    time_dep = phase.evolver isa Pair
    if time_dep
        evolver = first(phase.evolver)
        coefs = last(phase.evolver)
    else
        evolver = phase.evolver
        coefs = nothing
    end
    state = sim.state
    # PreMPO adapts the evolver to the representation of the state, and handles the vector
    # form of a time dependent evolver
    pre = PreMPO(state, evolver)
    done, _ = resume_sweeps!(sim.checkpoint)
    algo = phase.algo
    if algo isa ApproxW
        state = approx_W(pre, duration, state;
            coefs, algo.n_hermitianize, nsweeps, algo.order, algo.w, time_start = sim.time, phase.limits,
            observer! = ApproxWObserver(sim, phase.measures, phase.measures_period),
            first_sweep = done + 1)
    else
        state = tdvp(pre, duration, state;
            coefs, algo.n_expand, algo.n_hermitianize, nsweeps, time_start = sim.time, phase.limits,
            observer! = TdvpObserver(sim, phase.measures, phase.measures_period),
            first_sweep = done + 1)
    end
    # a phase cut short by a checkpoint stops at the time it actually reached
    c = sim.checkpoint
    return Simulation(sim, state, c.stopping ? c.last.time : time_stop)
end


function run_phase(sim::Simulation, phase::Gates)
  log_msg(sim, "Applying $(length(prodsubs(phase.gates))) gates")
  return apply(phase.gates, sim; phase.limits)
end


function run_phase(sim::Simulation, phase::GroundState)
    done, e = resume_sweeps!(sim.checkpoint)
    # a search whose checkpoint fell on its last sweep has only its last line left to write,
    # with the energy the checkpoint recorded
    if done < phase.nsweeps
        log_msg(sim, "Optimizing state with $(phase.nsweeps - done) sweeps of Dmrg")
        e, sim = dmrg(phase.hamiltonian, sim; phase.nsweeps, first_sweep = done + 1,
            phase.limits, phase.noise,
            observer! = DmrgObserver(sim, phase.measures, phase.measures_period, phase.tolerance,
                                     done; phase.nsweeps, energy = e))
    end
    log_msg(sim, "Done, dmrg final energy is $e")
    return sim
end

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
    if isnothing(pos)
        return partial_trace(sim, keep; keepers = true)
    else
        return partial_trace(sim, pos)
    end
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
    done, e = resume_sweeps!(sim.checkpoint)
    # as for GroundState, a search whose checkpoint fell on its last sweep only writes its
    # last line
    if done < phase.nsweeps
        log_msg(sim, "Searching for steady state with $(phase.nsweeps - done) sweeps of Dmrg")
        e, sim = steady_state(phase.lindbladian, sim;
            phase.nsweeps, first_sweep = done + 1, phase.limits, phase.mpo_limits, phase.mpo_algo,
            observer! = DmrgObserver(sim, phase.measures, phase.measures_period, phase.tolerance,
                                     done; phase.nsweeps, energy = e))
    end
    log_msg(sim, "Done, dmrg final value is $e (0 for steady state)")
    return sim
end