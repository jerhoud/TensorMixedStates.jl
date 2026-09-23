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

run_phase(sim::Simulation, sd::SimData) =
    log_phase(sim, sd.phases)

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
    else
        if phase.state isa State
            # `type` is what the phase was asked for, so a State given in the other
            # representation is converted rather than silently kept as it is
            state = as_representation(sim, R, phase.state)
        elseif isnothing(phase.system)
            error("CreateState needs a system or a State object")
        else
            state = State{R}(phase.system, phase.state)
        end
        if phase.randomize ≠ 0
            state = RandomState(state, phase.randomize)
        end
    end
    return Simulation(sim, state)
end


function run_phase(sim::Simulation, phase::ToMixed)
    if sim.state isa State{Mixed}
        log_msg(sim, "State is already in mixed representation")
    else
        log_msg(sim, "Creating mixed representation with $(length(sim)) sites")
        sim = truncate(mix(sim); phase.limits)
        log_msg(sim, "State is now in mixed representation")
    end
    return sim
end


function run_phase(sim::Simulation, phase::Evolve)
    nsweeps = Int(round(phase.duration / phase.time_step))
    duration = phase.time_step * nsweeps
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
    algo = phase.algo
    if algo isa ApproxW
        state = approx_W(pre, duration, state;
            coefs, algo.n_hermitianize, nsweeps, algo.order, algo.w, time_start = sim.time, phase.limits,
            observer! = ApproxWObserver(sim, phase.measures, phase.measures_period),
            first_sweep = first_sweep!(sim.checkpoint))
    else
        state = tdvp(pre, duration, state;
            coefs, algo.n_expand, algo.n_hermitianize, nsweeps, time_start = sim.time, phase.limits,
            observer! = TdvpObserver(sim, phase.measures, phase.measures_period),
            first_sweep = first_sweep!(sim.checkpoint))
    end
    # a phase cut short by a checkpoint stops at the time it actually reached
    return Simulation(sim, state, sim.checkpoint.stopping ? sim.checkpoint.simtime : time_stop)
end


function run_phase(sim::Simulation, phase::Gates)
  log_msg(sim, "Applying $(length(prodsubs(phase.gates))) gates")
  return apply(phase.gates, sim; phase.limits)
end


function run_phase(sim::Simulation, phase::GroundState)
    done = first_sweep!(sim.checkpoint) - 1
    log_msg(sim, "Optimizing state with $(phase.nsweeps - done) sweeps of Dmrg")
    if done ≥ phase.nsweeps
        return sim
    end
    e, sim = dmrg(phase.hamiltonian, sim; phase.nsweeps, first_sweep = done + 1,
        phase.limits, phase.noise,
        observer! = DmrgObserver(sim, phase.measures, phase.measures_period, phase.tolerance, done))
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
    done = first_sweep!(sim.checkpoint) - 1
    log_msg(sim, "Searching for steady state with $(phase.nsweeps - done) sweeps of Dmrg")
    if done ≥ phase.nsweeps
        return sim
    end
    e, sim = steady_state(phase.lindbladian, sim;
        phase.nsweeps, first_sweep = done + 1, phase.limits, phase.mpo_limits, phase.mpo_algo,
        observer! = DmrgObserver(sim, phase.measures, phase.measures_period, phase.tolerance, done))
    log_msg(sim, "Done, dmrg final value is $e (0 for steady state)")
    return sim
end