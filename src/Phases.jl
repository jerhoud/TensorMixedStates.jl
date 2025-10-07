run_phase(sim::Simulation, sd::SimData) =
    log_phase(sim, sd.phases)

function run_phase(sim::Simulation, phase::CreateState{R}) where R
    if !isnothing(phase.seed)
        Random.seed!(phase.seed)
    end
    if isnothing(phase.state)
        if phase.randomize == 0
            error("CreateState without state nor randomize: no state created !")
        else
            state = RandomState{R}(phase.system, phase.randomize)
        end
    else
        if phase.state isa State
            state = phase.state
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
    if state isa State{Pure} && evolver isa IndexedOp{Mixed}
        error("Evolving error: state must be in mixed representation to use this evolver")
    elseif state isa State{Mixed} && evolver isa IndexedOp{Pure}
        evolver = Evolver(evolver)
    end
    pre = PreMPO(state, evolver)
    algo = phase.algo
    if algo isa ApproxW
        state = approx_W(pre, duration, state;
            coefs, algo.n_hermitianize, nsweeps, algo.order, algo.w, time_start = sim.time, phase.limits,
            observer! = ApproxWObserver(sim, phase.measures, phase.measures_period))
    else
        state = tdvp(pre, duration, state;
            coefs, algo.n_expand, algo.n_hermitianize, nsweeps, time_start = sim.time, phase.limits,
            observer! = TdvpObserver(sim, phase.measures, phase.measures_period))
    end
    return Simulation(sim, state, time_stop)
end


function run_phase(sim::Simulation, phase::Gates)
  log_msg(sim, "Applying $(length(prodsubs(phase.gates))) gates")
  return apply(phase.gates, sim; phase.limits)
end


function run_phase(sim::Simulation, phase::GroundState)
    log_msg(sim, "Optimizing state with $(phase.nsweeps) sweeps of Dmrg")
    e, sim = dmrg(phase.hamiltonian, sim; phase.nsweeps, phase.limits,
        observer! = DmrgObserver(sim, phase.measures, phase.measures_period, phase.tolerance))
    log_msg(sim, "Done, dmrg final energy is $e")
    return sim
end

run_phase(sim::Simulation, phase::SaveState) =
    save_state(phase.file, phase.statename, sim.state)
    
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

function run_phase(sim::Simulation, phase::SteadyState)
    log_msg(sim, "Searching for steady state with $(phase.nsweeps) sweeps of Dmrg")
    if sim.state isa State{Pure}
        error("state must be in mixed representation for computing steady state")
    end
    l = make_mpo(sim.state, phase.lindbladian)
    l2 = apply(dag(l), l; phase.mpo_limits.cutoff, phase.mpo_limits.maxdim)
    e, sim = dmrg(l2, sim; phase.nsweeps, phase.limits,
        observer! = DmrgObserver(sim, phase.measures, phase.measures_period, phase.tolerance))
    log_msg(sim, "Done, dmrg final value is $e (0 for steady state)")
    return sim
end