run_phase(sim::Simulation, sd::SimData) =
    log_phase(sim, sd.phases)

function run_phase(sim::Simulation, phase::CreateState)
    if isnothing(phase.state)
        if phase.randomize == 0
            error("CreateState without state nor randomize: no state created !")
        else
            state = random_state(phase.type, phase.system, phase.randomize)
        end
    else
        state = State(phase.type, phase.system, phase.state)
        if phase.randomize ≠ 0
            state = random_state(state, phase.randomize)
        end
    end
    return Simulation(sim, state)
end


function run_phase(sim::Simulation, phase::ToMixed)
    if sim.state.type == Mixed 
        log_msg(sim, "State is already in mixed representation")
    else
        log_msg(sim, "Creating mixed representation with $(length(sim)) sites")
        sim = truncate(mix(sim); phase.limits.cutoff, phase.limits.maxdim)
        log_msg(sim, "State is now in mixed representation")
    end
    return sim
end


function run_phase(sim::Simulation, phase::Evolve)
    nsweeps = Int(round(phase.duration / phase.time_step))
    duration = phase.time_step * nsweeps
    time_stop = sim.time + duration
    log_msg(sim, "Evolving state from simulation time $(sim.time) to $(time_stop)")
    time_dep = isa(phase.evolver, Pair)
    if time_dep
        evolver = first(phase.evolver)
        coefs = last(phase.evolver)
    else
        evolver = phase.evolver
        coefs = nothing
    end
    state = sim.state
    tp = state.type
    if multipleLit(evolver)
        error("Evolver cannot contain multiple site operators")
    end
    if tp == Pure && dissipLit(evolver)
        error("Evolving error: state must be in mixed representation to use dissipators")
    end
    pre = PreMPO(state, evolver)
    algo = phase.algo
    if algo isa ApproxW
        state = approx_W(pre, duration, state;
            coefs, n_symmetrize = phase.symmetrize, nsweeps, algo.order, algo.w, time_start = sim.time,
            phase.limits.cutoff, phase.limits.maxdim,
            observer! = ApproxWObserver(sim, phase.measures, phase.measures_period))
    else
        state = tdvp(pre, duration, state;
            coefs, n_expand = phase.exapnd, n_symmetrize = phase.symmetrize, nsweeps, time_start = sim.time,
            phase.limits.cutoff, phase.limits.maxdim,
            observer! = TdvpObserver(sim, phase.measures, phase.measures_period))
    end
    return Simulation(sim, state, time_stop)
end


function run_phase(sim::Simulation, phase::Gates)
  log_msg(sim, "Applying $(length(phase.gates.ls)) gates")
  return apply(phase.gates, sim; phase.limits.cutoff, phase.limits.maxdim)
end


function run_phase(sim::Simulation, phase::Dmrg)
    log_msg(sim, "Optimizing state with $(phase.nsweeps) sweeps of Dmrg")

    if multipleLit(phase.hamiltonian) || dissipLit(phase.hamiltonian)
        error("Dmrg hamiltonian cannot contain multiple site operators nor dissipators")
    end
    e, sim = dmrg(phase.hamiltonian, sim; phase.nsweeps, phase.maxdim, phase.cutoff,
        observer! = DmrgObserver(sim, phase.measures, phase.measures_period, phase.tolerance))
    log_msg(sim, "Done, dmrg final energy is $e")
    return sim
end

run_phase(sim::Simulation, phase::SaveState) =
    save_state(phase.file, phase.statename, sim.state)
    
run_phase(sim::Simulation, phase::LoadState) =
    Simulation(sim, truncate(load_state(phase.file, phase.statename); phase.limits.cutoff, phase.limits.maxdim))
