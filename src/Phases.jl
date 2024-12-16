function run_phase(sim::Simulation, sd::SimData)
    for phase in sd.phases
        sim = log_phase(sim, phase)
    end
    return sim
end

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
            state = random_state(phase.type, state, phase.randomize)
        end
    end
    return Simulation(sim, state)
end


function run_phase(sim::Simulation, phase::ToMixed)
    if sim.type == Mixed 
        log_msg(sim, "State is already in mixed representation")
    else
        log_msg(sim, "Creating mixed representation with $(length(sim)) sites")
        sim = mix(sim; phase.limits...)
        log_msg(sim, "State is now in mixed representation")
    end
    return sim
end


function run_phase(sim::Simulation, phase::Evolve)
    time_start = phase.time_start
    if isnothing(time_start)
        time_start = sim.time
    end
    nsweeps = round(phase.time / phase.time_step)
    duration = phase.time_step * nsweeps
    time_stop = time_start + duration
    log_msg(sim, "Evolving state from simulation time $(time_start) to $(time_stop)")
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
            coefs, n_correct = phase.corrections, nsweeps, algo.order, algo.w, time_start,
            phase.limits.cutoff, phase.limits.maxdim,
            observer! = ApproxWObserver(sim, phase.measures, phase.measures_periodicity))
    else
        state = tdvp(pre, duration, state;
            coefs, n_expand = phase.corrections, nsweeps, time_start,
            phase.limits.cutoff, phase.limits.maxdim,
            observer! = TdvpWObserver(sim, phase.measures, phase.measures_periodicity))
    end
end




function run_phase(phase::Gates)
  log_msg("Applying $(length(phase.gates.ls)) gates")
  global sim_state = apply_gate(sim_type, sim_state, phase.gates; phase.limits.cutoff, phase.limits.maxdim)
end


function run_phase(phase::Dmrg)
  log_msg("Optimizing state with $(phase.nsweep) sweeps of Dmrg")
  if !simpleLit(phase.hamiltonian)
    die("Dmrg hamiltonian cannot contain multiple site operators")
  end
  mpo = make_mpo(tp, siteinds(sim_state), insertFfactors(reorder(phase.hamiltonian)))
  log_msg("MPO optimizer has maxdim $(maxlinkdim(mpo)) and uses $(Base.format_bytes(Base.summarysize(mpo)))")
  step_counter = 0
  optimize_time = time()
  dl = length(phase.maxdim)
  maxdim = 1
  while step_counter < phase.nsweep
    if step_counter < dl
      maxdim = phase.maxdim[1 + step_counter]
    end
    enrgy, st = dmrg(mpo, sim_state; nsweeps = 1, maxdim, cutoff=phase.cutoff)
    global sim_state = st
    step_counter += 1
    log_msg("after $step_counter dmrg sweep(s) state energy is $enrgy")
    if phase.output_periodicity ≠ 0 && mod(step_counter, phase.output_periodicity) == 0
      write_output(phase.data_output)
      tm = time()
      elapsed = round(tm - optimize_time; digits = 3)
      log_msg("$(phase.output_periodicity) steps took $elapsed seconds")
      optimize_time = tm
    end
  end
end
