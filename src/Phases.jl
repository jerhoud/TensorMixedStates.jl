function run_phases(phase)
  start_phase(phase.name)
  time_data = @timed begin
    run_phase(phase)
    write_output(phase.output)
  end
  end_phase(phase.name, time_data)
end

run_phases(phases::Union{Vector, Tuple}) = 
  foreach(run_phases, phases)

function run_phase(phase::SaveState)
  log_msg("saving state to file $(phase.file)")
  file = h5open(phase.file, "cw")
  write(file, phase.object_name, sim_state)
  close(file)
  log_msg("done")
end

function run_phase(phase::LoadState)
  log_msg("loading state from file $(phase.file)")
  cd(start_dir) do
    file = h5open(phase.file, "r")
    global sim_state = read(file, phase.object_name, MPS)
    if phase.limits ≠ no_limits
      truncate!(mps; phase.limits.cutoff, phase.limits.maxdim)
    end
    close(file)
  end
  if isnothing(findfirst("Mixed", string(tags(siteind(sim_state, 1)))))
    log_msg("loaded pure state with $(length(sim_state)) sites from file $(phase.file)")
    global sim_type = Pure
  else
    log_msg("loaded mixed state with $(length(sim_state)) sites from file $(phase.file)")
    global sim_type = Mixed
  end
  clear_mpos()
end


function run_phase(phase::CreateState)
  global sim_type = phase.type
  sites = [ st(phase.type) for st in phase.system ]
  if phase.randomize == 0
    if isnothing(phase.state)
      die("CreateState without state nor randomize: no state created !")
    else
      state = phase.state
      if phase.state isa String
        state = fill(state, length(sites))
      end
      if length(sites) ≠ length(state)
        die("Incompatible size between system($(length(sites))) and state($(length(state))) in CreateState")
      end
    end
    global sim_state = MPS(ComplexF64, sites, state)
  elseif sim_type == Mixed
      global sim_state = random_mixed_state(sites, phase.randomize)
  elseif isnothing(phase.state)
      global sim_state = random_mps(ComplexF64, sites; linkdims = phase.randomize)
    else
      global sim_state = random_mps(ComplexF64, sites, state; linkdims = phase.randomize)
  end
  clear_mpos()
end

function run_phase(phase::ToMixed)
  if sim_type == Mixed 
    log_msg("State is already in mixed representation")
  else
    log_msg("Creating mixed representation with $(length(sim_state)) sites")
    global sim_state = Pure2Mixed(sim_state; phase.limits.cutoff, phase.limits.maxdim)
    global sim_type = Mixed 
    clear_mpos()
    log_msg("State is now in mixed representation")
  end
end

function run_phase(phase::Gates)
  global sim_operator_mode = GateMode
  log_msg("Applying $(length(phase.gates.ls)) gates")
  global sim_state = apply_gate(sim_state, phase.gates; phase.limits.cutoff, phase.limits.maxdim)
end

function run_phase(phase::Tdvp)
  global sim_operator_mode = EvolveMode
  stop_sim_time = sim_time + phase.duration
  log_msg("Evolving state with Tdvp from simulation time $(sim_time) to $(stop_sim_time)")
  evolver = OpSum()
  if !isnothing(phase.hamiltonian)
    evolver += Lit_to_OpSum(-im * phase.hamiltonian)
  end
  if !isnothing(phase.dissipator)
    if sim_type == Pure
      die("Tdvp evolving error: state must be in mixed representation to use dissipators")
    end
    evolver += Lit_to_OpSum(phase.dissipator)
  end
  mpo = MPO(evolver, siteinds(sim_state))
  log_msg("MPO evolver has maxdim $(maxlinkdim(mpo)) and uses $(Base.format_bytes(Base.summarysize(mpo)))")
  step_counter = 0
  evolve_time = time()
  while sim_time < stop_sim_time
    global sim_state = tdvp(mpo, phase.tau, sim_state; nsteps = 1, time_step = phase.tau)
    step_counter += 1
    global sim_time += phase.tau
    if phase.expand ≠ 0 && mod(step_counter, phase.expand) == 0
      tm = time()
      global sim_state = expand(sim_state, mpo; alg="global_krylov")
      elapsed = round(time() - tm; digits = 3)
      log_msg("Krylov expand step took $elapsed seconds")
    end
    if phase.output_periodicity ≠ 0 && mod(step_counter, phase.output_periodicity) == 0
      write_output(phase.data_output)
      tm = time()
      elapsed = round(tm - evolve_time; digits = 3)
      log_msg("$(phase.output_periodicity) steps took $elapsed seconds")
      evolve_time = tm
    end
  end
end

function run_phase(phase::Dmrg)
  global sim_operator_mode = EvolveMode
  log_msg("Optimizing state with $(phase.nsweep) sweeps of Dmrg")
  hamiltonian = Lit_to_OpSum(phase.hamiltonian)
  mpo = MPO(hamiltonian, siteinds(sim_state))
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
