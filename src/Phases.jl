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
  log("saving state to file $(phase.file)")
  file = h5open(phase.file, "cw")
  write(file, phase.object_name, sim_state)
  close(file)
  log("done")
end

function run_phase(phase::LoadState)
  log("loading state from file $(phase.file)")
  cd(start_dir) do
    file = h5open(phase.file, "r")
    global sim_state = read(file, phase.object_name, MPS)
    if phase.limits ≠ no_limits
      truncate!(mps; phase.limits.cutoff, phase.limits.maxdim)
    end
    close(file)
  end
  if isnothing(findfirst("Mixed", string(tags(siteind(sim_state, 1)))))
    log("loaded pure state with $(length(sim_state)) sites from file $(phase.file)")
    global sim_type = Pure
  else
    log("loaded mixed state with $(length(sim_state)) sites from file $(phase.file)")
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
  elseif isnothing(phase.state)
      global sim_state = random_mps(ComplexF64, sites; linkdims = phase.randomize)
    else
      global sim_state = random_mps(ComplexF64, sites, state; linkdims = phase.randomize)
  end
  clear_mpos()
end

function run_phase(phase::ToMixed)
  if sim_type == Mixed 
    log("State is already in mixed representation")
  else
    log("Creating mixed representation with $(length(sim_state)) sites")
    global sim_state = Pure2Mixed(sim_state; phase.limits.cutoff, phase.limits.maxdim)
    global sim_type = Mixed 
    clear_mpos()
    log("State is now in mixed representation")
  end
end

function run_phase(phase::Gates)
  log("Applying $(length(phase.gates.ls)) gates")
  global sim_state = apply_gate(sim_type, sim_state, phase.gates; phase.limits.cutoff, phase.limits.maxdim)
end

function run_phase(phase::Tdvp)
  stop_sim_time = sim_time + phase.duration
  log("Evolving state with Tdvp from simulation time $(sim_time) to $(stop_sim_time)")
  evolver = OpSum()
  if !isnothing(phase.hamiltonian)
    evolver += Lit_to_OpSum(phase.hamiltonian, if sim_type == Pure "" else "oper" end)
  end
  if !isnothing(phase.dissipator)
    if sim_type == Pure
      die("Tdvp evolving error: state must be in mixed representation to use dissipators")
    end
    evolver += Lit_to_OpSum(phase.dissipator)
  end
  mpo = MPO(evolver, siteinds(sim_state))
  log("MPO evolver has maxdim $(maxlinkdim(mpo)) and uses $(Base.format_bytes(Base.summarysize(mpo)))")
  output_counter = 0
  evolve_time = time()
  while sim_time < stop_sim_time
    global sim_state = tdvp(mpo, phase.tau, sim_state; nsteps = 1, time_step = phase.tau)
    output_counter += 1
    global sim_time += phase.tau
    if output_counter == phase.output_periodicity
      write_output(phase.data_output)
      tm = time()
      elapsed = round(tm - evolve_time; digits = 3)
      log("$(phase.output_periodicity) steps took $elapsed seconds")
      evolve_time = tm
      output_counter = 0
    end
  end
end

function run_phase(phase::Dmrg)
  if sim_type == Mixed
    die("Cannot perform dmrg on state with mixed representation")
  end
  log("Optimizing state with $(phase.nsweep) sweeps of Dmrg")
  hamiltonian = Lit_to_OpSum(phase.hamiltonian)
  mpo = MPO(hamiltonian, siteinds(sim_state))
  log("MPO optimizer has maxdim $(maxlinkdim(mpo)) and uses $(Base.format_bytes(Base.summarysize(mpo)))")
  output_counter = 0
  optimize_time = time()
  dl = length(phase.maxdim) 
  i = 1
  while i <= phase.nsweep
    global sim_state = dmrg(mpo, sim_state; nsweeps = 1, maxdim = phase.maxdim[if i <= dl i else dl end], cutoff=phase.cutoff)
    i += 1
    output_counter += 1
    if output_counter == phase.output_periodicity
      write_output(phase.data_output)
      tm = time()
      elapsed = round(tm - optimize_time; digits = 3)
      log("$(phase.output_periodicity) steps took $elapsed seconds")
      optimize_time = tm
      output_counter = 0
    end
  end

  hamiltonian::Union{ProdLit, SumLit}
  cutoff::Float64
  maxdim::Union{Int, Vector{Int}}
  output_periodicity::Int = 1
  data_output::Output = no_output
end