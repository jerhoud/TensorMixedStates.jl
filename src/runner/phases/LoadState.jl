# The phases reading and writing states: LoadState, from a file, SaveState, to a file, and
# CopyState, between the names of the simulation.

export LoadState, SaveState, CopyState

"""
    in_dir(sim, file)

the path of `file`, relative to the directory of the simulation
"""
in_dir(sim::Simulation, file::String) = joinpath(getfield(sim, :run).dir, file)

"""
    LoadState(; file, statename = "state", limits, on = :state, name, time_start,
                final_measurements)

a phase that loads a state from an HDF5 file written by `SaveState` or `save_state`, see
`load_state`, and truncates it to `limits`, unless they are the default `Limits()`.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `file`: the HDF5 file, taken in the directory of the simulation when it is relative: a state
  another simulation saved is found under `../othername/`
- `statename`: the name of the state in the file
- `limits`: the truncation of the state once loaded, see `Limits`
- `on`: the name the state takes in the simulation, `:state` for the current state, another
  for a reference state, `c.ref` to the measurements, see `Fidelity`

# Examples

    LoadState(file = "ground.h5")
    LoadState(file = "../other/ground.h5", statename = "gs", on = :ground)
"""
@kwdef struct LoadState <: AbstractPhase
    name::String = "Loading state"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    file::String
    statename::String = "state"
    limits::Limits = Limits()
    on::Symbol = :state
end

provided_values(phase::LoadState) = phase.on === :state ? () : (phase.on,)

function run_phase(phase::LoadState, sim::Simulation, ::Int)
    st = load_state(in_dir(sim, phase.file), phase.statename)
    # not truncated at the default limits, which a state of an extension need not support
    return with_state(sim, phase.on, phase.limits == Limits() ? st : truncate(st; phase.limits))
end

"""
    SaveState(; file, statename = "state", on = :state, name, time_start, final_measurements)

a phase that saves a state in an HDF5 file, see `save_state`, replacing a state saved there
under the same name.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `file`: the HDF5 file, taken in the directory of the simulation when it is relative
- `statename`: the name of the state in the file
- `on`: the state saved, `:state` the current one, or one the simulation named

# Examples

    SaveState(file = "states.h5", statename = "after_quench")
"""
@kwdef struct SaveState <: AbstractPhase
    name::String = "Saving state"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    file::String
    statename::String = "state"
    on::Symbol = :state
end

function run_phase(phase::SaveState, sim::Simulation, ::Int)
    if !isempty(getfield(sim, :run).dir) && normpath(phase.file) in simulation_files
        error("cannot save to $(phase.file), a file of the simulation directory: choose " *
              "another name")
    end
    save_state(in_dir(sim, phase.file), phase.statename, state_named(sim, phase.on))
    return sim
end

"""
    CopyState(from => to; move = false, name, time_start, final_measurements)

a phase that copies the state named `from` under the name `to`, `:state` naming the current
state: `CopyState(:state => :ground)` keeps the current state as `:ground`, which the
measurements read as `c.ground`, see `Fidelity`, and `CopyState(:ground => :state)` makes it
the current state again. With `move = true`, `from` is removed, which frees its memory and its
place in the checkpoints. A copy costs nothing, the two names sharing the state.

# Examples

    CopyState(:state => :ground)
    CopyState(:ground => :state; move = true)
"""
@kwdef struct CopyState <: AbstractPhase
    name::String = "Copying state"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    names::Pair{Symbol, Symbol}
    move::Bool = false
    function CopyState(name, time_start, final_measurements, names, move)
        from, to = names
        for n in (from, to)
            if n ≠ :state && n in reserved_values
                error("CopyState cannot name a state $n, a name the run reserves")
            end
        end
        if from == to
            error("CopyState copies a state under another name, not under $to")
        end
        return new(name, time_start, final_measurements, names, move)
    end
end

CopyState(names::Pair{Symbol, Symbol}; kwargs...) = CopyState(; names, kwargs...)

provided_values(phase::CopyState) = last(phase.names) === :state ? () : (last(phase.names),)

function run_phase(phase::CopyState, sim::Simulation, ::Int)
    from, to = phase.names
    sim = with_state(sim, to, state_named(sim, from))
    if !phase.move
        return sim
    elseif from === :state
        return Simulation(sim, nothing)
    end
    values = copy(getfield(sim, :values))
    delete!(values, from)
    return Simulation(sim.state, sim.time, values, getfield(sim, :run))
end
