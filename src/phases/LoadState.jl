# The LoadState phase, which gives the simulation a state read from a file.

export LoadState

"""
    LoadState(; file, statename = "state", limits, name, time_start, final_measurements)

a phase that loads the state from an HDF5 file written by `SaveState` or `save_state`, see
`load_state`, and truncates it to `limits`, unless they are the default `Limits()`.

# Fields

- `name`, `time_start`, `final_measurements`: the fields every phase has, see `AbstractPhase`
- `file`: the name of the HDF5 file to read from, taken in the simulation directory when it is
  relative: a state another simulation saved is found under `../othername/`
- `statename`: the name under which the state is stored in the file
- `limits`: the truncation applied to the state once loaded, see `Limits` (default `Limits()`,
  none)

# Examples

    LoadState(file = "myfile.h5")
    LoadState(file = "myfile.h5", statename = "after_evolution")
"""
@kwdef struct LoadState <: AbstractPhase
    name::String = "Loading state"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    file::String
    statename::String = "state"
    limits::Limits = Limits()
end

function run_phase(sim::Simulation, phase::LoadState)
    st = load_state(phase.file, phase.statename)
    # a state read back is not truncated without limits of its own, which a representation of
    # one's own need not support
    return Simulation(sim, phase.limits == Limits() ? st : truncate(st; phase.limits))
end

creates_state(::LoadState) = true
