# The SaveState phase, which saves the state of the simulation in a file.

export SaveState

"""
    SaveState(; file, statename = "state", name, time_start, final_measurements)

a phase that saves the state in an HDF5 file, see `save_state`. Several states can be saved in
one file under different `statename`; saving under a name already in the file replaces it.

# Fields

- `name`, `time_start`, `final_measurements`: the fields every phase has, see `AbstractPhase`
- `file`: the name of the HDF5 file to write to, taken in the simulation directory when it is
  relative
- `statename`: the name under which the state is stored in the file

# Examples

    SaveState(file = "myfile.h5")
    SaveState(file = "myfile.h5", statename = "after_evolution")
"""
@kwdef struct SaveState <: AbstractPhase
    name::String = "Saving state"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    file::String
    statename::String = "state"
end

function run_phase(sim::Simulation, phase::SaveState)
    save_state(phase.file, phase.statename, sim)
    return sim
end
