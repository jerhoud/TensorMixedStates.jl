# The ToMixed phase, which takes the state to the mixed representation.

export ToMixed

"""
    ToMixed(; limits, name, time_start, final_measurements)

a phase that switches the state to the mixed representation and truncates it to `limits`; a
state already mixed is only truncated.

# Fields

- `name`, `time_start`, `final_measurements`: the fields every phase has, see
  `TensorMixedStates.AbstractPhase`
- `limits`: constraints on the mixed state, see `Limits` (default `Limits()`)

# Examples

    ToMixed()
    ToMixed(limits = Limits(cutoff = 1e-10, maxdim = 10))
"""
@kwdef struct ToMixed <: AbstractPhase
    name::String = "Switching to mixed state representation"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    limits::Limits = Limits()
end

function run_phase(sim::Simulation, phase::ToMixed)
    if sim.state isa State{Mixed}
        log_message(sim, "State is already in mixed representation")
        sim = truncate(sim; phase.limits)
    else
        log_message(sim, "Creating mixed representation with $(length(sim)) sites")
        sim = truncate(mix(sim); phase.limits)
        log_message(sim, "State is now in mixed representation")
    end
    return sim
end
