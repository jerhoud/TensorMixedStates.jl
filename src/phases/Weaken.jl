# The Weaken phase, which takes the state down to a lower level of conservation.

export Weaken

"""
    Weaken(; target, name, time_start, final_measures)

a phase that takes the state down to a lower level of conservation, see `weaken`. A phase may
evolve under a strong symmetry, which every dissipator commuting with the charge allows, and
the next one continue under a weak one, where a jump that moves the charge becomes possible.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `AbstractPhase`
- `target`: what the state must still conserve, as `weaken` takes it (default `nothing`, one
  level down: every strong quantity made weak or, when none is strong, every quantity
  dropped)

# Examples

    Weaken()
    Weaken(target = (strong(Ntot), 2Sz))
    Weaken(target = ())                     # no charges at all
"""
@kwdef struct Weaken <: AbstractPhase
    name::String = "Weakening the symmetries"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    target = nothing
end

function run_phase(sim::Simulation, phase::Weaken)
    from = symmetries(sim.state.system)
    sim = isnothing(phase.target) ? weaken(sim) : weaken(sim, phase.target)
    log_msg(sim, "Symmetries weakened from $from to $(symmetries(sim.state.system))")
    return sim
end
