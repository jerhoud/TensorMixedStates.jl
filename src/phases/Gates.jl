# The Gates phase, which applies gates to the state.

export Gates

"""
    Gates(; gates, limits, name, time_start, final_measurements)

a phase that applies gates to the state.

# Fields

- `name`, `time_start`, `final_measurements`: the fields every phase has, see `AbstractPhase`
- `gates`: the gates to apply
- `limits`: the truncations made while applying a gate of several sites, see `apply` (default
  `Limits()`)

# Examples

    Gates(gates = controlled(X)(1, 3)*controlled(Z)(2, 4), limits = Limits(cutoff=1e-10, maxdim = 20))
"""
@kwdef struct Gates <: AbstractPhase
    name::String = "Applying gates"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    gates::IndexedOp
    limits::Limits = Limits()
end

function run_phase(sim::Simulation, phase::Gates)
    log_msg(sim, "Applying gates")
    return apply(phase.gates, sim; phase.limits)
end
