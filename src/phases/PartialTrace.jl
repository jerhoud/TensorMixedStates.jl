# The PartialTrace phase, which traces out part of the sites.

export PartialTrace

"""
    PartialTrace(; positions, keep = false, name, time_start, final_measurements)

a phase that traces out the sites at `positions` or, with `keep = true`, all the others, see
`partial_trace`.

# Fields

- `name`, `time_start`, `final_measurements`: the fields every phase has, see `AbstractPhase`
- `positions`: the sites to trace out, or to keep
- `keep`: whether `positions` are the sites to keep rather than those to trace out (default
  `false`)

The state must be mixed, see `ToMixed`, and conserve nothing strongly, see `Weaken`. The sites
kept make a new system, numbered from 1 in their order, which the operators of the phases
after this one refer to.

# Examples

    PartialTrace(positions = [2, 3, 6])
    PartialTrace(positions = [1, 4, 5], keep = true)
"""
@kwdef struct PartialTrace <: AbstractPhase
    name::String = "Computing partial trace"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    positions::Vector{Int}
    keep::Bool = false
end

run_phase(sim::Simulation, phase::PartialTrace) = partial_trace(sim, phase.positions; phase.keep)
