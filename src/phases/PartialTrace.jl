# The PartialTrace phase, which traces out part of the sites.

export PartialTrace

"""
    PartialTrace(; trace_positions | keep_positions, name, time_start, final_measures)

a phase that traces out part of the sites, given by exactly one of `trace_positions` and
`keep_positions`.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `AbstractPhase`
- `trace_positions`: the sites to trace out
- `keep_positions`: the sites to keep, all the others being traced out

The state must be mixed, see `ToMixed`, and conserve nothing strongly, see `Weaken`. The sites
kept make a new system, numbered from 1 in their order, which the operators of the phases
after this one refer to.

# Examples

    PartialTrace(trace_positions = [2, 3, 6])
    PartialTrace(keep_positions = [1, 4, 5])
"""
@kwdef struct PartialTrace <: AbstractPhase
    name::String = "Computing partial trace"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    trace_positions::Union{Nothing, Vector{Int}} = nothing
    keep_positions::Union{Nothing, Vector{Int}} = nothing
end

function run_phase(sim::Simulation, phase::PartialTrace)
    pos = phase.trace_positions
    keep = phase.keep_positions
    if isnothing(pos) == isnothing(keep)
        error("PartialTrace requires one and only one of trace_positions and keep_positions")
    end
    return isnothing(pos) ? partial_trace(sim, keep; keepers = true) : partial_trace(sim, pos)
end
