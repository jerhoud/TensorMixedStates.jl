# The phases transforming the state at once: ToMixed, Weaken, PartialTrace and Gates.

export ToMixed, Weaken, PartialTrace, Gates

"""
    ToMixed(; limits, name, time_start, final_measurements)

a phase that takes the state to the mixed representation and truncates it to `limits`, see
`mix`; a state already mixed is only truncated.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `limits`: the truncation of the mixed state, see `Limits`

# Examples

    ToMixed(limits = Limits(cutoff = 1e-10, maxdim = 10))
"""
@kwdef struct ToMixed <: AbstractPhase
    name::String = "Switching to mixed state representation"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    limits::Limits = Limits()
end

function run_phase(phase::ToMixed, sim::Simulation, ::Int)
    if sim.state isa State{Mixed}
        log_message(sim, "The state is already mixed")
    end
    return Simulation(sim, truncate(mix(sim.state); phase.limits))
end

"""
    Weaken(; target, name, time_start, final_measurements)

a phase that takes the state down to a lower level of conservation, see `weaken`: a phase may
evolve under a strong symmetry, and the next one under a weak one, where a jump that moves the
charge becomes possible.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `target`: what the state must still conserve, as `weaken` takes it (default `nothing`, one
  level down)

# Examples

    Weaken()
    Weaken(target = ())
"""
@kwdef struct Weaken <: AbstractPhase
    name::String = "Weakening the symmetries"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    target = nothing
end

function run_phase(phase::Weaken, sim::Simulation, ::Int)
    from = symmetries(sim.state.system)
    st = isnothing(phase.target) ? weaken(sim.state) : weaken(sim.state, phase.target)
    log_message(sim, "Symmetries weakened from $from to $(symmetries(st.system))")
    return Simulation(sim, st)
end

"""
    PartialTrace(; positions, keep = false, name, time_start, final_measurements)

a phase that traces out the sites at `positions`, or with `keep = true` all the others, see
`partial_trace`. The sites kept make a new system, numbered from 1 in their order, which the
operators of the phases after this one refer to.

# Examples

    PartialTrace(positions = [2, 3, 6])
    PartialTrace(positions = [1, 4, 5], keep = true)
"""
@kwdef struct PartialTrace <: AbstractPhase
    name::String = "Computing partial trace"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    positions::Vector{Int} = nothing
    keep::Bool = false
    # the positions, which have no default, asked for by an error of their own: Julia attaches no
    # hint to the error of a keyword missing, which a keyword of 1.x left unread would raise
    function PartialTrace(name, time_start, final_measurements, positions, keep)
        if isnothing(positions)
            error("PartialTrace needs positions, the sites to trace out, or with keep = true " *
                  "those to keep")
        end
        return new(name, time_start, final_measurements, positions, keep)
    end
end

run_phase(phase::PartialTrace, sim::Simulation, ::Int) =
    Simulation(sim, partial_trace(sim.state, phase.positions; phase.keep))

"""
    Gates(; gates, limits, name, time_start, final_measurements)

a phase that applies gates to the state, see `apply`.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `gates`: the gates, an operator or a channel placed on sites
- `limits`: the truncation, see `Limits`

# Examples

    Gates(gates = controlled(X)(1, 3) * controlled(Z)(2, 4), limits = Limits(maxdim = 20))
"""
@kwdef struct Gates <: AbstractPhase
    name::String = "Applying gates"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    gates::IndexedOp
    limits::Limits = Limits()
end

run_phase(phase::Gates, sim::Simulation, ::Int) =
    Simulation(sim, apply(phase.gates, sim.state; phase.limits))
