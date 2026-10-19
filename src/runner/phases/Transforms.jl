# The phases transforming the state: ToMixed, Weaken and PartialTrace at once, Gates step by step.

export ToMixed, Weaken, PartialTrace, Gates

"""
    ToMixed(; limits, on, name, time_start, final_measurements)

a phase that takes the state to the mixed representation and truncates it to `limits`, see
`mix`; a state already mixed is only truncated.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `limits`: the truncation of the mixed state, see `Limits`
- `on`: the name of the state the phase acts on, see `CopyState` (default `:state`)

# Examples

    ToMixed(limits = Limits(cutoff = 1e-10, maxdim = 10))
"""
@kwdef struct ToMixed <: AbstractPhase
    name::String = "Switching to mixed state representation"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    limits::Limits = Limits()
    on::Symbol = :state
end

function run_phase(phase::ToMixed, sim::Simulation, ::Int)
    st = state_named(sim, phase.on)
    if st isa State{Mixed}
        log_message(sim, "The state is already mixed")
    end
    return with_state(sim, phase.on, truncate(mix(st); phase.limits))
end

"""
    Weaken(; target, on, name, time_start, final_measurements)

a phase that takes the state down to a lower level of conservation, see `weaken`: a phase may
evolve under a strong symmetry, and the next one under a weak one, where a jump that moves the
charge becomes possible.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `target`: what the state must still conserve, as `weaken` takes it (default `nothing`, one
  level down)
- `on`: the name of the state the phase acts on, see `CopyState` (default `:state`)

# Examples

    Weaken()
    Weaken(target = ())
"""
@kwdef struct Weaken <: AbstractPhase
    name::String = "Weakening the symmetries"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    target = nothing
    on::Symbol = :state
end

function run_phase(phase::Weaken, sim::Simulation, ::Int)
    state = state_named(sim, phase.on)
    from = symmetries(state.system)
    st = isnothing(phase.target) ? weaken(state) : weaken(state, phase.target)
    log_message(sim, "Symmetries weakened from $from to $(symmetries(st.system))")
    return with_state(sim, phase.on, st)
end

"""
    PartialTrace(; positions, keep = false, on, name, time_start, final_measurements)

a phase that traces out the sites at `positions`, or with `keep = true` all the others, see
`partial_trace`. The sites kept make a new system, numbered from 1 in their order, which the
operators of the phases after this one refer to. `on` is the name of the state it acts on, see
`CopyState` (default `:state`).

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
    on::Symbol = :state
    # the positions, which have no default, asked for by an error of their own: Julia attaches no
    # hint to the error of a keyword missing, which a keyword of 1.x left unread would raise
    function PartialTrace(name, time_start, final_measurements, positions, keep, on)
        if isnothing(positions)
            error("PartialTrace needs positions, the sites to trace out, or with keep = true " *
                  "those to keep")
        end
        return new(name, time_start, final_measurements, positions, keep, on)
    end
end

run_phase(phase::PartialTrace, sim::Simulation, ::Int) =
    with_state(sim, phase.on,
               partial_trace(state_named(sim, phase.on), phase.positions; phase.keep))

"""
    Gates(; gates, limits, on, name, time_start, final_measurements)

a phase that applies gates to the state, see `apply`, by its steps, a batch of gates or an MPO,
between which the run can write a checkpoint or stop.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `gates`: the gates, an operator or a channel placed on sites
- `limits`: the truncation, see `Limits`
- `on`: the name of the state the phase acts on, see `CopyState`, or a vector of names, each
  state acted on alike (default `:state`)

# Examples

    Gates(gates = controlled(X)(1, 3) * controlled(Z)(2, 4), limits = Limits(maxdim = 20))
"""
@kwdef struct Gates <: AbstractPhase
    name::String = "Applying gates"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    gates::IndexedOp
    limits::Limits = Limits()
    on::Union{Symbol, Vector{Symbol}} = :state
end

"""
    GatesSteps(coefs, steps, limits, on)

a `Gates` prepared: for each state it acts on, the number it is multiplied by and the steps of
`apply`, see `apply_steps`, made in turn, its limits and the names of the states
"""
struct GatesSteps
    coefs::Vector{Number}
    steps::Vector{Vector}
    limits::Limits
    on::Vector{Symbol}
end

function prepare_phase(phase::Gates, sim::Simulation)
    names = state_names(phase.on)
    prepared = [ apply_steps(phase.gates, state_named(sim, name)) for name in names ]
    return GatesSteps(first.(prepared), last.(prepared), phase.limits, names)
end

# a product that vanishes, of no step, takes one to give the null state
step_count(g::GatesSteps) = max(1, maximum(length, g.steps))

# each state multiplied by its number at its last step, as `applied` does
function run_step(g::GatesSteps, sim::Simulation, k::Int)
    for (name, coef, steps) in zip(g.on, g.coefs, g.steps)
        state = state_named(sim, name)
        st = state.state
        if coef == 0
            if k == 1
                st = 0 * st
            end
        else
            if k ≤ length(steps)
                st = applied_step(st, steps[k], truncation(g.limits))
            end
            if k == max(1, length(steps)) && coef ≠ 1
                st = coef * st
            end
        end
        sim = with_state(sim, name, State(state, st))
    end
    return sim
end
