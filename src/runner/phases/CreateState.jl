# The CreateState phase, which gives the simulation a state built from a description, drawn at
# random, or given.

export CreateState

"""
    CreateState(; type, system, state, randomize, on, name, time_start, final_measurements)
    CreateState{Pure|Mixed}(n, site, state; options...)
    CreateState{Pure|Mixed}(sites, state; options...)

a phase that creates the state of the simulation.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `type`: the representation of the state, `Pure()` or `Mixed()`
- `system`: the `System` of the state, unused when `state` is a `State`
- `state`: a description of the state, see `State`, or a `State`, mixed if `type` asks for it
- `randomize`: the link dimension of a random state to draw with the generator of the run,
  see `random_state` (default 0, none): a pure state randomized from `state`, a mixed one drawn
  from a purification starting from it, or without `state` any state of the system
- `on`: the name the state created takes, see `CopyState` (default `:state`)

# Examples

    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up")
    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up", randomize = 50)
    CreateState{Mixed}(4, Fermion(conserve = N), ["Occ", "Emp", "Occ", "Emp"]; randomize = 16)
    CreateState{Pure}(10, Qubit(), "Up")
"""
@kwdef struct CreateState{R <: PM} <: AbstractPhase
    name::String = "Creating state"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    type::R
    system::Union{Nothing, System} = nothing
    state = nothing
    randomize::Int = 0
    on::Symbol = :state
end

CreateState{R}(n::Int, site::AbstractSite, state; kwargs...) where R =
    CreateState(; type = R(), system = System(n, site), state, kwargs...)
CreateState{R}(sites::Vector{<:AbstractSite}, state; kwargs...) where R =
    CreateState(; type = R(), system = System(sites), state, kwargs...)

provided_values(phase::CreateState) = phase.on === :state ? () : (phase.on,)

"""
    as_representation(R, state)

the `State` handed to a `CreateState` in representation `R`: a pure state is mixed if `R` is
`Mixed`, a mixed state refused if `R` is `Pure`
"""
as_representation(::Type{R}, state::State{R}) where R = state
as_representation(::Type{Mixed}, state::State{Pure}) = mix(state)
as_representation(::Type{Pure}, ::State{Mixed}) =
    error("CreateState was asked for a pure state but given a mixed one")

function run_phase(phase::CreateState{R}, sim::Simulation, ::Int) where R
    rng = sim.rng
    n = phase.randomize
    state = if isnothing(phase.state)
        if n == 0 || isnothing(phase.system)
            error("CreateState needs a state, or a system and randomize to draw one")
        end
        random_state(R, phase.system, n; rng)
    elseif phase.state isa State
        # a random mixed state is drawn from a purification, which a State does not give
        if R === Mixed && n ≠ 0
            error("CreateState cannot randomize a State into a mixed state: give a description " *
                  "of the state, whose purification it is drawn from")
        end
        st = as_representation(R, phase.state)
        n == 0 ? st : random_state(st, n; rng)
    elseif isnothing(phase.system)
        error("CreateState needs a system or a State")
    elseif n == 0
        State{R}(phase.system, phase.state)
    elseif R === Mixed
        random_state(Mixed, phase.system, phase.state, n; rng)
    else
        random_state(State{Pure}(phase.system, phase.state), n; rng)
    end
    return with_state(sim, phase.on, state)
end
