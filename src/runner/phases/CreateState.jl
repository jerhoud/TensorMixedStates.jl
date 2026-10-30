# The CreateState phase, which gives the simulation a state built from a description, drawn at
# random, or given.

export CreateState

"""
    CreateState(; type, system, state, randomize, on, name, time_start, final_measurements)
    CreateState(state::State; options...)
    CreateState{Pure|Mixed}(system, state; options...)
    CreateState{Pure|Mixed}(n, site, state; options...)
    CreateState{Pure|Mixed}(sites, state; options...)

a phase that creates the state of the simulation.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `type`: the representation of the state, `Pure` or `Mixed`, by default that of `state` when
  it is a `State`
- `system`: the `System` of the state, unused when `state` is a `State`
- `state`: a description of the state, see `State`, or a `State`, mixed if `type` asks for it
- `randomize`: the link dimension of a random state to draw with the generator of the run,
  see `random_state` (default 0, none): a pure state randomized from `state`, a mixed one drawn
  from a purification starting from it, or without `state` any state of the system
- `on`: the name the state created takes, see `CopyState` (default `:state`)

# Examples

    CreateState(type = Pure, system = System(10, Qubit()), state = "Up", randomize = 50)
    CreateState(graph_state(System(10, Qubit()), complete_graph(10)))
    CreateState{Mixed}(4, Fermion(conserve = N), ["Occ", "Emp", "Occ", "Emp"]; randomize = 16)
    CreateState{Pure}(10, Qubit(), "Up")
"""
struct CreateState{R <: PM} <: AbstractPhase
    name::String
    time_start::Union{Nothing, Number}
    final_measurements
    type::Type{R}
    system::Union{Nothing, System}
    state
    randomize::Int
    on::Symbol
end

function CreateState(; name::String = "Creating state",
                     time_start::Union{Nothing, Number} = nothing, final_measurements = [],
                     type = nothing, system::Union{Nothing, System} = nothing, state = nothing,
                     randomize::Int = 0, on::Symbol = :state)
    if type isa PM
        r = nameof(typeof(type))
        error("CreateState takes the type of its state, type = $r, where 1.x took type = $r()")
    elseif isnothing(type) && !(state isa State)
        error("CreateState needs the type of its state, type = Pure or type = Mixed, unless its " *
              "state is a State")
    end
    R = something(type, state isa State{Pure} ? Pure : Mixed)
    return CreateState{R}(name, time_start, final_measurements, R, system, state, randomize, on)
end

CreateState(state::State; kwargs...) = CreateState(; state, kwargs...)
CreateState{R}(system::System, state; kwargs...) where R =
    CreateState(; type = R, system, state, kwargs...)
CreateState{R}(n::Int, site::AbstractSite, state; kwargs...) where R =
    CreateState(; type = R, system = System(n, site), state, kwargs...)
CreateState{R}(sites::Vector{<:AbstractSite}, state; kwargs...) where R =
    CreateState(; type = R, system = System(sites), state, kwargs...)

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
    log_message(sim, "Created $state")
    return with_state(sim, phase.on, state)
end
