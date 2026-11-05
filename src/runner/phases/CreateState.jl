# The CreateState phase, which gives the simulation a state built from a description, drawn at
# random, or given.

export CreateState

"""
    CreateState(; type, system, state, randomize, backend, on, name, time_start,
                  final_measurements)
    CreateState(state; options...)
    CreateState{Pure|Mixed}(system, state; options...)
    CreateState{Pure|Mixed}(n, site, state; options...)
    CreateState{Pure|Mixed}(sites, state; options...)

a phase that creates the state of the simulation.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `type`: the representation of the state, `Pure` or `Mixed`, by default that of `state` when
  it is a state; the type holding it is `backend`
- `system`: the `System` of the state, unused when `state` is a state
- `state`: a description of the state, see `State`, or a state, mixed if `type` asks for it
- `randomize`: the link dimension of a random state to draw with the generator of the run,
  see `random_state` (default 0, none): a pure state randomized from `state`, a mixed one drawn
  from a purification starting from its local states, which `Mixture` and `FullyMixed` do not
  give, or without `state` any state of the system, drawn as a `State` only
- `backend`: the type of states holding the state, `State` (MPS) or `ExactState`, by default that
  of the simulation, see `SimData`
- `on`: the name the state created takes, see `CopyState` (default `:state`)

# Examples

    CreateState(type = Pure, system = System(10, Qubit()), state = "Up", randomize = 50)
    CreateState(graph_state(System(10, Qubit()), complete_graph(10)))
    CreateState{Mixed}(4, Fermion(conserve = N), ["Occ", "Emp", "Occ", "Emp"]; randomize = 16)
    CreateState{Pure}(10, Qubit(), "Up")
    CreateState{Mixed}(4, Qubit(), "Up"; backend = ExactState)
"""
struct CreateState{R <: PM} <: AbstractPhase
    name::String
    time_start::Union{Nothing, Number}
    final_measurements
    type::Type{R}
    system::Union{Nothing, System}
    state
    randomize::Int
    backend::Union{Nothing, Type}
    on::Symbol
end

function CreateState(; name::String = "Creating state",
                     time_start::Union{Nothing, Number} = nothing, final_measurements = [],
                     type = nothing, system::Union{Nothing, System} = nothing, state = nothing,
                     randomize::Int = 0, backend = nothing, on::Symbol = :state)
    if type isa PM
        r = nameof(typeof(type))
        error("CreateState takes the type of its state, type = $r, where 1.x took type = $r()")
    elseif isnothing(type) && !(state isa AbstractState)
        error("CreateState needs the type of its state, type = Pure or type = Mixed, unless " *
              "it is given a state")
    end
    R = isnothing(type) ? representation(state) : type
    return CreateState{R}(name, time_start, final_measurements, R, system, state, randomize,
                          check_backend(backend), on)
end

# the fields of 2.0, before backend
CreateState{R}(name, time_start, final_measurements, type, system, state, randomize,
               on) where R =
    CreateState{R}(name, time_start, final_measurements, type, system, state, randomize, nothing,
                   on)
CreateState(name, time_start, final_measurements, type::Type{R}, system, state, randomize,
            on) where {R <: PM} =
    CreateState{R}(name, time_start, final_measurements, type, system, state, randomize, nothing,
                   on)

CreateState(state::AbstractState; kwargs...) = CreateState(; state, kwargs...)
CreateState{R}(system::System, state; kwargs...) where R =
    CreateState(; type = R, system, state, kwargs...)
CreateState{R}(n::Int, site::AbstractSite, state; kwargs...) where R =
    CreateState(; type = R, system = System(n, site), state, kwargs...)
CreateState{R}(sites::Vector{<:AbstractSite}, state; kwargs...) where R =
    CreateState(; type = R, system = System(sites), state, kwargs...)

provided_values(phase::CreateState) = phase.on === :state ? () : (phase.on,)

"""
    as_representation(R, state; limits = Limits())

`state` in representation `R`, for `CreateState` and `ConvertState`: a pure state is mixed if `R`
is `Mixed`, a mixed state refused if `R` is `Pure`; truncated to `limits` unless they are the
default
"""
as_representation(::Type{R}, state::State{R}; limits::Limits = Limits()) where R =
    limits == Limits() ? state : truncate(state; limits)
as_representation(::Type{Mixed}, state::State{Pure}; limits::Limits = Limits()) =
    mix(state; limits)
as_representation(::Type{Pure}, ::State{Mixed}; limits::Limits = Limits()) =
    error("a mixed state cannot be made pure: ask for Mixed rather than Pure")

"""
    as_backend(backend, state; limits = Limits())

`state` held by `backend`, a type of states as `State` or `ExactState`, in the same
representation, for `CreateState`, `ConvertState` and `LoadState`; truncated to `limits` by a
backend that approximates, unless they are the default. Two types of one's own convert through
`State`, or by a method naming both.
"""
as_backend(::Type{T}, state::T; limits::Limits = Limits()) where {T <: AbstractState} = state
as_backend(::Type{State}, state::State; limits::Limits = Limits()) =
    limits == Limits() ? state : truncate(state; limits)
as_backend(::Type{State}, state::AbstractState; limits::Limits = Limits()) = State(state; limits)
as_backend(::Type{T}, state::State; limits::Limits = Limits()) where {T <: AbstractState} =
    T(state; limits)
as_backend(::Type{T}, state::AbstractState; limits::Limits = Limits()) where {T <: AbstractState} =
    T(State(state); limits)

"""
the type of states holding `state`, without its representation, as `ExactState`
"""
backend_of(state::AbstractState) = Base.typename(typeof(state)).wrapper

function run_phase(phase::CreateState{R}, sim::Simulation, ::Int) where R
    rng = sim.rng
    n = phase.randomize
    backend = something(phase.backend, run_backend(sim))
    # an MPS drawn and converted is not the uniform draw an exact state is to have
    if n ≠ 0 && backend !== State
        error("CreateState draws a random state as a State only: draw it with backend = State " *
              "and convert it with ConvertState")
    end
    state = if isnothing(phase.state)
        if n == 0 || isnothing(phase.system)
            error("CreateState needs a state, or a system and randomize to draw one")
        end
        random_state(R, phase.system, n; rng)
    elseif phase.state isa AbstractState
        # a random mixed state is drawn from a purification, which a state given does not give
        if R === Mixed && n ≠ 0
            error("CreateState cannot randomize a state into a mixed state: give its local " *
                  "states, whose purification it is drawn from")
        end
        # the backend first, then the representation, as ConvertState does
        st = as_representation(R, as_backend(backend, phase.state))
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
    state = as_backend(backend, state)
    log_message(sim, "Created $state")
    return with_state(sim, phase.on, state)
end
