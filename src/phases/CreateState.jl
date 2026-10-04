# The CreateState phase, which gives the simulation a state built from a description, drawn at
# random, or given.

export CreateState

"""
    CreateState(; type, system, state, randomize, seed, name, time_start, final_measurements)
    CreateState{Pure|Mixed}(n, site, state; options...)
    CreateState{Pure|Mixed}(sites, state; options...)

a phase that creates the state of the simulation. The first phase of a simulation is this one
or `LoadState`.

# Fields

- `name`, `time_start`, `final_measurements`: the fields every phase has, see `AbstractPhase`
- `type`: the representation of the state, `Pure()` or `Mixed()`, or a `Representation` an
  extension defines, together with the method of `run_phase` creating its state
- `system`: the `System` of the state, unused when `state` is a `State`
- `state`: a description of the state, or a `State`, which is mixed if `type` asks for it (a
  mixed one cannot be made pure)
- `randomize`: the link dimension of a random state to create (default 0, none). With a
  `state`, a pure state is randomized from it and a mixed one drawn from a purification
  starting from it, see `RandomState`; a `State` can only be randomized into a pure state
- `seed`: the seed the global random generator is given at the start of the phase (default
  `nothing`, none). A checkpoint does not save the generator, so a run resumed after this
  phase does not set it again

# Examples

    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up")
    CreateState(type = Mixed(), system = System(3, Qubit()), state = ["Up", "Dn", "Up"])
    CreateState(type = Pure(), system = System(10, Qubit()), randomize = 50)
    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up", randomize = 50)
    CreateState{Mixed}(4, Fermion(conserve = N), ["Occ", "Emp", "Occ", "Emp"]; randomize = 16)
    CreateState{Pure}(10, Qubit(), "Up")                                      # simple form
    CreateState{Mixed}([Qubit(), Boson(4), Fermion()], ["Up", "2", "Occ"])    # other simple form
"""
@kwdef struct CreateState{R <: Representation} <: AbstractPhase
    name::String = "Creating state"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    type::R
    system::Union{Nothing, System} = nothing
    state = nothing
    randomize::Int = 0
    seed::Union{Nothing, Int} = nothing
end

CreateState{R}(n, site, state; kwargs...) where R =
    CreateState(;type = R(), system = System(n, site), state, kwargs...)
CreateState{R}(sites, state; kwargs...) where R =
    CreateState(;type = R(), system = System(sites), state, kwargs...)

"""
    as_representation(sim, R, ::State)

the given state in representation `R`, for a `CreateState` handed a `State` rather than a
description. A pure state is mixed if `R` is `Mixed`; a mixed state is refused if `R` is
`Pure`, since it holds no purification to go back to.
"""
as_representation(::Simulation, ::Type{R}, state::State{R}) where R = state
function as_representation(sim::Simulation, ::Type{Mixed}, state::State{Pure})
    log_msg(sim, "Creating mixed representation with $(length(state)) sites")
    return mix(state)
end
as_representation(::Simulation, ::Type{Pure}, ::State{Mixed}) =
    error("CreateState was asked for a pure state but given a mixed one, which cannot be " *
          "turned back into a pure state")

function run_phase(sim::Simulation, phase::CreateState{R}) where {R <: PM}
    if !isnothing(phase.seed)
        Random.seed!(phase.seed)
    end
    if isnothing(phase.state)
        if phase.randomize == 0
            error("CreateState without state nor randomize: no state created !")
        elseif isnothing(phase.system)
            error("CreateState needs a system to create a random state")
        else
            state = RandomState{R}(phase.system, phase.randomize)
        end
    elseif phase.state isa State
        # a random mixed state is drawn from a purification, which a State does not give
        if R === Mixed && phase.randomize ≠ 0
            error("CreateState cannot randomize a State into a mixed state: give a description " *
                  "of the state, whose purification it is drawn from")
        end
        # `type` is what the phase was asked for, so a State given in the other
        # representation is converted rather than silently kept as it is
        state = as_representation(sim, R, phase.state)
        if phase.randomize ≠ 0
            state = RandomState(state, phase.randomize)
        end
    elseif isnothing(phase.system)
        error("CreateState needs a system or a State object")
    elseif phase.randomize == 0
        state = State{R}(phase.system, phase.state)
    elseif R === Mixed
        # a mixed state is drawn from the states its purification starts from, the only way
        # there is on a system that conserves something
        state = RandomState{Mixed}(phase.system, phase.state, phase.randomize)
    else
        state = RandomState(State{Pure}(phase.system, phase.state), phase.randomize)
    end
    return Simulation(sim, state)
end

creates_state(::CreateState) = true
phase_system(p::CreateState) = p.state isa State ? p.state.system : p.system
