export RandomState

"""
    random_state(Pure()|Mixed(), ::System, linkdims::Int)
    random_state(Pure()|Mixed(), ::Int, ::AbstractSite, linkdims::Int)
    random_state(Pure()|Mixed(), ::Vector{<:AbstractSite}, linkdims::Int)
    random_state(Pure()|Mixed(), ::State, linkdims::Int)

Return a random state, with the specified link dimension.
If a State is given, randomize the given state.
"""
struct RandomState{R <: PM} end

function RandomState{Pure}(system::System, linkdims::Int)
    st = random_mps(ComplexF64, system.pure_indices; linkdims)
    return State{Pure}(system, st)
end

function RandomState{Mixed}(system::System, linkdims::Int)
    n = length(system)
    super = system ⊗ system
    super_rand = mix(RandomState{Pure}(super, floor(Int, sqrt(linkdims))))
    return partial_trace(super_rand, collect(1:n); keepers = true)
end

RandomState{R}(size::Int, site::AbstractSite, linkdims::Int) where R =
    RandomState{R}(System(size, site), linkdims)

RandomState{R}(sites::Vector{<:AbstractSite}, linkdims::Int) where R =
    RandomState{R}(System(sites), linkdims)

function RandomState(state::State{Pure}, linkdims::Int)
    st = random_mps(ComplexF64, state.system.pure_indices, state.state; linkdims)
    return State{Pure}(system, st)
end

RandomState(::State{Mixed}, ::Int) =
    error("randomizing a mixed state is not implemented")
