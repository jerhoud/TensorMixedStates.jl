export random_state

"""
    random_state(Pure()|Mixed(), ::System, linkdims::Int)
    random_state(Pure()|Mixed(), ::Int, ::AbstractSite, linkdims::Int)
    random_state(Pure()|Mixed(), ::Vector{<:AbstractSite}, linkdims::Int)
    random_state(Pure()|Mixed(), ::State, linkdims::Int)

Return a random state, with the specified link dimension.
If a State is given, randomize the given state.
"""
function random_state(::Pure, system::System, linkdims::Int)
    st = random_mps(ComplexF64, system.pure_indices; linkdims)
    return State(Pure(), system, st)
end

function random_state(::Mixed, system::System, linkdims::Int)
    n = length(system)
    super = system ⊗ system
    super_rand = mix(random_state(Pure(), super, floor(Int, sqrt(linkdims))))
    return partial_trace(super_rand, collect(1:n); keepers = true)
end

random_state(type::PM, size::Int, site::AbstractSite, linkdims::Int) =
    random_state(type, System(size, site), linkdims)

random_state(type::PM, sites::Vector{<:AbstractSite}, linkdims::Int) =
    random_state(type, System(sites), linkdims)

function random_state(::Pure, state::State, linkdims::Int)
    st = random_mps(ComplexF64, state.system.pure_indices, state.state; linkdims)
    return State(Pure(), system, st)
end

random_state(::Mixed, state::State, linkdims::Int) =
    error("randomizing a mixed state is not implemented")

random_state(state::State, linkdims::Int) =
    random_state(state.type, state, linkdims)