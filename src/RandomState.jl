export random_state

"""
    random_state(Pure|Mixed, ::System, linkdims::Int)
    random_state(Pure|Mixed, ::State, linkdims::Int)

Return a random state, with the specified link dimension.
If a State is given, randomize the given state.
"""
function random_state(::Pure, system::System, linkdims::Int)
    st = random_mps(ComplexF64, system.pure_indices; linkdims)
    return State(Pure(), system, st)
end

function random_state(::Mixed, system::System, linkdims::Int)
    n = length(system)
    sites = system.sites
    pind = system.pure_indices
    mind = system.mixed_indices
    super = Sytem([ sites ; sites], [pind ; sim.(pind)], [mind ; sim.(mind)])
    super_pure = random_state(Pure(), super, floor(Int, sqrt(linkdims)))
    super_mixed = mix(super_pure)
    
    t = ITensor(1)
    for i in 2n:-1:n+1
        t *= tensor_obs(super_mixed, i)
    end
    super_mixed.state[n] *= t
    return State(Mixed(), system, MPS(super_mixed.state[1:n]))
end

function random_state(::Pure, state::State, linkdims::Int)
    st = random_mps(ComplexF64, state.system.pure_indices, state.state; linkdims)
    return State(Pure(), system, st)
end

random_state(::Mixed, state::State, linkdims::Int) =
    error("randomizing a mixed state is not implemented")

random_state(state::State, linkdims::Int) =
    random_state(state.type, state, linkdims)