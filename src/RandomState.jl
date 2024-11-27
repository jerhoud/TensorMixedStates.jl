export random_state

"""
    random_state(Pure|Mixed, ::System, linkdims::Int; start_time=0.0)
    random_state(Pure|Mixed, ::State, linkdims::Int; start_time=0.0)

Return a random state, with the specified link dimension.
If a State is given, randomize the given state.
The start_time argument sets the initial state simulation time.
"""
function random_state(::TPure, system::System, linkdims::Int; start_time::Float64=0.0)
    st = random_mps(ComplexF64, system.pure_sites; linkdims)
    return State(Pure, system, st, start_time)
end

function random_state(::TPure, state::State, linkdims::Int; start_time::Float64=state.time)
    if state.type ≠ Pure
        error("cannot randomize a mixed state into a pure state")
    end
    st = random_mps(ComplexF64, state.system.pure_sites, state.state; linkdims)
    return State(Pure, system, st, start_time)
end

function random_state(::TMixed, system::System, linkdims::Int; start_time::Float64=0.0)
    n = length(system)
    psites = system.pure_sites
    msites = system.mixed_sites
    super = System([psites; sim.(psites)], [msites, sim.(msites)])
    super_pure = random_state(Pure, super, 1 + linkdims ÷ 2)
    super_mixed = truncate(mix_state(super_pure), maxdim = linkdims, cutoff = 0)
    t = ITensor(1)
    for i in 2n:-1:n+1
        t *= mixed_obs(super_mixed, i)
    end
    super_mixed.state[n] *= t
    return State(Mixed, system, MPS(super_mixed.state[1:n]), start_time)
end

random_state(::TMixed, state::State, linkdims::Int; start_time::Float64=state.tim) =
    error("randomizing a mixed state is not implemented")
