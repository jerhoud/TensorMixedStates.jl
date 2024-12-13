import Base: length, copy
import ITensors: truncate
import ITensorMPS: maxlinkdim, apply

export State, mix, truncate, maxlinkdim, memory_usage, apply

"""
    struct State
    State(Pure|Mixed, system, states)
    State(::State, ::MPS)


Represent the complete state of the simulated quantum system

# Fields
- `type::Union{Type{Pure}, Type{Mixed}}`: pure or mixed representation
- `system::System`: system description
- `state::MPS`: system state

# Examples
    State(Pure, system, "Up")
    State(Mixed, system, ["Up", "Dn", "Up"])

"""
struct State
    type::Union{TPure, TMixed}
    system::System
    state::MPS
end

"""
    length(::State)

Return the number of sites in the state
"""
length(state::State) = length(state.system)

"""
    maxlinkdim(::State)

Return the maximum link dimension in the state
"""
maxlinkdim(state::State) = maxlinkdim(state.state)

function combinerto(i1::Index, i2::Index, i3::Index)
    c = combiner(i1, i2)
    i = combinedind(c)
    replaceind(c, i, i3)
end

make_one_state(::TPure, system::System, i::Int, st::String) =
    state(system.pure_sites[i], st)

function make_one_state(::TMixed, system::System, i::Int, st::String)
    k = system.mixed_sites[i]
    try
        return state(k, st)
    catch   
        j = system.pure_sites[i]
        s = state(j, st)
        return s * dag(s') * combinerto(j', j, k)
    end
end

function make_state(type::Union{TPure, TMixed}, system::System, states::Vector{String})
    n = length(system)
    st = MPS(n)
    if n == 1
        st[1] = make_one_state(type, system, 1, states[1])
        return st
    end
    l = [ Index(1; tags="Link,l=$n") for n in 1:n ]
    st[1] =  make_one_state(type, system, 1, states[1]) * ITensor(1, l[1]) 
    for i in 2:n - 1
      st[i] = make_one_state(type, system, i, states[i]) * ITensor(1, l[i-1]) * ITensor(1, l[i])
    end
    st[n] = make_one_state(type, system, n, states[n]) * ITensor(1, l[n-1])
    return st
end

function State(type::Union{TPure, TMixed}, system::System, states::Vector{String})
    if length(system) ≠ length(states)
        error("incompatible sizes between system($(length(system))) and state($(length(states)))")
    end
    return State(type, system, make_state(type, system, states))
end

State(type::Union{TPure, TMixed}, system::System, state::String) =
    State(type, system, fill(state, length(system.pure_sites)))

State(state::State, st::MPS) =
    State(state.type, state.system, st)


"""
    mix(::State)

Transform a pure representation into a mixed representation
"""
function mix(state::State)
    if state.type == Mixed
        return state
    end
    n = length(state)
    system = state.system
    st = dense(state.state)
    v = Vector{ITensor}(undef, n)
    comb = ITensor(1)
    for (i, t) in enumerate(st)
        idx = system.pure_sites[i]
        midx = system.mixed_sites[i]
        mt = t * dag(t') * combinerto(idx', idx, midx)
        mt *= comb
        if (i < n)
            rlink = commonind(t, st[i+1])
            d = dim(rlink)
            comb = combinerto(rlink', rlink, Index(d*d, "Link,l=$i"))
        else
            comb = ITensor(1)
        end
        mt *= comb
        v[i] = mt
    end
    return State(Mixed, system, MPS(v))
end


"""
    truncate(::State; maxdim::Int, cutoff::Number)::State

Apply the truncation to the given state, in particular we get maxlinkdim(state)<=maxdim
"""
truncate(state::State; kwargs...) =
    State(state.type, state.system, truncate(state.state; kwargs...))

copy(state::State) = State(state, copy(state.state))

apply(mpo::MPO, state::State; kwargs...) =
    State(state, apply(mpo, state.state; kwargs...))