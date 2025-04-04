export State, maxlinkdim, Limits

"""
    struct PreObs
    PreObs(::State)

A data structure to hold preprocessing data for observable expectation computations.
"""
struct PreObs
    loc::Vector{ITensor}
    left::Vector{ITensor}
    right::Vector{ITensor}
    trace::Vector{Number}
end
PreObs() = PreObs([], [], [], [])

"""
A type to hold MPS limits

# Fields
- `cutoff`: the cutoff under which singular values are neglected
- `maxdim`: the maximum bond dimension
"""
@kwdef struct Limits
    cutoff::Float64 = 0.
    maxdim::Int = typemax(Int)
end
  
"""
    struct State
    State(Pure()|Mixed(), system, states)
    State(::State, ::MPS)


Represent the complete state of the simulated quantum system

# Fields
- `type::Union{Pure, Mixed}`: pure or mixed representation
- `system::System`: system description
- `state::MPS`: system state

# Examples
    State(Pure(), system, "Up")
    State(Mixed(), system, ["Up", "Dn", "Up"])
    State(Mixed(), system, "FullyMixed")

# Operations
states can be added, substracted and multiplied by numbers

"""
struct State
    type::PM
    system::System
    state::MPS
    preops::PreOps
    State(type::PM, system::System, state::MPS) =
        new(type, system, state, PreObs())
end

"""
    length(::State)
    length(::Simulation)

Return the number of sites in the state
"""
length(state::State) = length(state.system)

"""
    maxlinkdim(::State)
    maxlinkdim(::Simulation)

Return the maximum link dimension in the state
"""
maxlinkdim(state::State) = maxlinkdim(state.state)

make_one_state(type::PM, system::System, i::Int, st) = 
    make_one_state(type, system[type, i], state(system[i], st))

make_one_state(::Pure, i::Index, v::Vector) = ITensor(v, i)
make_one_state(::Pure, ::Index, ::Matrix) = error("cannot use a mixed local state to create a pure local state")
make_one_state(::Mixed, i::Index, v::Vector) = ITensor(v * v', i)
make_one_state(::Mixed, i::Index, m::Matrix) = ITensor(m, i)


function make_state(type::PM, system::System, states::Vector)
    n = length(system)
    st = MPS(n)
    if n == 1
        st[1] = make_one_state(type, system, 1, states[1])
    else
        l = [ Index(1; tags="Link,l=$n") for n in 1:n-1 ]
        st[1] =  make_one_state(type, system, 1, states[1]) * ITensor(1, l[1]) 
        for i in 2:n - 1
        st[i] = make_one_state(type, system, i, states[i]) * ITensor(1, l[i-1]) * ITensor(1, l[i])
        end
        st[n] = make_one_state(type, system, n, states[n]) * ITensor(1, l[n-1])
    end
    return st
end

function State(type::PM, system::System, states::Vector)
    if length(system) ≠ length(states)
        error("incompatible sizes between system ($(length(system))) and states ($(length(states)))")
    end
    return State(type, system, make_state(type, system, states))
end

State(type::PM, system::System, state) =
    State(type, system, fill(state, length(system)))

State(type::PM, system::System, state::Vector{<:Number}) =
    State(type, system, fill(state, length(system)))

State(state::State, st::MPS) =
    State(state.type, state.system, st)

(a::Number * b::State) = State(b, a * b.state)
(a::State * b::Number) = b * a
(a::State / b::Number) = inv(b) * a
(-a::State) = -1 * a

+(a::State...; limits=Limits()) = State(a[1], +((i.state for i in a)...; limits.cutoff, limits.maxdim))
-(a::State, b::State; limits=Limits()) = State(a, -(a.state, b.state; limits.cutoff, limits.maxdim))

"""
    mix(::State)
    mix(::Simulation)

Transform a pure representation into a mixed representation
"""
mix(state::State) =
    mix(state.type, state)

mix(::Mixed, state::State) = state

function mix(::Pure, state::State)
    n = length(state)
    system = state.system
    st = dense(state.state)
    v = Vector{ITensor}(undef, n)
    comb = ITensor(1)
    for (i, t) in enumerate(st)
        idx = system[Pure(), i]
        midx = system[Mixed(), i]
        mt = t * dag(t') * combinerto(idx', idx, midx)
        mt *= comb
        if i < n
            rlink = commonind(t, st[i+1])
            d = dim(rlink)
            comb = combinerto(rlink', rlink, Index(d*d, "Link,l=$i"))
        else
            comb = ITensor(1)
        end
        mt *= comb
        v[i] = mt
    end
    return State(Mixed(), system, MPS(v))
end


"""
    truncate(::State; limits::Limits)
    truncate(::Simulation; limits::Limits)

Apply the truncation to the given state
"""
truncate(state::State; kwargs...) =
    State(state.type, state.system, truncate(state.state; limits.cutoff, limits.maxdim, limits.mindim))

apply(mpo::MPO, state::State; kwargs...) =
    State(state, apply(mpo, state.state; kwargs...))