export State, maxlinkdim, Limits

"""
    struct PreObs

A data structure to hold preprocessing data for observable expectation computations.
Used by `State`
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

Either field may be given one value per sweep, as a vector, for a phase that sweeps:
`Limits(cutoff = 1e-14, maxdim = [2, 4, 8])` starts small and lets the state grow. A
schedule shorter than the number of sweeps is continued with its last value, which is
what ITensor does too.

# Examples

    Limits(cutoff = 1e-14, maxdim = 100)
    Limits(cutoff = 1e-14, maxdim = [10, 20, 50, 100])
"""
@kwdef struct Limits
    cutoff::Union{Float64, Vector{Float64}} = 0.
    maxdim::Union{Int, Vector{Int}} = typemax(Int)
end

"""
    sweep_value(x, sweep)
    sweep_limits(::Limits, sweep)

the value a per sweep schedule takes on the given sweep, and the `Limits` holding those
values. A plain value covers every sweep, and a schedule shorter than the number of sweeps
is continued with its last value, as ITensor does with its own.
"""
sweep_value(x, ::Int) = x
sweep_value(x::Vector, sweep::Int) = x[min(sweep, length(x))]
sweep_limits(l::Limits, sweep::Int) =
    Limits(sweep_value(l.cutoff, sweep), sweep_value(l.maxdim, sweep))

"""
    sweep_due(period, sweep)

whether something asked for every `period` sweeps is due on this one.

A period below one means never. That is what `0` was already taken to mean for `n_expand`,
`n_hermitianize` and `checkpoint_interval`, and the rule is extended to any value below
one so that a negative period is not silently read as `mod(sweep, -2)`, which is zero on
every second sweep. `mod(sweep, 0)` would raise a division by zero outright, which is what
`measures_period = 0` used to do.

The three sweep counters of the library go through this, and `checkpoint_due` applies the
same rule to the interval in seconds of the checkpointer.
"""
sweep_due(period::Int, sweep::Int) = period ≥ 1 && mod(sweep, period) == 0
  
"""
    struct State{R <: PM}
    State{R}(::System, states)
    State{R}(::Int, ::AbstractSite, state)
    State{R}(::Vector{<:AbstractSite}, state)
    State(::State, ::MPS)

represent the complete state of the simulated quantum system

# Type parameter

- `R` is `Pure` or `Mixed` and represent the type of representation used

# Fields

- `system::System`: system description
- `state::MPS`: system state
- `preobs::PreObs`: preprocessing data for computing observables

# Examples

    State{Pure}(system, "Up")
    State{Mixed}(system, ["Up", "Dn", "Up"])
    State{Mixed}(system, "FullyMixed")
    State{Pure}(system, [1, 0])
    State{Pure}(10, Qubit(), "Up")
    State{Mixed}([Qubit(), Boson(4), Fermion()], ["Up", "2", "Occ"])
    State(state, mps)        # a new state with the same system but a new mps

# Operations

states can be added, subtracted and multiplied by numbers

"""
struct State{R <: PM}
    system::System
    state::MPS
    preobs::PreObs
    State{R}(system::System, state::MPS) where {R <: PM} =
        new{R}(system, state, PreObs())
end

show(io::IO, s::State{R}) where R =
    print(io, "State{$R}($(s.system), (maxlinkdim = $(maxlinkdim(s.state))))")

"""
    length(::State)

return the number of sites in the state
"""
length(state::State) = length(state.system)

"""
    maxlinkdim(::State)

return the maximum link dimension in the state
"""
maxlinkdim(state::State) = maxlinkdim(state.state)


"""
    make_one_state(type::R, system::System, i::Int, st) where {R <: PM}

return the ITensor of the local state `st` at site `i` of the system
"""
make_one_state(type::R, system::System, i::Int, st) where {R <: PM} = 
    make_one_state(type, SysIndex{R}(system, i), state(system[i], st))

make_one_state(::Pure, i::Index, v::Vector) = ITensor(v, i)
make_one_state(::Pure, ::Index, ::Matrix) = error("cannot use a mixed local state to create a pure local state")
make_one_state(::Mixed, i::Index, v::Vector) = ITensor(v * v', i)
make_one_state(::Mixed, i::Index, m::Matrix) = ITensor(m, i)

"""
    make_state(type::R, system::System, states::Vector) where {R <: PM}

return the MPS of the local states `states` for the system
"""
function make_state(type::PM, system::System, states::Vector)
    n = length(system)
    st = MPS(n)
    if n == 1
        st[1] = make_one_state(type, system, 1, states[1])
    else
        l = [ Index(1; tags="Link,l=$k") for k in 1:n-1 ]
        st[1] =  make_one_state(type, system, 1, states[1]) * ITensor(1, l[1]) 
        for i in 2:n - 1
        st[i] = make_one_state(type, system, i, states[i]) * ITensor(1, l[i-1]) * ITensor(1, l[i])
        end
        st[n] = make_one_state(type, system, n, states[n]) * ITensor(1, l[n-1])
    end
    return st
end

function State{R}(system::System, states::Vector) where R
    if length(system) ≠ length(states)
        error("incompatible sizes between system ($(length(system))) and states ($(length(states)))")
    end
    return State{R}(system, make_state(R(), system, states))
end

State{R}(system::System, state) where R =
    State{R}(system, fill(state, length(system)))

State{R}(system::System, state::Union{Vector{<:Number}, Matrix}) where R =
    State{R}(system, fill(state, length(system)))

State(state::State{R}, st::MPS) where R =
    State{R}(state.system, st)

# a `System` draws ITensor indices of its own, so two states built on two systems cannot be
# contracted together even when they describe the very same sites. This is what puts one on
# the system of the other, and what `inner` and the fidelities point at when they refuse a
# pair of states
"""
    State(::System, ::State)

the same state on the given system, whose sites must be the ones the state was built on.

This is the mirror of `State(state, mps)`: that one keeps the system and takes a new mps,
this one keeps the mps and takes a new system. It is what makes a state read from disk, or
built before a run, comparable with the state of that run, since `inner` and the fidelities
require their two arguments to share a system.

# Examples

    ref = State(sim.state.system, load_state("ground.h5", "gs"))
"""
function State(system::System, st::State{R}) where R
    if system.sites ≠ st.system.sites
        error("cannot put a state on a system of other sites: the state was built on " *
              "$(st.system.sites) and the system given has $(system.sites)")
    end
    return State{R}(system, replace_siteinds(st.state, SysIndex{R}(system, 1:length(system))))
end

State{R}(size::Int, site::AbstractSite, state) where R =
    State{R}(System(size, site), state)

State{R}(sites::Vector{<:AbstractSite}, state) where R =
    State{R}(System(sites), state)

(a::Number * b::State) = State(b, a * b.state)
(a::State * b::Number) = b * a
(a::State / b::Number) = inv(b) * a
(-a::State) = -1 * a

+(a::State{R}, b::State{R}; limits::Limits=Limits()) where R =
    State(a, +(a.state, b.state; limits.cutoff, limits.maxdim))
-(a::State{R}, b::State{R}; limits::Limits=Limits()) where R =
    State(a, -(a.state, b.state; limits.cutoff, limits.maxdim))

"""
    mix(::State)

transform a pure representation into a mixed representation
"""
mix(state::State{Mixed}) = state

function mix(state::State{Pure})
    n = length(state)
    system = state.system
    st = dense(state.state)
    v = Vector{ITensor}(undef, n)
    comb = ITensor(1)
    for (i, t) in enumerate(st)
        idx = SysIndex{Pure}(system, i)
        midx = SysIndex{Mixed}(system, i)
        mt = t * dag(t') * combinerto(midx, idx, idx')
        mt *= comb
        if i < n
            rlink = commonind(t, st[i+1])
            d = dim(rlink)
            comb = combinerto(Index(d*d, "Link,l=$i"), rlink, rlink')
        else
            comb = ITensor(1)
        end
        mt *= comb
        v[i] = mt
    end
    return State{Mixed}(system, MPS(v))
end


"""
    truncate(::State; limits::Limits)

apply the truncations to the given state
"""
truncate(state::State{R}; limits::Limits) where R =
    State{R}(state.system, truncate(state.state; limits.cutoff, limits.maxdim))
