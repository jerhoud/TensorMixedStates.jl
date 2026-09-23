export State, mix, maxlinkdim, Limits

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
    qlinks::Vector{Index}
end
PreObs() = PreObs([], [], [], [], [])

"""
A type to hold MPS limits

# Fields
- `cutoff`: the cutoff under which singular values are neglected
- `maxdim`: the maximum bond dimension
- `mindim`: the minimum bond dimension, `0` meaning no minimum

Any field may be given one value per sweep, as a vector, for a phase that sweeps:
`Limits(cutoff = 1e-14, maxdim = [2, 4, 8])` starts small and lets the state grow. A
schedule shorter than the number of sweeps is continued with its last value, which is
what ITensor does too.

# Examples

    Limits(cutoff = 1e-14, maxdim = 100)
    Limits(cutoff = 1e-14, maxdim = [10, 20, 50, 100])
    Limits(cutoff = 1e-14, maxdim = 100, mindim = 10)
"""
@kwdef struct Limits
    cutoff::Union{Float64, Vector{Float64}} = 0.
    maxdim::Union{Int, Vector{Int}} = typemax(Int)
    mindim::Union{Int, Vector{Int}} = 0
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
    Limits(sweep_value(l.cutoff, sweep), sweep_value(l.maxdim, sweep),
           sweep_value(l.mindim, sweep))

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
    make_one_state(type, SysIndex{Pure}(system, i), SysIndex{Mixed}(system, i),
                   state(system[i], st), st, system[i])

make_one_state(::Pure, i::Index, ::Index, v::Vector, what, site::AbstractSite) =
    charged_state(() -> ITensor(v, i), i, what, site)
make_one_state(::Pure, ::Index, ::Index, ::Matrix, _, _) =
    error("cannot use a mixed local state to create a pure local state")
make_one_state(::Mixed, i::Index, k::Index, v::Vector, what, site::AbstractSite) =
    make_one_state(Mixed(), i, k, v * v', what, site)
# the density matrix is laid on the ket and the bra and only then gathered, never written
# straight onto the mixed index: combining charged indices merges and sorts their sectors,
# so the flat order of the mixed basis is not the order of the matrix
function make_one_state(::Mixed, i::Index, k::Index, m::Matrix, what, site::AbstractSite)
    b, c = mixer(i, k, site)
    return charged_state(() -> op_on_sites(m, [i], [dag(b')]), i, what, site) * c
end

"""
    state_links(ts)

the link indices of the product state made of the tensors `ts`.

Each one carries the charge of everything to its left, so that the flux of the whole state
is its sector. They are built from the right and daggered, which is the arrangement
ITensorMPS uses for its own product states and what makes the pieces contract. Without
charges they are the indices of dimension one they always were.
"""
function state_links(ts::Vector{ITensor})
    n = length(ts)
    if !hasqns(ts[1])
        return [ Index(1; tags = "Link,l=$k") for k in 1:n-1 ]
    end
    l = Vector{Index}(undef, n - 1)
    q = sum(flux, ts[1:n-1])
    for k in n-1:-1:1
        l[k] = dag(Index(q => 1; tags = "Link,l=$k"))
        q -= flux(ts[k])
    end
    return l
end

"""
    make_state(type::R, system::System, states::Vector) where {R <: PM}

return the MPS of the local states `states` for the system
"""
function make_state(type::PM, system::System, states::Vector)
    n = length(system)
    ts = [ make_one_state(type, system, i, states[i]) for i in 1:n ]
    st = MPS(n)
    if n == 1
        st[1] = ts[1]
    else
        l = state_links(ts)
        # `onehot` rather than `ITensor(1, l)`: the latter asks for a tensor of zero flux,
        # which a link carrying a charge has no block for
        st[1] = ts[1] * onehot(l[1] => 1)
        for i in 2:n-1
            st[i] = ts[i] * onehot(dag(l[i-1]) => 1) * onehot(l[i] => 1)
        end
        st[n] = ts[n] * onehot(dag(l[n-1]) => 1)
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
    State(a, +(a.state, b.state; limits.cutoff, limits.maxdim, limits.mindim))
-(a::State{R}, b::State{R}; limits::Limits=Limits()) where R =
    State(a, -(a.state, b.state; limits.cutoff, limits.maxdim, limits.mindim))

"""
    mix(::State)

transform a pure representation into a mixed representation
"""
mix(state::State{Mixed}) = state

function mix(state::State{Pure})
    n = length(state)
    system = state.system
    # densified so that a tensor left diagonal by a decomposition becomes an ordinary one,
    # which a charged state must not be: `dense` would throw its sectors away
    st = hasqns(state.state) ? state.state : dense(state.state)
    # the bra is a starred copy of the whole tensor, links included: a link carries the same
    # charges as the sites, and starring one half of a tensor while leaving the other alone
    # would have the two count in different ways. One copy per index, kept here, so that the
    # link starred on the right of a site is the same index on the left of the next
    starred = relabeller(i -> star(i, strong_names(system)))
    bra(t) = relabel(t, starred)
    v = Vector{ITensor}(undef, n)
    left = ITensor(1)
    for (i, t) in enumerate(st)
        idx = SysIndex{Pure}(system, i)
        midx = SysIndex{Mixed}(system, i)
        mt = t * dag(bra(t')) * combinerto(midx, idx, dag(starred(idx'))) * left
        if i < n
            rlink = commonind(t, st[i+1])
            # the combined link is taken from the combiner rather than named in advance:
            # combining charged indices merges and sorts their sectors, and only the
            # combiner knows which ones come out and in what order
            right = combiner(rlink, dag(starred(rlink')); tags = "Link,l=$i")
            mt *= right
            # the two ends of a link point in opposite directions, so the site on its right
            # gets the daggered combiner. Without charges this is the same tensor
            left = dag(right)
        end
        v[i] = mt
    end
    return State{Mixed}(system, MPS(v))
end


"""
    weaken(::State)
    weaken(::State, target)

the same state on a system conserving less, `target` naming what it must still conserve in
the vocabulary `conserve` takes. Without a target every strong quantity is asked for weakly,
or, when none is strong, every quantity is dropped: repeating it walks strong, then weak,
then nothing, and stops there.

The levels hold the same physics and differ in how the tensors are cut into blocks. Keeping
the charge of the ket apart from that of the bra cuts them finest and confines the state to a
single sector; keeping only the difference lets it spread over sectors and makes the trace a
product of one vector per site again; keeping nothing gives plain tensors.

Weakening is a step of a simulation in its own right: a phase may evolve under a strong
symmetry, which every dissipator commuting with the charge allows, and the next one continue
under a weak one, where a jump that moves the charge becomes possible. There is no way back,
the finer blocks not being recoverable from the coarser ones, and a target asking for more
than the state has is refused.

The system it builds is a new one, so two states weakened apart live on two systems and have
to be put on one another before `inner` will compare them.

A `Simulation` is weakened the same way, through its state.

# Examples

    weaken(state)                          # one rung down
    weaken(state, (strong(Ntot), 2Sz))     # exactly these
    weaken(state, ())                      # no charges at all
    weaken(state, symmetries(system))      # the identity
"""
function weaken(state::State{R}, target::Conserved) where R
    system = state.system
    source = symmetries(system)
    check_target(source, target, "this system")
    if target.names == source.names
        return state
    end
    weak = weaken(system, target)
    n = length(state)
    collapse, drop = transitions(source, target)
    if R === Pure
        # the pure index keeps one block per basis state, in the order of the basis, so its
        # flat order survives both the relabelling and the densifying and the tensors only
        # have to be put on the indices of the new system
        st = is_charged(weak) ?
            MPS([ relabel(t, relabeller(i -> weak_index(i, collapse, drop))) for t in state.state ]) :
            dense(state.state)
        return State{Pure}(weak, replace_siteinds(st, SysIndex{Pure}(weak, 1:n)))
    end
    if !is_charged(weak)
        # an ITensor holds either charged indices or plain ones, so the last rung densifies
        # first and then permutes, the two orders having nothing in common
        return State{Mixed}(weak,
            MPS([ dense(state.state[i]) * dense_map(system, weak, i) for i in 1:n ]))
    end
    relab = relabeller(i -> weak_index(i, collapse, drop))
    return State{Mixed}(weak,
        MPS([ relabel(state.state[i], relab) * weak_map(system, weak, i, relab)
              for i in 1:n ]))
end

weaken(state::State{R}, spec) where R = weaken(state, Conserved(spec_names(spec)))

weaken(state::State) = weaken(state, one_step_down(symmetries(state.system)))

"""
    truncate(::State; limits::Limits)

apply the truncations to the given state
"""
truncate(state::State{R}; limits::Limits) where R =
    State{R}(state.system, truncate(state.state; limits.cutoff, limits.maxdim, limits.mindim))
