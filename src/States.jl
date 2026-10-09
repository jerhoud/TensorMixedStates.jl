# States, MPS on a system in the pure or the mixed representation: product states built from
# local states, the truncation limits, mixing a pure state and weakening a state.

export State, mix, maxlinkdim, Limits

"""
    struct PreObs

the caches a `State` keeps for measurements, each filled on first use: the local tensors,
the environments from the left and from the right, the trace, and the weakened form of a
state conserving something strongly. `lock` is held while they are filled, so that several
threads may measure one state; it is reentrant, filling a cache filling others.
"""
struct PreObs
    loc::Vector{ITensor}
    left::Vector{ITensor}
    right::Vector{ITensor}
    trace::Vector{Number}
    weak::Vector{Any}
    lock::ReentrantLock
end
PreObs() = PreObs([], [], [], [], [], ReentrantLock())

"""
    Limits(; cutoff = 0., maxdim = typemax(Int), mindim = 1)

the truncation limits of an MPS.

# Fields
- `cutoff`: the largest truncation error allowed, the weight of the discarded singular values,
  the sum of their squares relative to that of all
- `maxdim`: the maximum bond dimension
- `mindim`: the minimum bond dimension; the default `1` means no minimum, and a smaller value
  is taken as `1`

For `tdvp`, `approx_W`, `dmrg`, `steady_state` and the phases built on them, any field may
be a vector, one value per sweep: `Limits(cutoff = 1e-14, maxdim = [2, 4, 8])` starts small
and lets the state grow. A vector shorter than the number of sweeps is continued with its last
value, as in ITensor. Elsewhere, as in `apply`, `truncate` or a sum of states, the fields are
single values.

# Examples

    Limits(cutoff = 1e-14, maxdim = 100)
    Limits(cutoff = 1e-14, maxdim = [10, 20, 50, 100])
    Limits(cutoff = 1e-14, maxdim = 100, mindim = 10)
"""
struct Limits
    cutoff::Union{Float64, Vector{Float64}}
    maxdim::Union{Int, Vector{Int}}
    mindim::Union{Int, Vector{Int}}
    # a bond has a dimension of one at least, which is what ITensors takes for no minimum:
    # told less, it truncates a spectrum of zeros past its first value, so that a gate of
    # several sites taking a state to zero, or the sum of two zero states, raised a BoundsError
    Limits(cutoff, maxdim, mindim) = new(cutoff, maxdim, at_least_one(mindim))
end

# a cutoff is a real number, and `cutoff = 0` has to be accepted as one: a field whose type
# is a union is not converted to, so `@kwdef` refused it
Limits(; cutoff = 0., maxdim = typemax(Int), mindim = 1) =
    Limits(float_cutoff(cutoff), maxdim, mindim)

# printed as the call that builds it, the fields left at their default omitted
function show(io::IO, l::Limits)
    default = Limits()
    given = [string(f, " = ", repr(getfield(l, f))) for f in fieldnames(Limits)
             if getfield(l, f) != getfield(default, f)]
    print(io, "Limits(", join(given, ", "), ")")
end

"""
    float_cutoff(x)

the cutoff `x`, a real number or a vector of them, as `Float64`.
"""
float_cutoff(x::Real) = Float64(x)
float_cutoff(x::AbstractVector{<:Real}) = Vector{Float64}(x)

"""
    at_least_one(m)

the minimum bond dimension `m`, or each value of a vector of them, raised to `1` if below.
"""
at_least_one(m::Int) = max(m, 1)
at_least_one(m::Vector{Int}) = max.(m, 1)

"""
    sweep_value(x, sweep)
    sweep_limits(::Limits, sweep)

the value of the per sweep schedule `x` on the given sweep, and the `Limits` holding those
values. A plain value covers every sweep; a vector shorter than the number of sweeps is
continued with its last value, as in ITensor.
"""
sweep_value(x, ::Int) = x
sweep_value(x::Vector, sweep::Int) = x[min(sweep, length(x))]
sweep_limits(l::Limits, sweep::Int) =
    Limits(sweep_value(l.cutoff, sweep), sweep_value(l.maxdim, sweep),
           sweep_value(l.mindim, sweep))

"""
    sweep_due(period, sweep)

whether something asked for every `period` sweeps is due on this one. A period below one
means never: `mod(sweep, 0)` would raise a division by zero, and a negative period would make
it due every `-period` sweeps. `checkpoint_due` applies the same rule.
"""
sweep_due(period::Int, sweep::Int) = period ≥ 1 && mod(sweep, period) == 0
  
"""
    struct State{R <: PM}
    State{R}(::System, states)
    State{R}(::Int, ::AbstractSite, state)
    State{R}(::Vector{<:AbstractSite}, state)
    State(::State, ::MPS)

the state of a quantum system, as an MPS: a wave function for `R = Pure`, a density matrix
for `R = Mixed`.

The local states are given as a vector, one per site, or as a single one for every site. A
local state is a name, the number of a basis state counted from 0, a vector of amplitudes, a
function of the site, the `AbstractSite` and not its position, giving one of those, or, in
mixed representation only, a density matrix.
A vector of numbers is always the amplitudes of one local state, given to every site: basis
numbers site by site are written `Any[0, 1, 0]` or `["0", "1", "0"]`.

# Fields

- `system::System`: the system
- `state::MPS`: the MPS
- `preobs::PreObs`: caches filled by measurements

# Examples

    State{Pure}(system, "Up")
    State{Mixed}(system, ["Up", "Dn", "Up"])
    State{Mixed}(system, "FullyMixed")
    State{Pure}(system, [1, 0])
    State{Pure}(10, Qubit(), "Up")
    State{Mixed}([Qubit(), Boson(4), Fermion()], ["Up", "2", "Occ"])
    State(state, mps)        # a new state with the same system but a new mps

# Operations

States can be added, subtracted, multiplied and divided by numbers. The two states of a sum
or a difference must be on the same system. It takes truncation limits as
`+(a, b; limits = Limits(maxdim = 100))`, and truncates nothing by default.
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

the number of sites of the state.
"""
length(state::State) = length(state.system)

"""
    maxlinkdim(::State)

the largest link dimension of the MPS of the state.
"""
maxlinkdim(state::State) = maxlinkdim(state.state)


"""
    make_one_state(type, system, i, st)

the ITensor of the local state `st` on site `i` of `system`, in the representation `type`,
`Pure()` or `Mixed()`. A vector gives a density matrix in mixed representation, and a density
matrix is refused in pure representation.
"""
make_one_state(type::R, system::System, i::Int, st) where {R <: PM} = 
    make_one_state(type, SysIndex{Pure}(system, i), SysIndex{Mixed}(system, i),
                   state(system[i], st), st, system[i])

make_one_state(::Pure, i::Index, ::Index, v::Vector, what, site::AbstractSite) =
    charged_state(v, [i], what, site)
make_one_state(::Pure, ::Index, ::Index, ::Matrix, _, _) =
    error("cannot use a mixed local state to create a pure local state")
make_one_state(::Mixed, i::Index, k::Index, v::Vector, what, site::AbstractSite) =
    make_one_state(Mixed(), i, k, v * v', what, site)
# the density matrix is laid on the ket and the bra and only then gathered, never written
# straight onto the mixed index: combining charged indices merges and sorts their sectors,
# so the flat order of the mixed basis is not the order of the matrix
function make_one_state(::Mixed, i::Index, k::Index, m::Matrix, what, site::AbstractSite)
    b, c = mixer(i, k, site)
    return charged_state(on_legs(m, [i], [dag(b')])..., what, site) * c
end

"""
    state_links(ts)

the link indices, of dimension one, of the product state of the tensors `ts`. With charges,
each carries the charge of the sites on its left, so that the flux of the state is its sector;
they are daggered, as ITensorMPS does for its own product states, so that the pieces contract.
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
    make_state(type, system, states)

the product state MPS of the local states `states`, one per site of `system`, in the
representation `type`, `Pure()` or `Mixed()`.
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

# the local state repeated on every site, in a list of type Any: `fill(3, n)` is a Vector{Int},
# which the method below took for the amplitudes of a single site
State{R}(system::System, state) where R =
    State{R}(system, Any[ state for _ in 1:length(system) ])

State{R}(system::System, state::Union{Vector{<:Number}, Matrix}) where R =
    State{R}(system, Any[ state for _ in 1:length(system) ])

State(state::State{R}, st::MPS) where R =
    State{R}(state.system, st)

"""
    State(::System, ::State)

the same state on the given system, whose sites must be those of the state.

Every `System` has ITensor indices of its own, so states built on two systems cannot be
contracted together even when their sites are the same, and `inner` and the fidelities refuse
them. This puts a state, read from disk or built before a run for instance, on the system of
another. It mirrors `State(state, mps)`, which keeps the system and takes a new MPS.

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

the state in mixed representation: the density matrix ``|\\psi\\rangle\\langle\\psi|`` of a
pure state ``|\\psi\\rangle``, or a mixed state unchanged.

# Examples

    mix(State{Pure}(10, Qubit(), "Up"))
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

the same state on a system conserving less, `target` naming what it must still conserve, as
`conserve` is given. Without a target, every strong quantity becomes weak or, when none is
strong, every quantity is dropped: repeated, it goes from strong to weak to nothing.

The levels hold the same physics and differ in how the tensors are cut into blocks. Strong
keeps the charges of ket and bra apart: the blocks are finest and the state lies in a single
sector. Weak keeps their difference: the state may spread over sectors. Nothing gives plain
tensors.

Weakening is a step of a simulation: a phase may evolve under a strong symmetry, which every
dissipator commuting with the charge allows, and the next one under a weak one, where a jump
moving the charge is possible.
There is no way back, since the finer blocks cannot be recovered, and a target asking for more
than the state has is refused.

The system built is a new one, so two states weakened separately must be put on one system,
see `State(::System, ::State)`, before `inner` compares them. A `Simulation` is weakened the
same way, through its state.

# Examples

    weaken(state)                          # one level down
    weaken(state, (strong(Ntot), 2Sz))     # exactly these
    weaken(state, ())                      # no charges at all
    weaken(state, symmetries(system))      # the identity
"""
function weaken(state::State{R}, target::Conserved) where R
    system = state.system
    source = symmetries(system)
    if target.names == source.names
        return state
    end
    weak = weaken(system, target)
    n = length(state)
    collapse, drop = transitions(source, target)
    relab = relabeller(i -> weak_index(i, collapse, drop))
    if R === Pure
        # the pure index keeps one block per basis state, in the order of the basis, so its
        # flat order survives both the relabelling and the densifying and the tensors only
        # have to be put on the indices of the new system. One relabeller for all the
        # tensors, or the two ends of a link would each get a new index of their own and the
        # state would come apart
        st = is_charged(weak) ? MPS([ relabel(t, relab) for t in state.state ]) :
                                dense(state.state)
        return State{Pure}(weak, replace_siteinds(st, SysIndex{Pure}(weak, 1:n)))
    end
    if !is_charged(weak)
        # an ITensor holds either charged indices or plain ones, so the last rung densifies
        # first and then permutes, the two orders having nothing in common
        return State{Mixed}(weak,
            MPS([ dense(state.state[i]) * dense_map(system, weak, i) for i in 1:n ]))
    end
    return State{Mixed}(weak,
        MPS([ relabel(state.state[i], relab) * weak_map(system, weak, i, relab)
              for i in 1:n ]))
end

weaken(state::State{R}, spec) where R = weaken(state, Conserved(spec_names(spec)))

weaken(state::State) = weaken(state, one_step_down(symmetries(state.system)))

"""
    truncate(::State; limits::Limits)

the state with its MPS truncated to `limits`.

# Examples

    truncate(state; limits = Limits(cutoff = 1e-10, maxdim = 50))
"""
truncate(state::State{R}; limits::Limits) where R =
    State{R}(state.system, truncate(state.state; limits.cutoff, limits.maxdim, limits.mindim))
