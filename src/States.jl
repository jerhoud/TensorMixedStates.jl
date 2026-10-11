# States, MPS on a system in the pure or the mixed representation: the truncation limits,
# product states built from local states, sums of states, mixing a pure state, weakening a
# state, and the norm, the trace and the inner product. The later layers go between the pure
# and the mixed representations with `combinerto`, `mixer`, `relabel`, `relabeller`,
# `element_maps` and `on_legs`, give limits to ITensors with `truncation`, and the measurements
# fill the cache `preobs` of a state with `cached!`, measuring on `weak_form`.

export State, Limits, mix, maxlinkdim, trace, norm, normalize, inner, dot

"""
    Limits(; cutoff = nothing, maxdim = nothing, mindim = nothing)

the truncation limits of an MPS, a field left `nothing` taking the default of the function the
limits are given to, which its docstring states.

- `cutoff`: the largest truncation error allowed, the sum of the squares of the discarded
  singular values relative to that of all, so that singular values up to `√cutoff` of the norm
  go: `1e-26` discards what rounding leaves, those below `1e-13` of the norm, `0` nothing
- `maxdim`: the maximum bond dimension
- `mindim`: the minimum bond dimension

For the evolutions and searches, any field may be a vector, one value per step or per sweep,
continued with its last value.

# Examples

    Limits(cutoff = 1e-14, maxdim = 100)
    Limits(maxdim = [10, 20, 50, 100])
"""
struct Limits
    cutoff::Union{Nothing, Float64, Vector{Float64}}
    maxdim::Union{Nothing, Int, Vector{Int}}
    mindim::Union{Nothing, Int, Vector{Int}}
    # at least one, which ITensors takes for no minimum: with less, it keeps nothing of a
    # spectrum of zeros
    Limits(cutoff, maxdim, mindim) = new(cutoff, maxdim, at_least_one(mindim))
end

"""
    at_least_one(m)

the minimum bond dimension `m`, or each value of a vector of them, raised to 1 if below
"""
at_least_one(m::Int) = max(m, 1)
at_least_one(m::Vector{Int}) = max.(m, 1)
at_least_one(::Nothing) = nothing

"""
    float_cutoff(x)
    int_dims(x, what)

the cutoff `x` as `Float64`, and the bond dimension `x` as `Int`, a float refused by a message
naming `what`, each possibly a vector, and `nothing` left as it is
"""
float_cutoff(x::Real) = Float64(x)
float_cutoff(x::AbstractVector{<:Real}) = Vector{Float64}(x)
float_cutoff(::Nothing) = nothing
int_dims(x::Integer, _) = Int(x)
int_dims(x::AbstractVector{<:Integer}, _) = Vector{Int}(x)
int_dims(::Nothing, _) = nothing
int_dims(x, what) = error("$what is an integer, or a vector of them, not $(repr(x))")

Limits(; cutoff = nothing, maxdim = nothing, mindim = nothing) =
    Limits(float_cutoff(cutoff), int_dims(maxdim, "maxdim"), int_dims(mindim, "mindim"))

"""
    truncation(limits; cutoff = rounding_tol^2, maxdim = typemax(Int), mindim = 1)

the keywords `cutoff`, `maxdim` and `mindim` of a truncation by ITensors, a field of `limits` left
`nothing` taking the default given, that of the function the limits are given to: by default,
the singular values below `rounding_tol` of the norm go, what rounding leaves
"""
truncation(l::Limits; cutoff = rounding_tol^2, maxdim = typemax(Int), mindim = 1) =
    (cutoff = something(l.cutoff, cutoff), maxdim = something(l.maxdim, maxdim),
     mindim = something(l.mindim, mindim))

(a::Limits == b::Limits) = (a.cutoff, a.maxdim, a.mindim) == (b.cutoff, b.maxdim, b.mindim)
hash(a::Limits, h::UInt) = hash((a.cutoff, a.maxdim, a.mindim), hash(Limits, h))

function show(io::IO, l::Limits)
    given = [ string(f, " = ", repr(getfield(l, f))) for f in fieldnames(Limits)
              if !isnothing(getfield(l, f)) ]
    print(io, "Limits(", join(given, ", "), ")")
end

"""
    abstract type AbstractState

the supertype of the states a simulation holds: `State`, and the states of a representation
an extension defines, which have a field `system`, the `System` of their physical sites.
"""
abstract type AbstractState end

"""
    PreObs

the caches a `State` keeps for measurements, filled on first use under `lock`: the local
tensors, the environments from the left and from the right, those from the left carrying the
Jordan-Wigner string, the trace, and the weakened form of a state conserving something strongly
"""
struct PreObs
    loc::Vector{ITensor}
    left::Vector{ITensor}
    right::Vector{ITensor}
    strung::Vector{ITensor}
    trace::Vector{Number}
    weak::Vector{Any}
    lock::ReentrantLock
end

PreObs() = PreObs([], [], [], [], [], [], ReentrantLock())

"""
    cached!(create!, v, state)

the cache `v` of `state`, filled by `create!(v, state)` on first use under the lock of
`PreObs`; once filled it never changes, and is read without the lock
"""
function cached!(create!, v, state)
    lock(state.preobs.lock) do
        if isempty(v)
            create!(v, state)
        end
    end
    return v
end

"""
    State{R}(system, states)
    State{R}(n::Int, site::AbstractSite, states)
    State{R}(sites::Vector{<:AbstractSite}, states)
    State(state, mps)

the state of a quantum system, as an MPS: a wave function for `R = Pure`, a density matrix for
`R = Mixed`. The local states are given as a vector, one per site, or as a single one for every
site: a name, `A => λ`, a vector of amplitudes, in mixed representation a density matrix, or a
function of the site giving one of those. A vector of numbers is the amplitudes of one local
state. `State(state, mps)` is a new state on the same system.

States can be added, subtracted, and multiplied or divided by numbers; a sum takes truncation
limits as `+(a, b; limits = Limits(maxdim = 100))`, by default `cutoff = 1e-26`, which discards
what rounding leaves.

# Examples

    State{Pure}(system, "Up")
    State{Mixed}(system, ["Up", "Dn", "Up"])
    State{Pure}(System(3, Qudit(4)), Level .=> [0, 3, 1])
    State{Mixed}([Qubit(), Boson(4), Fermion()], ["Up", "2", "Occ"])
"""
struct State{R <: PM} <: AbstractState
    system::System
    state::MPS
    preobs::PreObs
    State{R}(system::System, state::MPS) where {R <: PM} = new{R}(system, state, PreObs())
end

State(state::State{R}, st::MPS) where R = State{R}(state.system, st)

show(io::IO, s::State{R}) where R =
    print(io, "State{$R}($(s.system), (maxlinkdim = $(maxlinkdim(s.state))))")

"""
    length(::AbstractState)

the number of sites of the state.
"""
length(state::AbstractState) = length(state.system)

"""
    maxlinkdim(::State)

the largest link dimension of the MPS of the state.
"""
maxlinkdim(state::State) = maxlinkdim(state.state)

############### Between the pure and the mixed representations ###############

"""
    combinerto(i, j...)

a combiner gathering the indices `j` into the index `i`
"""
function combinerto(i::Index, j::Index...)
    c = combiner(j...; tags = "")
    return replaceind(c, combinedind(c), i)
end

"""
    mixer(j, k, site)

the pair `(b, c)`: `b` the index carrying the bra of the ket index `j` on `site`, and `c` the
combiner gathering `j` and `dag(b')` into the mixed index `k`. Everything crossing between the
pure and the mixed representations goes through it, so that the pairs are flattened in the
same order everywhere.
"""
function mixer(j::Index, k::Index, site::AbstractSite)
    b = bra_index(j, site)
    return b, combinerto(k, j, dag(b'))
end

"""
    ket_bras(system, i)

the function giving, for `(m, n)`, the element ``|m\\rangle\\langle n|`` of site `i`,
vectorised on its mixed index
"""
function ket_bras(system::System, i::Int)
    j = siteind(Pure, system, i)
    b, c = mixer(j, siteind(Mixed, system, i), system[i])
    return (m, n) -> onehot(j => m, dag(b') => n) * c
end

"""
    on_legs(m, outs, ins)

the matrix `m` of the combined space of several sites reshaped with one axis per index, and
those indices, `outs` then `ins`, each reversed: the array and the indices a tensor is built
from
"""
function on_legs(m::Matrix, outs, ins)
    idx = [ reverse(outs) ; reverse(ins) ]
    return (reshape(m, ntuple(k -> dim(idx[k]), length(idx))), idx)
end

"""
    relabel(t, f)

the tensor `t` with each of its indices passed through `f`, its storage untouched
"""
relabel(t::ITensor, f) = ITensors.setinds(t, map(f, inds(t)))

"""
    relabeller(f)

a function relabelling an index by `f`, memoised on its identity and prime level and keeping
its direction: the two ends of a link, one daggered, must be relabelled to the same index
"""
function relabeller(f)
    seen = Dict{Tuple{ITensors.IDType, Int}, Index}()
    return function (i::Index)
        x = get!(() -> f(i), seen, (id(i), plev(i)))
        return dir(x) == dir(i) ? x : dag(x)
    end
end

############### Product states ###############

"""
    charged_state(a, inds, what, site)

the tensor of the local state `what`, of array `a`, on the indices `inds` of `site`, refused
when it has no definite charge
"""
function charged_state(a::AbstractArray, inds, what, site::AbstractSite)
    if !has_definite_flux(a, inds)
        error("the state $(repr(what)) of site $(typeof(site)) spreads over several charges " *
              "of $(conserved(site)), so it has none of its own and cannot be used where " *
              "that is conserved")
    end
    return charged_itensor(a, inds)
end

"""
    make_one_state(R, system, i, state)

the ITensor of the local state `state` on site `i` of `system`, in the representation `R`
"""
function make_one_state(::Type{Pure}, system::System, i::Int, x)
    v = local_state(system[i], x)
    if !(v isa Vector)
        hint = x isa Pair ? ": A => λ gives a vector only for a diagonal operator and a " *
                            "simple eigenvalue" : ""
        error("$(repr(x)) is a density matrix on $(system[i]), which a pure state cannot take$hint")
    end
    return charged_state(v, [siteind(Pure, system, i)], x, system[i])
end

# laid on the ket and the bra, then combined: combining charged indices sorts their sectors,
# so the mixed basis is not in matrix order
function make_one_state(::Type{Mixed}, system::System, i::Int, x)
    v = local_state(system[i], x)
    m = v isa Vector ? v * v' : v
    j = siteind(Pure, system, i)
    b, c = mixer(j, siteind(Mixed, system, i), system[i])
    return charged_state(on_legs(m, [j], [dag(b')])..., x, system[i]) * c
end

"""
    state_links(ts)

the link indices, of dimension one, of the product state of the tensors `ts`, each carrying the
charge of the sites on its left, daggered so that the pieces contract
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
    make_state(R, system, states)

the product state MPS of the local states `states`, one per site of `system`
"""
function make_state(::Type{R}, system::System, states::Vector) where R
    n = length(system)
    ts = [ make_one_state(R, system, i, states[i]) for i in 1:n ]
    st = MPS(n)
    if n == 1
        st[1] = ts[1]
    else
        l = state_links(ts)
        # `onehot`: `ITensor(1, l)` asks for a zero flux, which a charged link has no block for
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
    return State{R}(system, make_state(R, system, states))
end

# a list of type Any: a vector of numbers is the amplitudes of a single local state
State{R}(system::System, state) where R = State{R}(system, Any[ state for _ in 1:length(system) ])

State{R}(system::System, state::Union{Vector{<:Number}, Matrix}) where R =
    State{R}(system, Any[ state for _ in 1:length(system) ])

State{R}(size::Int, site::AbstractSite, state) where R = State{R}(System(size, site), state)

State{R}(sites::Vector{<:AbstractSite}, state) where R = State{R}(System(sites), state)

############### Operations ###############

"""
    same_system(a, b)

refuse two states on different systems
"""
same_system(a::State, b::State) =
    if a.system ≠ b.system
        error("the two states are on different systems, $(a.system) and $(b.system)")
    end

(a::Number * b::State) = State(b, a * b.state)
(a::State * b::Number) = b * a
(a::State / b::Number) = inv(b) * a
-(a::State) = -1 * a

"""
    aligned_links(a, b)

the MPS `b` with each link pointing as that of `a`: a link reversed has its charges negated,
so that both tensors keep their flux, with no data moved. A direct sum asks for it.
"""
function aligned_links(a::MPS, b::MPS)
    b = copy(b)
    for k in 1:length(b)-1
        la, lb = linkind(a, k), linkind(b, k)
        if hasqns(lb) && dir(la) ≠ dir(lb)
            f = Index(id(lb), [ -q => d for (q, d) in space(lb) ], dir(dag(lb)), tags(lb), plev(lb))
            flip(i) = i == lb ? (dir(i) == dir(lb) ? f : dag(f)) : i
            b[k] = relabel(b[k], flip)
            b[k+1] = relabel(b[k+1], flip)
        end
    end
    return b
end

# a direct sum, exact, then truncated with respect to the result: the default algorithm of
# ITensorMPS cuts on the density matrices of the terms, which loses what is left of two close
# states once subtracted
function +(a::State{R}, b::State{R}; limits::Limits = Limits()) where R
    same_system(a, b)
    s = +(a.state, aligned_links(a.state, b.state); alg = "directsum")
    return State(a, truncate(s; truncation(limits)...))
end

-(a::State{R}, b::State{R}; limits::Limits = Limits()) where R = +(a, -b; limits)

"""
    truncate(::State; limits::Limits)

the state with its MPS truncated to `limits`, by default `cutoff = 1e-26`, which discards what
rounding leaves.

# Examples

    truncate(state; limits = Limits(cutoff = 1e-10, maxdim = 50))
"""
truncate(state::State; limits::Limits) =
    State(state, truncate(state.state; truncation(limits)...))

"""
    mix(::State)

the state in mixed representation: the density matrix ``|\\psi\\rangle\\langle\\psi|`` of a
pure state, or a mixed state unchanged.
"""
mix(state::State{Mixed}) = state

function mix(state::State{Pure})
    n = length(state)
    system = state.system
    # `dense` makes ordinary a tensor left diagonal by a decomposition, but would drop the
    # sectors of a charged state
    st = hasqns(state.state) ? state.state : dense(state.state)
    # the bra stars the whole tensor, links included, or links and sites would count charges
    # differently; one relabeller, so that a link starred on the right of a site is the same
    # index on the left of the next
    names = strong_names(system)
    starred = relabeller(i -> star(i, names))
    bra(t) = relabel(t, starred)
    v = Vector{ITensor}(undef, n)
    left = ITensor(1)
    for (i, t) in enumerate(st)
        idx = siteind(Pure, system, i)
        midx = siteind(Mixed, system, i)
        mt = t * dag(bra(t')) * combinerto(midx, idx, dag(starred(idx'))) * left
        if i < n
            rlink = commonind(t, st[i+1])
            # the combined link is the one of the combiner, which sorts charged sectors
            right = combiner(rlink, dag(starred(rlink')); tags = "Link,l=$i")
            mt *= right
            left = dag(right)
        end
        v[i] = mt
    end
    return State{Mixed}(system, MPS(v))
end

"""
    element_map(to, from, i, f)

the tensor pairing each element ``|m\\rangle\\langle n|`` of site `i` of the system `from`,
passed through `f`, with the same element of `to`
"""
function element_map(to::System, from::System, i::Int, f)
    d = dim(siteind(Pure, from, i))
    to_element, from_element = ket_bras(to, i), ket_bras(from, i)
    return sum(to_element(m, n) * dag(f(from_element(m, n))) for m in 1:d, n in 1:d)
end

"""
    element_maps(to, from, f)

`element_map` on every site, built once for the sites alike and carried onto the indices of
the others
"""
function element_maps(to::System, from::System, f)
    # the mixed index of site `i` in `to`, and the one `f` puts in place of that of `from`
    ends(i) = (siteind(Mixed, to, i), only(inds(f(onehot(siteind(Mixed, from, i) => 1)))))
    first_alike = Dict{Any, Int}()
    maps = ITensor[]
    for i in 1:length(to)
        key = (space(siteind(Pure, to, i)), space(siteind(Pure, from, i)),
               strong_names(to[i]), strong_names(from[i]))
        k = get!(first_alike, key, i)
        push!(maps, k == i ? element_map(to, from, i, f) : replaceinds(maps[k], ends(k), ends(i)))
    end
    return maps
end

"""
    pure_map(old, new)

the tensor carrying the pure index `old` of a site onto `new`, basis state by basis state
"""
pure_map(old::Index, new::Index) = sum(onehot(new => m) * dag(onehot(old => m)) for m in 1:dim(old))

"""
    weaken(::State)
    weaken(::State, target)

the same state on a system conserving less, `target` naming what it must still conserve, as
`conserve` is given. Without a target, every strong quantity becomes weak or, when none is
strong, every quantity is dropped. There is no way back.

# Examples

    weaken(state)                          # one level down
    weaken(state, (strong(Ntot), 2Sz))     # exactly these
    weaken(state, ())                      # no charges at all
"""
function weaken(state::State{R}, target::Conserved) where {R <: PM}
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
        if !is_charged(weak)
            return State{Pure}(weak, replace_siteinds(dense(state.state), siteinds(Pure, weak)))
        end
        # the relabelling keeps the basis order; only a site left conserving nothing needs a
        # map, the new system giving it a single block where the relabelling keeps one per
        # basis state
        return State{Pure}(weak, MPS(map(1:n) do i
            t = relabel(state.state[i], relab)
            old, new = relab(siteind(Pure, system, i)), siteind(Pure, weak, i)
            return space(old) == space(new) ? replaceind(t, old, new) : t * pure_map(old, new)
        end))
    end
    if !is_charged(weak)
        # densifying lays the blocks out in the order of the charges, which the maps permute
        # to the order of a plain combiner
        return State{Mixed}(weak, MPS([ dense(state.state[i]) * m
                                        for (i, m) in enumerate(element_maps(weak, system, dense)) ]))
    end
    # beyond renaming, the maps regroup the blocks, which the two mixed indices order
    # differently
    return State{Mixed}(weak, MPS([ relabel(state.state[i], relab) * m
                                    for (i, m) in enumerate(element_maps(weak, system, t -> relabel(t, relab))) ]))
end

weaken(state::State, spec) = weaken(state, Conserved(spec_names(spec)))

weaken(state::State) = weaken(state, one_step_down(symmetries(state.system)))

"""
    norm(::State)

the norm of the state: one for a normalized pure state, and on a mixed one the Hilbert-Schmidt
norm ``\\sqrt{\\mathrm{tr}(\\rho^\\dagger\\rho)}``, not the trace.
"""
norm(state::State) = norm(state.state)

"""
    weak_form(state)

the state measurements run on: `state` itself, or, for a mixed state conserving something
strongly, `weaken(state)`, computed once and kept with it. Under a strong symmetry an element
|m⟩⟨m| has no charge of its own, and weakening changes no expectation value.
"""
weak_form(state::State{Pure}) = state

weak_form(state::State{Mixed}) =
    if isempty(strong_names(state.system))
        state
    else
        only(cached!((w, s) -> push!(w, weaken(s)), state.preobs.weak, state))
    end

"""
    trace(::State)

the trace of the density matrix, or the squared norm of a pure state, computed once and kept
with the state.
"""
trace(state::State{Pure}) =
    only(cached!((t, s) -> push!(t, inner(s.state, s.state)), state.preobs.trace, state))

function trace(state::State{Mixed})
    t = cached!(state.preobs.trace, state) do t, state
        s = weak_form(state)
        x = ITensor(1.)
        for i in 1:length(s)
            e = ket_bras(s.system, i)
            x = x * s.state[i] * dag(sum(e(m, m) for m in 1:dim(s.system[i])))
        end
        push!(t, scalar(x))
    end
    return only(t)
end

"""
    normalize(::State)

the state rescaled to norm one on a pure representation, to a trace of real part one on a
mixed one.
"""
normalize(state::State{Pure}) = State(state, normalize(state.state))
normalize(state::State{Mixed}) = State(state, state.state / real(trace(state)))

"""
    inner(a::State, b::State)
    dot(a::State, b::State)

the inner product of two states of one system, ``\\langle a | b \\rangle`` on pure
representations, ``\\mathrm{tr}(a^\\dagger b)`` on mixed ones, not normalized.
"""
function inner(a::State{R}, b::State{R}) where R
    same_system(a, b)
    return dot(a.state, b.state)
end

inner(::State{Pure}, ::State{Mixed}) =
    error("no inner product between a pure and a mixed representation: mix the pure one")
inner(::State{Mixed}, ::State{Pure}) =
    error("no inner product between a pure and a mixed representation: mix the pure one")

dot(a::State, b::State) = inner(a, b)
