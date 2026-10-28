# States, MPS on a system in the pure or the mixed representation: the truncation limits,
# product states built from local states, sums of states, mixing a pure state, weakening a
# state, the norm, the trace and the inner product, and the product of two states and the
# permutation of the sites of one. The later layers go between the pure and the mixed
# representations with `combinerto`, `mixer`, `mixed_elements`, `relabel`, `relabeller`,
# `element_maps` and `on_legs`, give limits to ITensors with `truncation`, and the measurements
# fill the cache `preobs` of a state with `cached!`, measuring on `weak_form` or `plain_form`.

export State, Limits, mix, maxlinkdim, trace, tr, norm, normalize, inner, dot, permute_sites

"""
    check_limit(x, ok, what)

refuse a field of `Limits`, a value or a vector of them, holding a value for which `ok` is
false, or a vector holding none, which a truncation would read past its end, `what` saying what
it takes; `nothing` is the default of the function
"""
function check_limit(x, ok, what::String)
    if x isa AbstractVector && isempty(x)
        error("Limits takes $what, or one per step or sweep, not an empty vector")
    end
    for v in something(x, ())
        if !ok(v)
            error("Limits takes $what, not $v")
        end
    end
    return nothing
end

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
    function Limits(cutoff, maxdim, mindim)
        check_limit(cutoff, x -> x ≥ 0, "a cutoff of 0 or more")
        check_limit(maxdim, x -> x ≥ 1, "a maxdim of 1 or more")
        check_limit(mindim, x -> x ≥ 0, "a mindim of 0 or more")
        # at least one, which ITensors takes for no minimum: with less, it keeps nothing of a
        # spectrum of zeros
        return new(cutoff, maxdim, at_least_one(mindim))
    end
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
    truncation(limits; sweeps = false, cutoff = rounding_tol^2, maxdim = typemax(Int), mindim = 1)

the keywords `cutoff`, `maxdim` and `mindim` of a truncation by ITensors, a field of `limits` left
`nothing` taking the default given, that of the function the limits are given to: by default,
the singular values below `rounding_tol` of the norm go, what rounding leaves. A field holding a
vector, a value per sweep, is passed on only with `sweeps = true`, for the sweeps of dmrg
"""
function truncation(l::Limits; sweeps::Bool = false, cutoff = rounding_tol^2,
                    maxdim = typemax(Int), mindim = 1)
    if !sweeps && any(x -> x isa AbstractVector, (l.cutoff, l.maxdim, l.mindim))
        error("limits holding vectors, a value per step or sweep, are for the evolutions and " *
              "searches only: $l")
    end
    return (cutoff = something(l.cutoff, cutoff), maxdim = something(l.maxdim, maxdim),
            mindim = something(l.mindim, mindim))
end

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
Jordan-Wigner string, the trace, the weakened form of a state conserving something strongly,
the state conserving nothing, and the adjoint of a density matrix
"""
struct PreObs
    loc::Vector{ITensor}
    left::Vector{ITensor}
    right::Vector{ITensor}
    strung::Vector{ITensor}
    trace::Vector{Number}
    weak::Vector{Any}
    plain::Vector{Any}
    adjoint::Vector{Any}
    lock::ReentrantLock
end

PreObs() = PreObs([], [], [], [], [], [], [], [], ReentrantLock())

"""
    cached!(create!, v, state)

the cache `v` of `state`, filled by `create!(v, state)` on first use under the lock of
`PreObs`; once filled it never changes, and is read without the lock. Emptied when the filling
fails, by an interruption included, rather than left half filled for the next reading
"""
function cached!(create!, v, state)
    lock(state.preobs.lock) do
        if isempty(v)
            try
                create!(v, state)
            catch
                empty!(v)
                rethrow()
            end
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
    mixed_elements(system, i)

the element of the density matrix of site `i` of `system` at each position of its mixed index,
numbered as in the matrices of superoperators, the ket faster
"""
function mixed_elements(system::System, i::Int)
    k = siteind(Mixed, system, i)
    e = ket_bras(system, i)
    d = dim(system[i])
    return invperm(vec([ findfirst(!iszero, Array(e(m, n), k)) for m in 1:d, n in 1:d ]))
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
        error("incompatible sizes between system ($(length(system))) and states " *
              "($(length(states)))")
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
        error("the two states are on different systems, $(system_difference(a.system, b.system))")
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

# A workaround of ITensors, which cannot factorise a tensor of charged indices that is exactly
# zero, having no block: when a result vanishes, the difference of two equal states or gates
# annihilating a state, the SVD, the QR or the eigen of ITensorMPS fails, and so does any later
# step on that null state. The functions below run such a step and give the null state in that
# case, so that the code is written as it would be without the bug.

"""
    vanished(st)

whether the MPS `st` is the null state for having a tensor of zero, as a gate of one site or a
null result of `guarded` leaves it
"""
vanished(st::MPS) = any(t -> iszero(norm(t)), st)

"""
    guarded(f, vanishes, st)

the MPS `f(st)`, a step of ITensorMPS on `st`, or the null state when the result vanishes:
when `f` fails on charged indices, `st` itself if it is already null, which ITensors fails on,
and `0 * st` if `vanishes(st)` tells exactly that the result is zero; any other failure is
raised again. `f` is given a copy, which a failure may leave half changed.
"""
function guarded(f, vanishes, st::MPS)
    try
        return f(copy(st))
    catch
        if any(hasqns, st)
            if vanished(st)
                return st
            elseif vanishes(st)
                return 0 * st
            end
        end
        rethrow()
    end
end

"""
    vanishes_exactly(st)

whether the MPS `st` is exactly zero, read on its last site once orthogonalised there, where a
truncation starts
"""
vanishes_exactly(st::MPS) = iszero(norm(orthogonalize(st, length(st))[end]))

# a direct sum, exact, then truncated with respect to the result: the default algorithm of
# ITensorMPS cuts on the density matrices of the terms, which loses what is left of two close
# states once subtracted; two equal states give the null state, see `guarded`
function +(a::State{R}, b::State{R}; limits::Limits = Limits()) where R
    same_system(a, b)
    if length(a) == 1
        # nothing to truncate, and the direct sum needs links
        return State(a, MPS([ a.state[1] + b.state[1] ]))
    end
    s = +(a.state, aligned_links(a.state, b.state); alg = "directsum")
    return State(a, guarded(st -> truncate(st; truncation(limits)...), vanishes_exactly, s))
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
    State(state, guarded(st -> truncate(st; truncation(limits)...), vanishes_exactly, state.state))

"""
    mix(::State; limits = Limits())

the state in mixed representation: the density matrix ``|\\psi\\rangle\\langle\\psi|`` of a
pure state, of the square of its bond dimension, or a mixed state unchanged; truncated to
`limits`, see `truncate`, unless they are the default `Limits()`.

# Examples

    ρ = mix(ψ; limits = Limits(maxdim = 100))
"""
mix(state::State{Mixed}; limits::Limits = Limits()) =
    limits == Limits() ? state : truncate(state; limits)

function mix(state::State{Pure}; limits::Limits = Limits())
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
    return mix(State{Mixed}(system, MPS(v)); limits)
end

"""
    element_map(to, from, i, f; swap = false)

the tensor pairing each element ``|m\\rangle\\langle n|`` of site `i` of the system `from`,
passed through `f`, with the same element of `to`, or with ``|n\\rangle\\langle m|`` when
`swap` is true
"""
function element_map(to::System, from::System, i::Int, f; swap::Bool = false)
    d = dim(siteind(Pure, from, i))
    to_element, from_element = ket_bras(to, i), ket_bras(from, i)
    return sum(to_element((swap ? (n, m) : (m, n))...) * dag(f(from_element(m, n)))
               for m in 1:d, n in 1:d)
end

"""
    element_maps(to, from, f; swap = false)

`element_map` on every site, built once for the sites alike and carried onto the indices of
the others
"""
function element_maps(to::System, from::System, f; swap::Bool = false)
    # the mixed index of site `i` in `to`, and the one `f` puts in place of that of `from`
    ends(i) = (siteind(Mixed, to, i), only(inds(f(onehot(siteind(Mixed, from, i) => 1)))))
    first_alike = Dict{Any, Int}()
    maps = ITensor[]
    for i in 1:length(to)
        key = (space(siteind(Pure, to, i)), space(siteind(Pure, from, i)),
               strong_names(to[i]), strong_names(from[i]))
        k = get!(first_alike, key, i)
        push!(maps, k == i ? element_map(to, from, i, f; swap) :
                             replaceinds(maps[k], ends(k), ends(i)))
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
        maps = element_maps(weak, system, dense)
        return State{Mixed}(weak, MPS([ dense(state.state[i]) * m for (i, m) in enumerate(maps) ]))
    end
    # beyond renaming, the maps regroup the blocks, which the two mixed indices order
    # differently
    maps = element_maps(weak, system, t -> relabel(t, relab))
    return State{Mixed}(weak, MPS([ relabel(state.state[i], relab) * m
                                    for (i, m) in enumerate(maps) ]))
end

weaken(state::State, spec) = weaken(state, Conserved(spec_names(spec)))

weaken(state::State) = weaken(state, one_step_down(symmetries(state.system)))

symmetries(state::AbstractState) = symmetries(state.system)

"""
    common_symmetries(a, b)

what the systems `a` and `b` both conserve, a quantity strong where it is strong on both
"""
function common_symmetries(a::System, b::System)
    strong_on_b = Dict(symmetries(b).names)
    return Conserved([ (n, s && strong_on_b[n]) for (n, s) in symmetries(a).names
                       if haskey(strong_on_b, n) ])
end

"""
    on_common_system(a, b)

the states `a` and `b` of the same sites on what both conserve, see `common_symmetries`, each
weakened as needed; states of other sites refused
"""
function on_common_system(a::State, b::State)
    if a.system ≠ b.system
        target = common_symmetries(a.system, b.system)
        a, b = weaken(a, target), weaken(b, target)
        same_system(a, b)
    end
    return a, b
end

"""
    norm(::State)

the norm of the state: one for a normalized pure state, and on a mixed one the Hilbert-Schmidt
norm ``\\sqrt{\\mathrm{tr}(\\rho^\\dagger\\rho)}``, not the trace.
"""
norm(state::State) = norm(state.state)

# LinearAlgebra takes any argument for norm, normalize and dot, and failed on iterate for a state
# of an extension that defines none of them, rather than saying so
norm(state::AbstractState) = throw(MethodError(norm, (state,)))

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
    plain_form(state)

the state on its system conserving nothing, see `weaken`, computed once and kept: where an
operator changing the charges by no definite amount is measured
"""
plain_form(state::State) =
    if !is_charged(state.system)
        state
    else
        only(cached!((p, s) -> push!(p, weaken(s, ())), state.preobs.plain, state))
    end

"""
    trace_at(system, i)

the tensor taking the trace of a density matrix over site `i` of `system`, on its mixed index
"""
function trace_at(system::System, i::Int)
    e = ket_bras(system, i)
    return dag(sum(e(m, m) for m in 1:dim(system[i])))
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
            x = x * s.state[i] * trace_at(s.system, i)
        end
        push!(t, scalar(x))
    end
    return only(t)
end

"""
    tr(::AbstractState)

the trace of the state, `trace` under the name `LinearAlgebra` gives it
"""
tr(state::AbstractState) = trace(state)

"""
    nonzero_trace(state, what)
    nonzero_norm(state, what)

the real part of the trace, or the norm, of `state`, by which `what` divides, refused when it
is zero, as for the null state an operator annihilating a state leaves
"""
function nonzero_trace(state::State, what)
    t = real(trace(state))
    if iszero(t)
        error("$what divides by the trace of the state, which is zero, as for the null state")
    end
    return t
end

function nonzero_norm(state::State, what)
    n = norm(state.state)
    if iszero(n)
        error("$what divides by the norm of the state, which is zero, as for the null state")
    end
    return n
end

"""
    traced(state, kept)

the mixed `state` with its sites outside the increasing positions `kept` traced out, on the
system of the sites kept, its indices densified when these conserve nothing; no fermionic sign
is taken, see `partial_trace`. Nothing is divided by the trace nor truncated: a cut would leave
a trace off the one of `state` and negative eigenvalues.
"""
function traced(state::State{Mixed}, kept::AbstractVector{Int})
    system = System(state.system.sites[kept])
    loc(i) = state.state[i] * trace_at(state.system, i)
    # the sites on either side of those kept traced out from the ends inwards, and those
    # between two sites kept gathered before going into the first
    ts = [ state.state[k] for k in kept ]
    ts[1] = foldl((t, i) -> t * loc(i), 1:first(kept)-1; init = ITensor(1.)) * ts[1]
    for j in 1:length(kept)-1
        if kept[j + 1] > kept[j] + 1
            ts[j] *= prod(loc(i) for i in kept[j]+1:kept[j+1]-1)
        end
    end
    ts[end] *= foldr((i, t) -> loc(i) * t, last(kept)+1:length(state); init = ITensor(1.))
    mps = MPS(ts)
    if hasqns(ts[1]) && !is_charged(system)
        mps = dense(mps)
    end
    return State{Mixed}(system, replace_siteinds(mps, siteinds(Mixed, system)))
end

"""
    normalize(::State)

the state rescaled to norm one on a pure representation, to a trace of real part one on a
mixed one.
"""
function normalize(state::State{Pure})
    nonzero_trace(state, "normalize")
    return State(state, normalize(state.state))
end

normalize(state::State{Mixed}) = State(state, state.state / nonzero_trace(state, "normalize"))

normalize(state::AbstractState) = throw(MethodError(normalize, (state,)))

"""
    inner(a::State, b::State)
    dot(a::State, b::State)

the inner product of two states of the same sites, taken on what both conserve,
``\\langle a | b \\rangle`` on pure representations, ``\\mathrm{tr}(a^\\dagger b)`` on mixed
ones, not normalized.
"""
function inner(a::State{R}, b::State{R}) where R
    a, b = on_common_system(a, b)
    return dot(a.state, b.state)
end

inner(::State{Pure}, ::State{Mixed}) =
    error("no inner product between a pure and a mixed representation: mix the pure one")
inner(::State{Mixed}, ::State{Pure}) =
    error("no inner product between a pure and a mixed representation: mix the pure one")

dot(a::State, b::State) = inner(a, b)

dot(a::AbstractState, b::AbstractState) = throw(MethodError(dot, (a, b)))

############### Products and permutations ###############

"""
    lifted(R, st, from, to, offset)

the tensors of the MPS `st` of a state of `from`, a system conserving nothing, laid on the sites
`offset + 1` on of `to`, which conserves something: each index given the single sector of charge
zero a site conserving nothing takes there, see `site_index`, the elements of a mixed index in
the order of `to`, see `mixed_elements`
"""
function lifted(::Type{R}, st::MPS, from::System, to::System, offset::Int) where R
    n = length(st)
    links = [ commonind(st[k], st[k+1]) for k in 1:n-1 ]
    charged = [ Index(QN() => dim(l); tags = "Link,l=$(offset + k)")
                for (k, l) in enumerate(links) ]
    return map(1:n) do k
        left, right = max(1, k - 1):k-1, k:min(k, n - 1)
        old = [ siteind(R, from, k) ; links[left] ; links[right] ]
        new = [ siteind(R, to, offset + k) ; dag.(charged[left]) ; charged[right] ]
        a = Array(st[k], old...)
        if R === Mixed
            elements = invperm(mixed_elements(from, k))[mixed_elements(to, offset + k)]
            a = a[elements, ntuple(_ -> Colon(), ndims(a) - 1)...]
        end
        return charged_itensor(a, new)
    end
end

"""
    on_system(state, system, offset)

the tensors of the MPS of `state` laid on the sites `offset + 1` on of `system`, which holds its
sites there, its links renamed, so that the two copies of a state laid twice stay apart
"""
function on_system(state::State{R}, system::System, offset::Int) where R
    from = state.system
    if is_charged(system) && !is_charged(from)
        return lifted(R, state.state, from, system, offset)
    end
    st = state.state
    n = length(st)
    new = Dict{ITensors.IDType, Index}(id(siteind(R, from, k)) => siteind(R, system, offset + k)
                                       for k in 1:n)
    for k in 1:n-1
        l = commonind(st[k], st[k+1])
        new[id(l)] = settags(sim(l), "Link,l=$(offset + k)")
    end
    relab = relabeller(i -> new[id(i)])
    return [ relabel(t, relab) for t in st ]
end

"""
    a ⊗ b
    tensor(a, b)

the product of the states `a` and `b`, on the system `a.system ⊗ b.system`, the sites of `b`
following those of `a`: exact, its bond dimension the larger of theirs, its norm or its trace
the product of theirs. A pure state with a mixed one is mixed, see `mix`. On fermions, the
operators of `a` stand on the left of those of `b`: the fermionic product of two states of
definite parity.

# Examples

    ψ ⊗ State{Pure}(System(2, Qubit()), "Up")       # two qubits more
    left ⊗ right                                     # two halves made apart
"""
function ⊗(a::State{R}, b::State{R}) where R
    system = a.system ⊗ b.system
    na = length(a)
    ts = [ on_system(a, system, 0) ; on_system(b, system, na) ]
    # of dimension one, carrying the charge of the sites of `a`, as the links of a product state,
    # see `state_links`
    l = if hasqns(ts[1])
        q = sum(t -> something(flux(t), QN()), ts[1:na]; init = QN())
        dag(Index(q => 1; tags = "Link,l=$na"))
    else
        Index(1; tags = "Link,l=$na")
    end
    ts[na] *= onehot(l => 1)
    ts[na + 1] *= onehot(dag(l) => 1)
    return State{R}(system, MPS(ts))
end

⊗(a::State{Pure}, b::State{Mixed}) = mix(a) ⊗ b
⊗(a::State{Mixed}, b::State{Pure}) = a ⊗ mix(b)

"""
    element_parities(R, system, i)

the fermionic parity of each value of the index of site `i` of `system`, see
`fermion_parities`, on a density matrix those of the kets and those of the bras
"""
element_parities(::Type{Pure}, system::System, i::Int) = (collect(fermion_parities(system[i])),)

function element_parities(::Type{Mixed}, system::System, i::Int)
    p = fermion_parities(system[i])
    d = length(p)
    e = mixed_elements(system, i)
    return ([ p[mod1(x, d)] for x in e ], [ p[cld(x, d)] for x in e ])
end

"""
    exchange_sign(R, system, x, y, phi)

the tensor `phi` of the neighbouring sites `x` and `y` of `system` with the sign of their
exchange, ``(-1)^{p_a p_b}`` for the values `a` of `x` and `b` of `y` of parities ``p_a`` and
``p_b``, on the kets and on the bras of a density matrix
"""
function exchange_sign(::Type{R}, system::System, x::Int, y::Int, phi::ITensor) where R
    signs(i, p) = charged_itensor(diagm(Float64[ isodd(v) ? -1 : 1 for v in p ]), [i', dag(i)])
    for (u, v) in zip(element_parities(R, system, x), element_parities(R, system, y))
        if any(isodd, u) && any(isodd, v)
            ux, vy = signs(siteind(R, system, x), u), signs(siteind(R, system, y), v)
            a, b = noprime(phi * ux), noprime(phi * vy)
            # (-1)^(pa pb) = (1 + (-1)^pa + (-1)^pb - (-1)^(pa + pb)) / 2
            phi = (phi + a + b - noprime(a * vy)) / 2
        end
    end
    return phi
end

"""
    null_mps(R, system)

the MPS of the null state of `system`, its links of dimension one, its tensors without
elements, which the null state an exchange of sites cannot factorize gives instead
"""
function null_mps(::Type{R}, system::System) where R
    n = length(system)
    ls = [ is_charged(system) ? Index(QN() => 1; tags = "Link,l=$k") : Index(1; tags = "Link,l=$k")
           for k in 1:n-1 ]
    return MPS([ ITensor(Float64, siteind(R, system, k), dag.(ls[max(1, k - 1):k-1])...,
                         ls[k:min(k, n - 1)]...) for k in 1:n ])
end

"""
    permute_sites(state, p; limits = Limits())

the state with its sites moved, site `k` to position `p[k]`, `p` a permutation of the
positions: an operator measures on it what `map_sites(k -> p[k], op)` measures on `state`. The
sites are moved by exchanges of neighbours, each truncated by `limits`, by default
`cutoff = 1e-26`, as many as the pairs of sites `p` puts in the other order. The bond dimension
may grow: two states interleaved, `a ⊗ b` with the sites of `b` moved between those of `a`,
have the product of their bond dimensions. On fermions, two odd states exchanged give a sign,
as `Swap`.

# Examples

    permute_sites(ψ, [1, 3, 5, 2, 4, 6])        # the halves of a chain of 6 sites interleaved
    permute_sites(ψ, n:-1:1)                    # a chain of n sites reversed
"""
function permute_sites(state::State{R}, p::AbstractVector{<:Integer};
                       limits::Limits = Limits()) where R
    n = length(state)
    if sort(p) ≠ 1:n
        error("permute_sites needs a permutation of the positions 1 to $n, not $p")
    end
    system = state.system
    new = System(system.sites[invperm(p)])
    st = copy(state.state)
    if iszero(norm(st))
        return State{R}(new, null_mps(R, new))
    end
    # `at[j]`, the site now at position `j`, sorted on the positions they go to by exchanges of
    # neighbours, one for each pair in the wrong order
    at = collect(1:n)
    trunc = truncation(limits)
    exchanged = true
    while exchanged
        exchanged = false
        for b in 1:n-1
            if p[at[b]] > p[at[b+1]]
                orthogonalize!(st, b)
                phi = exchange_sign(R, system, at[b], at[b+1], st[b] * st[b+1])
                # as swapbondsites, which takes sites of different spaces, our sign laid before:
                # ITensors knows no fermion on our indices
                setindex!(st, phi, b:b+1; perm = [2, 1], trunc...)
                at[b], at[b+1] = at[b+1], at[b]
                exchanged = true
            end
        end
    end
    sites = Dict(id(siteind(R, system, at[j])) => siteind(R, new, j) for j in 1:n)
    relab = relabeller(i -> get(sites, id(i), i))
    return State{R}(new, MPS([ relabel(t, relab) for t in st ]))
end
