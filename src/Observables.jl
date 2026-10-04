# What is computed on a state: traces, norms, inner products and fidelities, expectation values
# and correlations, entropies, partial traces and sampling.

export trace, trace2, norm, normalize, hermitianize, hermiticity, renyi2
export inner, dot, fidelity, hs_fidelity
export expect, expect1, expect2
export entanglement_entropy, entanglement_by_sector, partial_trace, mutual_info_renyi2, sample, variance

"""
    weak_form(state)

the state measurements run on: `state` itself, or, for a mixed state conserving something
strongly, `weaken(state)`, computed once and kept with it.

Under a strong symmetry the ket and the bra carry charges of their own, so the trace is not a
product of one vector per site; under the weak one it is. Weakening is exact, so no
expectation value changes, and an operator measurable on the strong state is measurable on
the weak one.
"""
weak_form(state::State{Pure}) = state

function weak_form(state::State{Mixed})
    if isempty(strong_names(state.system))
        return state
    end
    w = state.preobs.weak
    if isempty(w)
        push!(w, weaken(state))
    end
    return only(w)
end

"""
    strong_measured(i)

raise the error for site `i` of a strongly conserving state measured without going through
`weak_form`, which is a bug: only `weak_form` gives such a state a trace with one vector per
site.
"""
strong_measured(i::Int) =
    error("bug: site $i of a strongly conserving state measured without going through weak_form")

"""
    on_trace(state, t, i)

the tensor `t` of a pure operator on site `i` carried onto the mixed index of the state
there, so that contracting it with the density matrix applies the operator and traces the
site out.
"""
function on_trace(state::State{Mixed}, t::ITensor, i::Int)
    s = state.system
    j = SysIndex{Pure}(s, i)
    k = SysIndex{Mixed}(s, i)
    b, c = mixer(j, k, s[i])
    if b !== j
        strong_measured(i)
    end
    # daggered so that the result meets the `k` of the state and not another copy of it
    return t * dag(c)
end

"""
    identity_at(state, i)

the tensor measuring the identity on site `i`, see `tensor_obs`: the identity of the site on a
pure state; on a mixed one, the trace over the site, which, contracted with the tensor of the
density matrix there, leaves only the links, the state conserving nothing strongly, see
`weak_form`. A placed identity has no site, while `expect1` and `expect2`, which close an
environment on a given site, need its tensor there: the site comes from them.
"""
function identity_at(state::State{Pure}, i::Int)
    j = SysIndex{Pure}(state.system, i)
    return delta(dag(j), j')
end

function identity_at(state::State{Mixed}, i::Int)
    j = SysIndex{Pure}(state.system, i)
    return on_trace(state, denseblocks(delta(dag(j), j')), i)
end

"""
    tensor_obs(state, o)

the tensor measuring the one site operator `o`, possibly times a coefficient, on `state`: the
tensor of the operator on a pure representation; on a mixed one, the tensor that, contracted
with the density matrix on that site, applies `o` and traces the site out, see `on_trace`.
"""
tensor_obs(state::State{Pure}, ind::AtIndex{Pure, 1}) =
    tensor(state.system, ind)

tensor_obs(state::State{Mixed}, ind::AtIndex{Pure, 1}) =
    on_trace(state, tensor(state.system, ind), only(ind.index))

# `(c * A)(i)` keeps its coefficient outside the AtIndex, so it has to be taken off here:
# everything downstream of `tensor_obs` works on the one site tensor alone
tensor_obs(state::State, a::ScalarOp{Pure, Indexed, 1}) =
    a.coef * tensor_obs(state, a.arg)

"""
    obs_at(state, op, i)

the tensor measuring the one site operator `op` on site `i`, see `tensor_obs`, a multiple of the
identity included, see `identity_at`: `(2X)^2`, simplified to `4Id`, went to `tensor_obs`, which
an identity, placed on no site, has no method of.
"""
obs_at(state::State, op::SimpleOp, i::Int) =
    scalararg(op) isa IdentityOp ? scalarcoef(op) * identity_at(state, i) : tensor_obs(state, op(i))

"""
    tensor_dag(state, i)

the tensor on site `i` of the adjoint of the density matrix, for a state conserving nothing
strongly, see `dag`.
"""
function tensor_dag(state::State, i::Int)
    s = state.system
    j = SysIndex{Pure}(s, i)
    k = SysIndex{Mixed}(s, i)
    c = last(mixer(j, k, s[i]))
    # the conjugate is spread back over the two indices, the ket and the bra are exchanged,
    # and the pair is gathered again: that transposition is what turns a conjugate into an
    # adjoint. The exchange is a renaming rather than a second combiner in the other order,
    # which would ask for a mixed index of the opposite charge
    return replaceinds(dag(state.state[i]) * c, (dag(j), j'), (dag(j'), j)) * c
end

"""
    cached!(create!, v, state)

the cache `v` of `state`, filled by `create!(v, state)` on first use
"""
function cached!(create!, v, state)
    if isempty(v)
        create!(v, state)
    end
    return v
end

"""
    create_loc!(l, state)

fill `l` with what `get_loc` returns for every site. A pure state has none: asking for it is
a bug.
"""
create_loc!(_, ::State{Pure}) = error("bug: get_loc on pure states")
function create_loc!(l, state::State{Mixed})
    n = length(state)
    resize!(l, n)
    for i in 1:n
        l[i] = identity_at(state, i) * state.state[i]
    end
    return l
end

"""
    get_loc(state, i)

the tensor of the density matrix on site `i` with that site traced out, see `identity_at`.
Those of all sites are computed on first use and kept with the state. Mixed states only.
"""
get_loc(state::State, i::Int) = cached!(create_loc!, state.preobs.loc, state)[i]

"""
    create_right!(r, state)

fill `r` with what `get_right` returns for every site. On a pure state, the sites from
`rightlim` on, being right orthogonal, are not contracted: their environment is the identity.
"""
function create_right!(r, state::State{Pure})
    st = state.state
    n = length(state)
    rl = ITensorMPS.rightlim(st) - 1
    resize!(r, n)
    r[n] = dag(st[n]')
    for i in n-1:-1:1
        rlink = commonind(st[i], st[i+1])
        v = if i >= rl
            delta(dag(rlink), rlink')
        else
            r[i+1] * identity_at(state, i+1) * st[i+1]
        end
        r[i] = v * dag(st[i]')
    end
    return r
end

function create_right!(r, state::State{Mixed})
    n = length(state)
    resize!(r, n)
    t = r[n] = ITensor(1.)
    for i in n-1:-1:1
        t = r[i] = get_loc(state, i+1) * t
    end
    return r
end

"""
    get_right(state, i)

the environment on the right of site `i`: the sites after `i` contracted, traced out on a
mixed state, ket with bra on a pure one, where the bra of site `i` is included as well.
Those of all sites are computed on first use and kept with the state.
"""
get_right(state::State, i::Int) = cached!(create_right!, state.preobs.right, state)[i]

"""
    create_trace!(t, state)

store the trace of `state` in `t[1]`, see `trace`.
"""
function create_trace!(t, state::State{Pure})
    resize!(t, 1)
    t[1] = scalar(get_right(state, 1) * identity_at(state, 1) * state.state[1])
    return t
end

function create_trace!(t, state::State{Mixed})
    resize!(t, 1)
    t[1] = scalar(get_loc(state, 1) * get_right(state, 1))
    return t
end

"""
    trace(::State)

the trace of the density matrix, which should be one. On a pure representation it is the
squared norm of the state, the trace of ``|\\psi\\rangle\\langle\\psi|``. It is computed
once and kept with the state.

It is a complex number whenever the tensors of the state are, after an evolution for
instance: its imaginary part is then numerical error, and `Trace` measures the real part.
"""
function trace(state::State)
    state = weak_form(state)
    return cached!(create_trace!, state.preobs.trace, state)[1]
end

"""
    extend_left!(l, state, i)

extend `l`, which holds what `get_left` returns for its first sites, the first one at least,
up to site `i`. On a pure state, the sites up to `leftlim`, being left orthogonal, are not
contracted: their environment is the identity.
"""
function extend_left!(l, state::State{Pure}, i::Int)
    st = state.state
    ll = ITensorMPS.leftlim(st) + 1
    j = length(l)
    resize!(l, i)
    for k in j+1:i
        llink = commonind(st[k-1], st[k])
        v = if k <= ll
            # what the branch below leaves: the link of the ket as st[k-1] holds it, and that
            # of the bra daggered and primed. The other way round, which only charges tell
            # apart, it did not contract with st[k]
            delta(llink, dag(llink)')
        else
            l[k-1] * identity_at(state, k-1) * dag(st[k-1]')
        end
        l[k] = v * st[k]
    end
    return l
end

function extend_left!(l, state::State{Mixed}, i::Int)
    st = state.state
    j = length(l)
    resize!(l, i)
    for k in j+1:i
        l[k] = l[k-1] * identity_at(state, k-1) * st[k]
    end
    return l
end

"""
    get_left(state, i)

the environment on the left of site `i`: the sites before `i` contracted, traced out on a mixed
state, ket with bra on a pure one, followed by the tensor of the state on `i`. It is not
divided by the trace, which the measurements built on it divide by when they normalize, so
that the same environments serve both. Computed up to `i` on demand and kept with the state.
"""
function get_left(state::State, i::Int)
    l = state.preobs.left
    if isempty(l)
        push!(l, state.state[1])
    end
    if length(l) < i
        extend_left!(l, state, i)
    end
    return l[i]
end

"""
    trace2(::State)

the purity ``\\mathrm{tr}(\\rho^2)`` of the density matrix normalised to trace one: 1 for a
pure state, less for a mixed one. On a pure representation it is 1 without computation.
"""
trace2(::State{Pure}) = 1.
trace2(state::State{Mixed}) = (norm(state.state) / real(trace(state))) ^ 2

"""
    norm(::State)

the norm of the state, which should be one on a pure representation. On a mixed one it is
the Hilbert-Schmidt norm of the density matrix, ``\\sqrt{\\mathrm{tr}(\\rho^\\dagger\\rho)}``,
and not its trace.
"""
norm(state::State) = norm(state.state)


"""
    normalize(::State)

the state rescaled to norm one on a pure representation, to trace one on a mixed one, the
trace being taken as its real part, see `expect`.
"""
normalize(state::State{Pure}) =
    State(state, normalize(state.state))
normalize(state::State{Mixed}) =
    State(state, state.state / real(trace(state)))

# LinearAlgebra has methods of these for any argument, which tried to iterate the state of a
# representation of one's own that has none of its own
norm(state::AbstractState) = throw(MethodError(norm, (state,)))
normalize(state::AbstractState) = throw(MethodError(normalize, (state,)))

"""
    check_same_system(a, b)

refuse two states that do not share their `System`. The indices of a `System` are drawn
afresh, so two systems have none in common even when their sites match, and `ITensorMPS`
would still contract the states, on a deprecated fallback that matches by position and warns.
"""
function check_same_system(a::State, b::State)
    if a.system !== b.system
        error("the two states do not share their System, so their indices differ. " *
              "Use State(system, state)")
    end
    return nothing
end

"""
    inner(a::State, b::State)
    dot(a::State, b::State)
    inner(a::State, op, b::State)
    dot(a::State, op, b::State)

the inner product of two states of the same system and representation, `dot` being an alias
of `inner`, and given an operator placed on sites `op`, the matrix element of `op` between
them, ``\\langle a | op | b \\rangle``, a superoperator on mixed representations.

On pure representations this is the overlap ``\\langle a | b \\rangle``, the first argument
conjugated. On mixed ones it is the Hilbert-Schmidt product ``\\mathrm{tr}(a^\\dagger b)``.
Neither is normalised: divide by the norms, or use `fidelity`.

The two states must share their `System`, see `State(::System, ::State)`. A pure and a mixed
representation are refused.

# Examples

    inner(state, ground_state)
    abs2(inner(a, b))              # the Loschmidt echo of a pure state
    inner(ground_state, X(1) * X(2), excited)
"""
function inner(a::State{R}, b::State{R}) where R
    check_same_system(a, b)
    return dot(a.state, b.state)
end

# Both orders are spelled out on purpose: a catch-all `inner(::State, ::State)` is what Julia
# picks over the diagonal method above, so it would capture the matching pairs as well.
"""
    different_representations()

refuse an inner product between a pure and a mixed representation: one is a vector of the
Hilbert space, the other of the space of operators on it, so there is no product to take.
Refused with a message rather than left to a `MethodError`, which names no way out.
"""
different_representations() =
    error("no inner product between a pure and a mixed representation. Mix the pure " *
          "one, or use fidelity")

inner(::State{Pure}, ::State{Mixed}) = different_representations()
inner(::State{Mixed}, ::State{Pure}) = different_representations()

dot(a::State, b::State) = inner(a, b)

function inner(a::State{R}, op::IndexedOp{R}, b::State{R}) where R
    check_same_system(a, b)
    # checked on the operator as it was written, as expect does
    check_indices(a.system, op)
    return inner(a.state', make_mpo(b, op), b.state)
end

dot(a::State, op::IndexedOp, b::State) = inner(a, op, b)
dot(a::AbstractState, b::AbstractState) = throw(MethodError(dot, (a, b)))

"""
    fidelity(a, b)

the fidelity of two states of the same system, between 0 and 1, whatever the norms and
traces of the arguments.

On two pure representations this is ``|\\langle a | b \\rangle|^2``; on a pure and a mixed
one, in either order, ``\\langle \\psi | \\rho | \\psi \\rangle``. Two mixed representations
are refused: their Uhlmann fidelity ``(\\mathrm{tr}\\sqrt{\\sqrt{\\rho}\\sigma\\sqrt{\\rho}})^2``
needs the spectrum of a density operator, out of reach for a matrix product state. Use
`hs_fidelity` there.

# Examples

    fidelity(state, ground_state)
    measures = "data" => Fidelity(ground_state)
"""
fidelity(a::State{Pure}, b::State{Pure}) =
    abs2(inner(a, b)) / (norm(a)^2 * norm(b)^2)

fidelity(p::State{Pure}, r::State{Mixed}) =
    real(inner(mix(p), r)) / (norm(p)^2 * real(trace(r)))

fidelity(r::State{Mixed}, p::State{Pure}) = fidelity(p, r)

# said here rather than left to a `MethodError`, which names no way out. The docstring
# above carries the reasoning, the message only has to point at it
fidelity(::State{Mixed}, ::State{Mixed}) =
    error("no fidelity between two mixed representations, it needs the spectrum of a " *
          "density operator. Use hs_fidelity")

"""
    hs_fidelity(a::State{Mixed}, b::State{Mixed})

the normalised Hilbert-Schmidt overlap of two mixed states of the same system,
``\\mathrm{tr}(ab)/\\sqrt{\\mathrm{tr}(a^2)\\mathrm{tr}(b^2)}``: the cosine between the two
density matrices seen as vectors, 1 exactly when they are proportional.

This is **not** the Uhlmann fidelity, which is out of reach for a matrix product state, see
`fidelity`, but a cheaper indicator of how close two mixed states are, for instance an
evolution and a reference density matrix.
"""
hs_fidelity(a::State{Mixed}, b::State{Mixed}) =
    real(inner(a, b)) / (norm(a) * norm(b))

"""
    dag(::State)

the state whose density matrix is the adjoint of the given one. A pure representation is
refused.
"""
dag(state::State{Pure}) =
    error("dag is meaningless on pure representations")
function dag(state::State{Mixed})
    n = length(state)
    s = state.system
    if isempty(strong_names(s))
        return State(state, MPS([ tensor_dag(state, i) for i in 1:n ]))
    end
    # `conj` rather than `dag`: the directions stay, only the charges are relabelled
    relab = relabeller(i -> adjoint_index(i, strong_names(s)))
    return State(state, MPS([ relabel(conj(state.state[i]), relab) * adj_map(s, i, relab)
                              for i in 1:n ]))
end


"""
    hermitianize(state [; limits])

the state whose density matrix is the Hermitian part ``(\\rho + \\rho^\\dagger)/2`` of that
of `state`, the sum being truncated according to `limits`, a `Limits` (default `Limits()`). A pure representation is
returned as it is.
"""
hermitianize(state::State{Pure}; limits::Limits=Limits()) =
    state
hermitianize(state::State{Mixed}; limits::Limits=Limits()) =
    State(state, 0.5*(+(state.state, dag(state).state;
                        limits.cutoff, limits.maxdim, limits.mindim)))


"""
    hermiticity(::State)

how Hermitian the density matrix is, as it should be, from 0 when it is anti-Hermitian to 1
when it is Hermitian:
``1/2 + \\mathrm{Re}\\,\\mathrm{tr}(\\rho^2)/(2\\,\\mathrm{tr}(\\rho^\\dagger\\rho))``. 1 on a
pure representation.
"""
hermiticity(::State{Pure}) = 1.
hermiticity(state::State{Mixed}) =
    0.5 + 0.5 * real(dot(state.state, dag(state).state)) / norm(state.state)^2

"""
    renyi2(::State)
    renyi2(::State, positions::AbstractVector{Int})
    renyi2(::State, cut::Int)

the Rényi entropy of order 2, ``-\\log \\mathrm{tr}(\\rho^2)``, of the state, 0 on a pure
representation.

Given positions, that of the state reduced to those sites, which on a pure state measures how
much they are entangled with the rest. A pure state is then mixed first, which is much more
expensive, since a partial trace needs a density matrix. Empty positions give 0, as the cut
0 does. A cut stands for the sites `1:cut`, from 0 to the number of sites; on a pure state it
is read off the entanglement spectrum, as cheap as `entanglement_entropy`.
"""
renyi2(::State{Pure}) = 0.
renyi2(state::State{Mixed}) = -log(trace2(state))

# a number read off a partial trace gives nothing away, so unlike `partial_trace` itself this
# goes through the weak form of a strongly conserving state. No site: the state reduced to
# nothing is its trace, a number, which has no entropy, where `partial_trace` refuses to give
# a state of no site
function renyi2(state::State{Mixed}, a::AbstractVector{Int})
    if isempty(a)
        return 0.0
    end
    return renyi2(partial_trace(weak_form(state), a; keepers = true))
end

# a subsystem of a pure state is not pure, so this is an entanglement measure rather than
# 0. There is no cheap route for an arbitrary subset, the same way mutual_info_renyi2 has
# none: a partial trace needs a density matrix.
function renyi2(state::State{Pure}, a::AbstractVector{Int})
    if isempty(a)
        return 0.0
    end
    return renyi2(mix(state), a)
end

# positions given in a vector of another element type, `[]` or `Any[1, 3]`, are taken as the
# integers they are, rather than met with no method at the first measurement
renyi2(state::State, a::AbstractVector) = renyi2(state, Vector{Int}(a))

"""
    unroll(x)

turn an array of results over the sites, or pairs of sites, each holding one value per
operator, into one such array per operator, arranged as the operators were given. An array
of numbers is returned as it is.
"""
unroll(x) =
    if x[1] isa Number
        x
    else
        v = [ map(t->t[i], x) for i in 1:length(x[1]) ]
        if x[1] isa Matrix
            return reshape(v, size(x[1]))
        else
            return v
        end
    end


"""
    Expector(pos, t)
    Expector()

an expectation value being contracted from left to right: `t` holds the sites up to `pos`,
the tensor of the state on `pos` included, with the site of `pos` left open for the operator
placed there, if only the identity, before moving on. `Expector()` is the empty one, at
`pos = 0`.
"""
struct Expector
    pos::Int
    t::ITensor
end

Expector() =
    Expector(0, ITensor())


"""
    zipto(state, a, i)

the expector `a` carried to site `i`, at or after `a.pos`: the sites in between are traced
out and the tensor of the state on `i` is added, its site left open. From the empty expector
this is `get_left(state, i)`.
"""
function zipto(state::State, a::Expector, i::Int)
    if a.pos == 0
        return Expector(i, get_left(state, i))
    elseif a.pos == i
        return a
    end
    return Expector(i, zip_between(state, a.t, a.pos, i))
end

"""
    zip_between(state, t, from, to)

the contraction `t` of an expector at `from`, its operator placed there, carried to `to`: the
site `from` closed, the sites in between traced out, and the tensor of the state on `to`
added, its site left open
"""
function zip_between(state::State{Pure}, t::ITensor, from::Int, to::Int)
    st = state.state
    t *= dag(st[from]')
    for k in from+1:to-1
        t *= st[k] * identity_at(state, k)
        t *= dag(st[k]')
    end
    return t * st[to]
end

function zip_between(state::State{Mixed}, t::ITensor, from::Int, to::Int)
    for k in from+1:to-1
        t *= get_loc(state, k)
    end
    return t * state.state[to]
end

"""
    zipend(state, a)

the expector `a` completed with the environment on the right of `a.pos`, see `get_right`.
"""
zipend(state::State, a::Expector) =
    Expector(a.pos, a.t * get_right(state, a.pos))


"""
    expectfactor(state, a, o)
    expectfactor(state, a, t, i)

the expector `a` carried to the site of the one site operator `o`, and multiplied there by
the tensor measuring it, see `tensor_obs`; or carried to site `i` and multiplied by `t`. A
`Multi_F` places `F` on each of its sites.
"""
expectfactor(state::State, a::Expector, o::AtIndex) =
    expectfactor(state, a, tensor_obs(state, o), only(o.index))

function expectfactor(state::State, a::Expector, t::ITensor, i::Int)
    a = zipto(state, a, i)
    Expector(a.pos, a.t * t)
end

function expectfactor(state::State, a::Expector, o::Multi_F)
    for i in o.start:o.stop
        a = expectfactor(state, a, F(i))
    end
    return a
end

"""
    times_piece(state, t, o, k)

the expector tensor `t`, carried to site `k`, times the piece `o` of a com measured there, see
`obs_at`. On a pure state the identity, the most frequent piece of a com, only renames the
index of the site, which is much cheaper than contracting a delta.
"""
times_piece(state::State{Pure}, t::ITensor, o::SimpleOp, k::Int) =
    if scalararg(o) isa IdentityOp
        scalarcoef(o) * prime(t, SysIndex{Pure}(state.system, k))
    else
        t * scalarcoef(o) * obs_at(state, scalararg(o), k)
    end

times_piece(state::State{Mixed}, t::ITensor, o::SimpleOp, k::Int) =
    t * scalarcoef(o) * obs_at(state, scalararg(o), k)


"""
    expect_norm(state, obs)
    expect_norm(state, coef, factors)
    expect_norm(state, coef, com)

`expect` without the simplification nor the normalization: ``\\mathrm{tr}(A\\rho)``, or
``\\langle\\psi|A|\\psi\\rangle``, of an operator already in the form `simplify` gives, or the
array of those of an array of them; given `coef` and `factors`, that of their product, and given
a com, that of the com times `coef`. A com is
contracted from left to right with one expector per channel: a channel opens on `get_left`,
goes from site to site by `zip_between` and closes on `get_right`. `measure` calls it on
operators `make_obs` has simplified once and for all.
"""
function expect_norm(state::State, coef::Number, subs::Vector{<:IndexedOp{Pure}})
    if coef == 0.
        return 0.
    end
    state = weak_form(state)
    # every expectation value reaches this leaf, `measure` included, which calls
    # `expect_norm` rather than `expect`. Checking here covers them all at the cost of a
    # few integer comparisons per term, nothing next to the contractions below
    foreach(o -> check_indices(state.system, o), subs)
    foreach(o -> check_one_site(o, "expect"), subs)
    # the identity has no site, and the trace is what it gives
    subs = filter(o -> !(o isa IdentityOp), subs)
    if isempty(subs)
        return coef * trace(state)
    end
    e = Expector()
    for o in subs
        e = expectfactor(state, e, o)
    end
    e = zipend(state, e)
    return coef * scalar(e.t)
end

function expect_norm(state::State, coef::Number, a::ComOp{Pure})
    state = weak_form(state)
    check_indices(state.system, a)
    open = Union{Nothing, ITensor}[]
    total = 0.
    for (j, (k, ps)) in enumerate(zip(com_sites(a), a.pieces))
        carried = [ isnothing(t) ? nothing : zip_between(state, t, k - 1, k) for t in open ]
        next = Vector{Union{Nothing, ITensor}}(nothing, get(a.linkdims, j, 0))
        for (l, r, o) in ps
            from = l == 0 ? get_left(state, k) : carried[l]
            if isnothing(from)
                continue
            end
            t = times_piece(state, from, o, k)
            if r == 0
                total += scalar(t * get_right(state, k))
            else
                # out of place: a piece may share its storage with `get_left`, cached with the
                # state, which adding in place would corrupt
                next[r] = isnothing(next[r]) ? t : next[r] + t
            end
        end
        open = next
    end
    return coef * total
end

expect_norm(state::State, coef::Number, a::IndexedOp{Pure}) =
    expect_norm(state, coef, prodsubs(a))

"""
    expect(state, obs; normalize = true)

the expectation value of `obs`, an operator placed on sites, or the array of those of an
array of them. It is divided by the trace of the state, its squared norm on a pure
representation, so the state need not be normalised. On a mixed representation the trace is
taken as its real part, its imaginary part being numerical error on a density matrix.

With `normalize = false` it is not divided: ``\\mathrm{tr}(A\\rho)``, or
``\\langle\\psi|A|\\psi\\rangle``, as it is. This is what an operator made non Hermitian by
`Left` or `Right` needs, ``B\\rho`` in a correlation at two times for instance, whose trace
has an imaginary part of its own or is zero.

`obs` is simplified first, so its factors may be given in any order and the Jordan-Wigner
strings of fermionic operators are inserted for you. An operator not placed on sites, `X`
rather than `X(1)`, is refused, `expect1` measuring it on every site, and so is a
superoperator.

A representation of one's own measures operators through `expect`, see
[Representations of one's own](@ref).

# Examples

    expect(state, X(1)*Y(2) + Y(1)*Z(3))
    expect(state, [X(1)*Y(2), X(3), Z(1)*X(2)])
    expect(state, C(3)*dag(C)(1))
"""
function expect(state::State, op::Op; normalize::Bool = true)
    # on the operator as it was written, as make_mpo and apply do: simplify places an identity
    # on the first site whatever site it was given, and its Jordan-Wigner strings would be
    # named in the message rather than what the caller wrote
    check_indices(state.system, op)
    v = expect_norm(state, simplify(op))
    return normalize ? v / real(trace(state)) : v
end

expect(state::AbstractState, ops::Union{AbstractArray, Tuple}; kwargs...) =
    map(ops) do o
        expect(state, o; kwargs...)
    end



expect_norm(state::State, p::IndexedOp{Pure}) =
    expect_norm(state, scalarcoef(p), scalararg(p))

expect_norm(state::State, op::SumOp{Pure, Indexed}) =
    sum(op.subs) do p
        expect_norm(state, p)
    end

# said here rather than left to the method below, which would try to iterate the operator
expect_norm(::State, a::Op{Mixed}) =
    error("expect takes an observable, and $a is a superoperator acting on a density matrix")

expect_norm(::State, a::GenericOp{Pure, N}) where N =
    error("expect takes an operator placed on sites, such as $(a((1:N)...)) rather than $a, " *
          "and expect1 measures a one site operator on every site")

expect_norm(state::State, op) =
    map(op) do o
        expect_norm(state, o)
    end

"""
    expect1_one(state, op, i, t)

the expectation value of the one site operator `op` on site `i`, or the array of those of an
array of operators, `t` being the state contracted on every site with that of `i` left open.
A fermionic operator is refused: its expectation value vanishes on any state of definite
fermion parity.
"""
expect1_one(state::State, op::SimpleOp, i::Int, t::ITensor) =
    if isfermionic(op)
        # this is a refusal, not a gap to fill: implementing it would add a path whose
        # only correct answer is zero, and would silently open the branch `expect2` picks
        # on the parity of its first operator
        error("expect1 does not take a fermionic operator: its expectation value is odd, " *
              "so it vanishes on any state of definite fermion parity. If you really " *
              "want it on a state that superposes parities, ask for it site by site " *
              "with expect(state, op(i))")
    else
        scalar(t * obs_at(state, op, i))
    end

expect1_one(state::State, ops, i::Int, t::ITensor) =
    map(ops) do o
        expect1_one(state, o, i, t)
    end


"""
    expect1(state, op)

the expectation values of the one site operator `op` on every site, as a vector indexed by
site, or, for an array of operators, an array of such vectors arranged as the operators. They
are normalised as by `expect`.

A fermionic operator is refused: its expectation value vanishes on any state of definite
fermion parity. On a state superposing parities, measure it site by site with
`expect(state, op(i))`.

# Examples

    expect1(state, X)
    expect1(state, [X, Y, Z])
"""
function expect1(state::State, op)
    state = weak_form(state)
    n = length(state)
    t = real(trace(state))
    r = [ expect1_one(state, op, i, get_left(state, i) * get_right(state, i) / t) for i in 1:n ]
    return unroll(r)
end


function expect2(state::State, ops::Vector{<:Tuple{SimpleOp, SimpleOp}})
    # the cross terms below pick their branch on the parity of the first operator alone,
    # and the sign a swap costs takes the second to have the same one. A pair mixing the
    # two is refused here rather than left to the diagonal, which rejects it today only
    # because `expect1_one` has no fermionic case of its own.
    for (o1, o2) in ops
        if isfermionic(o1) ≠ isfermionic(o2)
            error("cannot correlate $o1 and $o2: one is fermionic and the other is not, " *
                  "so their product is odd and has no expectation value. Both operators " *
                  "of a pair must have the same fermionic parity")
        end
    end
    state = weak_form(state)
    oplist = [first.(ops) ; last.(ops)]
    need_fermionic = any(isfermionic, oplist)
    need_non_fermionic = any(x->!isfermionic(x), oplist)
    n = length(state)
    # `Matrix{Any}` is said out loud rather than left to the bare `Matrix(undef, ...)`,
    # which means the same thing without showing it. The element type is not known before a
    # cell is computed, since it follows what `scalar` returns for this state, and `unroll`
    # rebuilds a concretely typed result at the end, so the untyped container stays internal
    r = Matrix{Any}(undef, n, n)
    scale = real(trace(state))
    for i in 1:n
        # divided by the trace once, every correlation of site i being built on it
        e = zipto(state, Expector(), i)
        lnf = Expector(e.pos, e.t / scale)
        t = zipend(state, lnf).t
        r[i, i] = map(ops) do (o1, o2)
            expect1_one(state, o1 * o2, i, t)
        end
        lf = lnf
        for j in i+1:n
            enf = need_non_fermionic ? zipend(state, zipto(state, lnf, j)) : lnf
            ef = need_fermionic ? zipend(state, zipto(state, lf, j)) : lf
            # the expectation value of a(i) * b(j), a and b of the same parity: a fermionic
            # a takes its string through the F of its site
            correlation(a, b) =
                isfermionic(a) ?
                    scalar(ef.t * tensor_obs(state, (a * F)(i)) * obs_at(state, b, j)) :
                    scalar(enf.t * obs_at(state, a, i) * obs_at(state, b, j))
            r[i, j] = map(((o1, o2),) -> correlation(o1, o2), ops)
            # swapping two fermionic operators costs a sign
            r[j, i] = map(ops) do (o1, o2)
                isfermionic(o1) ? -correlation(o2, o1) : correlation(o2, o1)
            end
            if j < n
                if need_non_fermionic
                    lnf = zipto(state, expectfactor(state, lnf, identity_at(state, j), j), j+1)
                end
                if need_fermionic
                    lf = zipto(state, expectfactor(state, lf, F(j)), j+1)
                end
            end
        end
    end
    unroll(r)
end

"""
    expect2(state, (o1, o2))
    expect2(state, [(o1, o2), ...])

the correlations of a pair of one site operators on every pair of sites: the matrix whose
entry `(i, j)` is the expectation value of `o1(i) * o2(j)`, the diagonal being that of
`(o1 * o2)(i)`. For a vector of pairs, a vector of such matrices. They are normalised as by
`expect`.

The two operators of a pair must both be fermionic or both not, their product being odd
otherwise; the Jordan-Wigner strings of fermionic ones are inserted for you.

# Examples

    expect2(state, (X, X))
    expect2(state, [(X, Y), (X, Z), (Y, Z)])
"""
expect2(state::State, ops::Tuple{SimpleOp, SimpleOp}) =
    expect2(state, [ops])[1]

"""
    variance(hamiltonian, ::State{Pure})
    variance(::MPO, ::State{Pure})

the variance of the energy, ``\\langle H^2 \\rangle - \\langle H \\rangle^2``, zero exactly
when the state is an eigenstate of the hamiltonian. The state need not be normalised, and a
mixed representation is refused.

This is the convergence check of a ground state search. The `tolerance` of `GroundState`
stops when the energy no longer progresses between two sweeps, which a search stuck in a
metastable state also does, with a large variance. Extrapolating the energy to zero variance
over several bond dimensions also gives an error bar.

``H^2`` is never formed: its terms are the products of pairs of those of `H`, and the product
of the MPO of `H` by itself has the square of its bond dimension. What is computed is
``\\langle H\\psi | H\\psi \\rangle``, with the MPO of `H` on either side, at the cost of a
`dmrg` sweep at the same bond dimension: this belongs in `final_measures`, or under a large
`measures_period`, rather than at every sweep. Pass an `MPO` to reuse one already built.

# Examples

    energy, gs = dmrg(hamiltonian, state; nsweeps = 10)
    variance(hamiltonian, gs)
    measures = "data" => Variance(hamiltonian)
"""
function variance(mpo::MPO, state::State{Pure})
    st = state.state
    n2 = real(dot(st, st))
    # the bra is primed because that is the form `inner` wants: contracting the mpo with
    # the ket leaves the site indices primed, and `inner(x, A, y)` with an unprimed `x`
    # takes a fallback that ITensorMPS deprecated and says it will turn into an error
    e = real(inner(prime(st), mpo, st)) / n2
    e2 = real(inner(mpo, st, mpo, st)) / n2
    return e2 - e^2
end

variance(h, state::State{Pure}) = variance(make_mpo(state, h), state)

# a hamiltonian on a density matrix does not have this reading, and the quantity that
# plays the part is already there: `steady_state` optimises on ``L^\\dagger L`` and returns
# ``\\|L\\rho\\|^2``, which is zero exactly when the state is stationary
variance(_, ::State{Mixed}) =
    error("variance needs a pure representation. On a mixed one steady_state returns " *
          "the equivalent for a Lindbladian")


"""
    cut_svd(state, pos)

the factors `U` and `S` of the singular value decomposition of the state across the cut
between sites `pos` and `pos + 1`, `U` holding the sites up to `pos`; a position the state does
not have is refused
"""
function cut_svd(state::State, pos::Int)
    n = length(state)
    if !(1 ≤ pos ≤ n)
        error("cannot compute the entanglement entropy at site $pos of a $n site state, " *
              "the cut is on the right of a site so pos must be between 1 and $n")
    end
    s = orthogonalize(state.state, pos)
    U, S = svd(s[pos], (linkinds(s, pos-1)..., siteinds(s, pos)...))
    return U, S
end

"""
    shannon_entropy(p)

the entropy ``-\\sum_i p_i \\log p_i`` of the probabilities `p`. A probability of exactly zero,
which a singular value kept by a `mindim` above the Schmidt rank gives, adds nothing, where
`0 * log(0)` would make it NaN.
"""
shannon_entropy(p) = -sum(x * log(x) for x in p if x > 0)

"""
    entropy_spectrum(S)

the entanglement entropy and spectrum read off the singular values `S` of a cut, as
`entanglement_entropy` returns them
"""
function entropy_spectrum(S::ITensor)
    # sorted, because ITensors reads the singular values off the diagonal, and on the block
    # sparse tensor of a state that conserves something that diagonal runs block by block:
    # the spectrum would come out grouped by sector rather than decreasing
    sp = sort([ S[i,i]^2 for i in 1:dim(S, 1) ]; rev = true)
    sp /= sum(sp)
    return (shannon_entropy(sp), sp)
end

"""
    entanglement_entropy(state, cut::Int)

the entanglement entropy across the cut between sites `cut` and `cut + 1`, together with the
spectrum it is computed from: the eigenvalues of the reduced density matrix of the sites up
to `cut`, that is the squared singular values of the cut normalised to sum to one, in
decreasing order.

`cut` runs from 1 to the number of sites; the last cut leaves nothing on its right and
always gives 0. On a mixed representation the same computation, on the density matrix seen
as a vector, gives the operator space entanglement entropy (OSEE).

# Examples

    ee, spectrum = entanglement_entropy(state, 3)   # cut between sites 3 and 4
"""
entanglement_entropy(state::State, cut::Int) = entropy_spectrum(cut_svd(state, cut)[2])

"""
    entanglement_by_sector(state::State{Pure}, cut::Int)

the entanglement across the cut between sites `cut` and `cut + 1`, resolved by the charge the
sites up to `cut` carry: a `Dict` from each charge, a `QN` as `flux` gives it, to the named
tuple `(weight, entropy, spectrum)`, `weight` being the probability of that charge and
`entropy` and `spectrum` those of the reduced density matrix restricted to it and normalised,
the spectrum in decreasing order.

The entanglement entropy is `Σ weight * entropy - Σ weight * log(weight)`, the second sum,
the number entropy, coming from the charge fluctuating across the cut. A state conserving
nothing has the single sector `QN()`. `QN` comes from ITensors, `using ITensors: QN`. A mixed
representation is refused.

# Examples

    using ITensors: QN
    sectors = entanglement_by_sector(state, 3)   # cut between sites 3 and 4
    sectors[QN("N", 2)].weight
"""
function entanglement_by_sector(state::State{Pure}, cut::Int)
    U, S = cut_svd(state, cut)
    u = commonind(U, S)
    if !hasqns(u)
        ee, sp = entropy_spectrum(S)
        return Dict(QN() => (weight = 1.0, entropy = ee, spectrum = sp))
    end
    # U carries no flux, so what flows into it through `u` is what the sites up to `cut` hold
    if flux(U) ≠ QN()
        error("bug: the left factor of a cut carries $(flux(U)), so its sectors cannot be read")
    end
    squares = Dict{QN, Vector{Float64}}()
    k = 0
    for (q, d) in space(u)
        append!(get!(squares, dir(u) == ITensors.In ? q : -q, Float64[]),
                [ S[k + i, k + i]^2 for i in 1:d ])
        k += d
    end
    total = sum(sum, values(squares))
    sectors = Dict{QN, NamedTuple{(:weight, :entropy, :spectrum),
                                  Tuple{Float64, Float64, Vector{Float64}}}}()
    for (q, sq) in squares
        w = sum(sq)
        # a sector the decomposition kept with nothing in it has no probability of occurring
        if w == 0
            continue
        end
        sp = sort(sq / w; rev = true)
        sectors[q] = (weight = w / total, entropy = shannon_entropy(sp), spectrum = sp)
    end
    return sectors
end

entanglement_by_sector(::State{Mixed}, ::Int) =
    error("entanglement_by_sector takes a pure state: on a mixed one the sectors of a link are " *
          "differences between the charges of the ket and of the bra")

"""
    check_positions(state, pos, what)

refuse a position of `pos` the state does not have, `what` naming the function given it
"""
function check_positions(state::State, pos, what)
    n = length(state)
    for p in pos
        if !(1 ≤ p ≤ n)
            error("$what was given site $p, which the state does not have: it has $n sites")
        end
    end
    return nothing
end

"""
    trace_signs(state, keep)

the state whose partial trace keeping the sites `keep` is the reduced state of `state`, its
fermionic signs included. In the Jordan-Wigner basis a fermion traced out has to be moved past
every fermion kept on its right, each one giving a sign: the reduced state is the trace of
``U \\rho U^\\dagger``, ``U = (-1)^{\\sum n_l n_k}``, with `l` running over the fermionic sites
traced out and `k` over the fermionic sites kept on their right, and `n` the parity of a site.
``U \\rho U^\\dagger`` is the product of the state by an MPO of bond dimension 2, which carries
the parity of the fermions traced out on the left: a site traced out projects on each parity,
on the side of the ket alone, being traced afterwards, and a site kept takes `Gate(F)` when that
parity is odd. A state needing no sign is given back as it is.
"""
function trace_signs(state::State{Mixed}, keep)
    sys = state.system
    n = length(sys)
    odd = [ matrix(F, sys[i]) != I for i in 1:n ]
    if !any(l -> odd[l] && l ∉ keep && any(k -> odd[k] && k > l, keep), 1:n)
        return state
    end
    links = [ is_charged(sys) ? Index([QN() => 2]; tags = "Parity,l=$i") : Index(2, "Parity,l=$i")
              for i in 0:n ]
    projs = [ Left((Id + F) / 2), Left((Id - F) / 2) ]
    mps = state.state
    ts = map(1:n) do i
        idx = SysIndex{Mixed}(sys, i)
        l, r = links[i], links[i+1]
        w = ITensor(idx', dag(idx), dag(l), r)
        id = delta(dag(idx), idx')
        for p in 1:2
            if !odd[i]
                add_block!(w, l, p, r, p, id, idx)
            elseif i in keep
                add_block!(w, l, p, r, p, p == 1 ? id : tensor(sys, Gate(F)(i)), idx)
            else
                for q in 1:2
                    add_block!(w, l, p, r, xor(p - 1, q - 1) + 1, tensor(sys, projs[q](i)), idx)
                end
            end
        end
        return noprime(w * mps[i])
    end
    # no fermion on the left of the first site, and either parity on the right of the last
    ts[1] *= onehot(links[1] => 1)
    ts[n] *= onehot(dag(links[n+1]) => 1) + onehot(dag(links[n+1]) => 2)
    for i in 1:n-1
        c = combiner(commonind(mps[i], mps[i+1]), links[i+1]; tags = "Link,l=$i")
        ts[i] *= c
        ts[i+1] *= dag(c)
    end
    return State(state, MPS(ts))
end

"""
    partial_trace(state, positions::AbstractVector{Int} [; keepers = false])

the state with the sites at `positions` traced out, or, with `keepers = true`, all the others.
The result is a mixed state on a new system made of the sites kept, in their order, and it
has the trace of `state`. It keeps the fermionic signs: an operator of the sites kept has on it
the expectation value it has on `state`, its Jordan-Wigner strings crossing the fermions traced
out. Tracing out a fermionic site with fermionic sites kept on its right may double the bond
dimension of the result.

A pure representation is refused, the reduced state being mixed in general: use `mix(state)`
first. So is a state conserving something strongly, the reduced state spreading over several
sectors: weaken it first. Tracing out every site is refused too.

# Examples

    partial_trace(state, [1, 2])
    partial_trace(state, 3:5; keepers = true)
"""
function partial_trace(state::State{Mixed}, pos::AbstractVector{Int}; keepers::Bool = false)
    if !isempty(strong_names(state.system))
        error("cannot trace out part of a state conserving something strongly, what is left " *
              "spreads over several sectors: weaken it first, giving up the strong symmetry")
    end
    n = length(state)
    # a position the state does not have would be silently ignored when tracing, the filter
    # below never meeting it, and would raise a BoundsError when keeping
    check_positions(state, pos, "partial_trace")
    if keepers
        keep = sort(unique(pos))
    else
        keep = filter(e->e ∉ pos, 1:n)
    end
    kn = length(keep)
    if kn == 0
        error("partial_trace cannot trace all sites of a state")
    end
    state = trace_signs(state, keep)
    mps = state.state
    sys = state.system
    j = 0
    t = Vector{ITensor}(undef, kn)
    for (i, k) in enumerate(keep)
        if i == 1
            # the sites on the left traced out, without the 1/trace the left environment of
            # expect carries: a partial trace keeps the trace of the state, and a traceless
            # state gave NaN
            x = ITensor(1.)
            for l in 1:k-1
                x *= get_loc(state, l)
            end
            t[1] = x * mps[k]
        else
            for l in j+1:k-1
                t[i-1] *= get_loc(state, l)
            end 
            t[i] = copy(mps[k])
        end
        j = k
    end
    t[kn] *= get_right(state, keep[kn])
    for i in 1:kn-1
        idx = commonind(t[i], t[i+1])
        jdx = settags(idx, "Link, l=$i")
        replaceind!(t[i], idx, jdx)
        replaceind!(t[i+1], idx, jdx)
    end
    kept = System(AbstractSite[ sys[k] for k in keep ],
                  Index[ SysIndex{Pure}(sys, k) for k in keep ],
                  Index[ SysIndex{Mixed}(sys, k) for k in keep ])
    return State{Mixed}(kept, MPS(t))
end

# a partial trace needs a density matrix: the reduced state of a subsystem is mixed in
# general, so there is nothing to hand back in pure representation. `renyi2` and
# `mutual_info_renyi2` mix on their own because they return a number; this one returns a
# state, and changing its representation behind the caller's back would be a surprise
partial_trace(::State{Pure}, ::AbstractVector{Int}; kwargs...) =
    error("partial_trace needs a mixed representation, use mix(state) first")

"""
    check_cut(state, cut, what)

refuse a cut the state does not have, `what` naming the function it was given to: `cut` counts
the sites on its left, from 0 to the number of sites, and at either end one part is empty.
"""
function check_cut(state::State, cut::Int, what)
    n = length(state)
    if !(0 ≤ cut ≤ n)
        error("$what was given the cut $cut, which the state does not have: " *
              "a cut counts the sites on its left, from 0 to $n")
    end
    return nothing
end

# a cut, the sites `1:cut`, as `mutual_info_renyi2` takes it. At either end one part is
# empty: nothing, which has no entropy, or the whole state
function renyi2(state::State, cut::Int)
    check_cut(state, cut, "renyi2")
    if cut == 0
        return 0.0
    elseif cut == length(state)
        return renyi2(state)
    end
    return renyi2(state, collect(1:cut))
end

# the sites on the left of a cut share their Schmidt spectrum with the rest: their entropy is
# read off it, as `mutual_info_renyi2` does, rather than from a partial trace of a mix
function renyi2(state::State{Pure}, cut::Int)
    check_cut(state, cut, "renyi2")
    if cut == 0 || cut == length(state)
        return 0.0
    end
    return -log(sum(abs2, last(entanglement_entropy(state, cut))))
end

"""
    mutual_info_renyi2(state::State, cut::Int)
    mutual_info_renyi2(state::State, a::AbstractVector{Int})

the Rényi-2 analogue of the mutual information between two parts of the state,
``S_2(A) + S_2(B) - S_2(A \\cup B)``, which unlike the mutual information can be negative on a
mixed state. Part A is given by its positions, or by a cut, sites `1:cut`, from 0 to the number
of sites; part B is the rest. Positions that are empty, or cover every site, give 0.

On a pure state and for a cut, it is read off the entanglement spectrum and costs no more
than `entanglement_entropy`. For positions a pure state is mixed first, which is much more
expensive.
"""
function mutual_info_renyi2(state::State, a::AbstractVector{Int})
    n = length(state)
    check_positions(state, a, "mutual_info_renyi2")
    # a part that is all the system or nothing shares nothing with the rest, which has no site
    # for partial_trace to keep
    k = length(unique(a))
    if k == 0 || k == n
        return 0.0
    end
    w = weak_form(state)
    return renyi2(partial_trace(w, a; keepers = true)) +
           renyi2(partial_trace(w, a; keepers = false)) -
           renyi2(w)
end

function mutual_info_renyi2(state::State, cut::Int)
    check_cut(state, cut, "mutual_info_renyi2")
    return mutual_info_renyi2(state, collect(1:cut))
end

# a pure state has no entropy of its own, and the two sides of a cut share their Schmidt
# spectrum, so the mutual information is just twice the renyi2 entropy of either side. The cut
# is checked here for the message to name this function
function mutual_info_renyi2(state::State{Pure}, cut::Int)
    check_cut(state, cut, "mutual_info_renyi2")
    return 2 * renyi2(state, cut)
end

# partial_trace needs a density matrix, there is no cheap route for an arbitrary subset
mutual_info_renyi2(state::State{Pure}, a::AbstractVector{Int}) =
    mutual_info_renyi2(mix(state), a)

# as for `renyi2`, positions in a vector of another element type are taken as integers
mutual_info_renyi2(state::State, a::AbstractVector) = mutual_info_renyi2(state, Vector{Int}(a))



"""
    draw(p, d, rnd)

the outcome among `0:d-1` that `rnd` falls on, `p(x)` being the probability of `x`: the first
whose cumulated probability exceeds `rnd`, the last one taking what the others leave
"""
function draw(p, d::Int, rnd)
    ptot = 0.
    for x in 0:d - 2
        ptot += p(x)
        if rnd < ptot
            return x
        end
    end
    return d - 1
end

"""
    sample(::State [; rng])
    sample(::State, pos::Int [; rng])

a random outcome of measuring the state in the computational basis, outcomes being numbered
from 0: a vector with one outcome per site, drawn with the correlations between sites, or,
given `pos`, the outcome of that site alone, the others being traced out. `rng` is the random
number generator, the global one by default.

# Examples

    sample(state)
    sample(state, 3)
"""
function sample(state::State{Pure}; rng = Random.default_rng())
    st = orthogonalize(state.state, 1)
    st[1] /= norm(st[1])
    return sample(rng, st) .- 1
end

function sample(state::State{Pure}, pos::Int; rng = Random.default_rng())
    st = orthogonalize(state.state, pos)
    t = st[pos] / norm(st[pos])
    r = rand(rng)
    ind = siteind(st, pos)
    # the conjugate leg, which on a charged site is the direction the tensor takes
    return draw(dim(ind), r) do x
        tx = t * onehot(dag(ind) => x + 1)
        real(scalar(tx * dag(tx)))
    end
end

function sample(state::State{Mixed}, pos::Int; rng = Random.default_rng())
    state = weak_form(state)
    sys = state.system
    l = get_left(state, pos) / real(trace(state))
    r = get_right(state, pos)
    d = dim(SysIndex{Pure}(sys, pos))
    rnd = rand(rng)
    return draw(x -> real(scalar(l * tensor_obs(state, Proj(x)(pos)) * r)), d, rnd)
end

function sample(state::State{Mixed}; rng = Random.default_rng())
    state = weak_form(state)
    sys = state.system
    n = length(state)
    result = Vector{Int}(undef, n)
    l = ITensor(1.)
    for pos in 1:n
        a = l * state.state[pos]
        r = get_right(state, pos)
        d = dim(SysIndex{Pure}(sys, pos))
        tot = real(scalar(a * identity_at(state, pos) * r))
        rnd = rand(rng) * tot
        x = draw(i -> real(scalar(a * tensor_obs(state, Proj(i)(pos)) * r)), d, rnd)
        result[pos] = x
        l = a * tensor_obs(state, Proj(x)(pos))
    end
    return result
end