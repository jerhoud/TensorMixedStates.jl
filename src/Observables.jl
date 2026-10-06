# What is computed on a state: traces, norms, inner products and fidelities, expectation values
# and correlations, entropies, partial traces and sampling.

export trace, trace2, norm, normalize, hermitianize, hermiticity, renyi2
export inner, dot, fidelity, hs_fidelity
export expect, expect1, expect2
export entanglement_entropy, entanglement_by_sector, partial_trace, mutual_info_renyi2, sample, variance

"""
    weak_form(state)

the state measurements run on: `state` itself, or, for a mixed state conserving something
strongly, `weaken(state)`, computed once and kept with it. Only the weak form has a trace
with one vector per site; weakening changes no expectation value.
"""
weak_form(state::State{Pure}) = state

function weak_form(state::State{Mixed})
    if isempty(strong_names(state.system))
        return state
    end
    w = state.preobs.weak
    lock(state.preobs.lock) do
        if isempty(w)
            push!(w, weaken(state))
        end
    end
    return only(w)
end

"""
    strong_measured(i)

raise the bug error for site `i` of a strongly conserving state measured without going
through `weak_form`
"""
strong_measured(i::Int) =
    error("bug: site $i of a strongly conserving state measured without going through weak_form")

"""
    on_trace(state, t, i)

the tensor `t` of a pure operator on site `i` carried onto the mixed index there, so that
contracting it with the density matrix applies the operator and traces the site out
"""
function on_trace(state::State{Mixed}, t::ITensor, i::Int)
    s = state.system
    j = SysIndex{Pure}(s, i)
    k = SysIndex{Mixed}(s, i)
    b, c = mixer(j, k, s[i])
    if b !== j
        strong_measured(i)
    end
    # daggered to meet the `k` of the state
    return t * dag(c)
end

"""
    identity_at(state, i)

the tensor measuring the identity on site `i`, see `tensor_obs`: the identity of the site on a
pure state, the trace over the site on a mixed one, conserving nothing strongly, see
`weak_form`
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

the tensor measuring the one site operator `o`, possibly times a coefficient, on `state`: that
of the operator on a pure representation, carried by `on_trace` on a mixed one
"""
tensor_obs(state::State{Pure}, ind::AtIndex{Pure, 1}) =
    tensor(state.system, ind)

tensor_obs(state::State{Mixed}, ind::AtIndex{Pure, 1}) =
    on_trace(state, tensor(state.system, ind), only(ind.index))

# `(c * A)(i)` keeps its coefficient outside the AtIndex
tensor_obs(state::State, a::ScalarOp{Pure, Indexed, 1}) =
    a.coef * tensor_obs(state, a.arg)

"""
    obs_at(state, op, i)

the tensor measuring the one site operator `op` on site `i`, see `tensor_obs`, a multiple of the
identity included, see `identity_at`
"""
obs_at(state::State, op::SimpleOp, i::Int) =
    scalararg(op) isa IdentityOp ? scalarcoef(op) * identity_at(state, i) : tensor_obs(state, op(i))

"""
    tensor_dag(state, i)

the tensor on site `i` of the adjoint of the density matrix, for a state conserving nothing
strongly, see `dag`
"""
function tensor_dag(state::State, i::Int)
    s = state.system
    j = SysIndex{Pure}(s, i)
    k = SysIndex{Mixed}(s, i)
    c = last(mixer(j, k, s[i]))
    # ket and bra exchanged by renaming: a combiner in the other order would ask for a mixed
    # index of the opposite charge
    return replaceinds(dag(state.state[i]) * c, (dag(j), j'), (dag(j'), j)) * c
end

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
    create_loc!(l, state)

fill `l` with what `get_loc` returns for every site, mixed states only
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

the tensor of the density matrix on site `i` with that site traced out, cached for all sites
on first use, mixed states only
"""
get_loc(state::State, i::Int) = cached!(create_loc!, state.preobs.loc, state)[i]

"""
    create_right!(r, state)

fill `r` with what `get_right` returns for every site; on a pure state, the right orthogonal
sites from `rightlim` on give the identity, uncontracted
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

the environment on the right of site `i`: the sites after `i` traced out on a mixed state,
ket with bra on a pure one, the bra of site `i` included; cached for all sites on first use
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

the trace of the density matrix, which should be one, the squared norm of the state on a pure
representation. It is complex when the tensors of the state are, after an evolution for
instance, its imaginary part being numerical error; `Trace` measures the real part.
"""
function trace(state::State)
    state = weak_form(state)
    return cached!(create_trace!, state.preobs.trace, state)[1]
end

"""
    extend_left!(l, state, i)

extend `l`, which holds what `get_left` returns for its first sites, the first one at least,
up to site `i`; on a pure state, the left orthogonal sites up to `leftlim` give the identity,
uncontracted
"""
function extend_left!(l, state::State{Pure}, i::Int)
    st = state.state
    ll = ITensorMPS.leftlim(st) + 1
    j = length(l)
    resize!(l, i)
    for k in j+1:i
        llink = commonind(st[k-1], st[k])
        v = if k <= ll
            # the directions the branch below leaves, the ket link as st[k-1] holds it and the
            # bra one daggered and primed: only charges tell the other way round apart
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

the environment on the left of site `i`: the sites before `i` traced out on a mixed state, ket
with bra on a pure one, times the tensor of the state on `i`, not divided by the trace;
computed up to `i` on demand and cached
"""
function get_left(state::State, i::Int)
    l = state.preobs.left
    # read under the lock too, the vector growing as it is extended
    return lock(state.preobs.lock) do
        if isempty(l)
            push!(l, state.state[1])
        end
        if length(l) < i
            extend_left!(l, state, i)
        end
        l[i]
    end
end

"""
    trace2(::State)

the purity ``\\mathrm{tr}(\\rho^2)`` of the density matrix normalised to trace one: 1 for a
pure state, less for a mixed one.
"""
trace2(::State{Pure}) = 1.
trace2(state::State{Mixed}) = (norm(state.state) / real(trace(state))) ^ 2

"""
    norm(::State)

the norm of the state, which should be one on a pure representation, and on a mixed one the
Hilbert-Schmidt norm ``\\sqrt{\\mathrm{tr}(\\rho^\\dagger\\rho)}``, not the trace.
"""
norm(state::State) = norm(state.state)


"""
    normalize(::State)

the state rescaled to norm one on a pure representation, to a trace of real part one on a
mixed one.
"""
normalize(state::State{Pure}) =
    State(state, normalize(state.state))
normalize(state::State{Mixed}) =
    State(state, state.state / real(trace(state)))

# the LinearAlgebra fallbacks would iterate a representation of one's own
norm(state::AbstractState) = throw(MethodError(norm, (state,)))
normalize(state::AbstractState) = throw(MethodError(normalize, (state,)))

"""
    check_same_system(a, b)

refuse two states that do not share their `System`, whose indices differ even when their
sites match
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

the inner product of two states, ``\\langle a | b \\rangle`` on pure representations,
``\\mathrm{tr}(a^\\dagger b)`` on mixed ones, not normalised (see `fidelity`); given an
operator placed on sites `op`, a superoperator on mixed representations, the matrix element
``\\langle a | op | b \\rangle``. `dot` is an alias of `inner`.

The two states must share their `System`, see `State(::System, ::State)`, and their
representation.

# Examples

    inner(state, ground_state)
    abs2(inner(a, b))              # the Loschmidt echo of a pure state
    inner(ground_state, X(1) * X(2), excited)
"""
function inner(a::State{R}, b::State{R}) where R
    check_same_system(a, b)
    return dot(a.state, b.state)
end

# both orders spelled out: Julia would pick a catch-all `inner(::State, ::State)` over the
# diagonal method above, for matching pairs as well
"""
    different_representations()

refuse an inner product between a pure and a mixed representation
"""
different_representations() =
    error("no inner product between a pure and a mixed representation. Mix the pure " *
          "one, or use fidelity")

inner(::State{Pure}, ::State{Mixed}) = different_representations()
inner(::State{Mixed}, ::State{Pure}) = different_representations()

dot(a::State, b::State) = inner(a, b)

function inner(a::State{R}, op::IndexedOp{R}, b::State{R}) where R
    check_same_system(a, b)
    # on the operator as written, see `expect`
    check_indices(a.system, op)
    return inner(a.state', make_mpo(b, op), b.state)
end

dot(a::State, op::IndexedOp, b::State) = inner(a, op, b)
dot(a::AbstractState, b::AbstractState) = throw(MethodError(dot, (a, b)))

"""
    fidelity(a, b)

the fidelity of two states of the same system, between 0 and 1, whatever their norms and
traces: ``|\\langle a | b \\rangle|^2`` on two pure representations,
``\\langle \\psi | \\rho | \\psi \\rangle`` on a pure and a mixed one, in either order. Two
mixed representations are refused, their Uhlmann fidelity being out of reach for a matrix
product state: use `hs_fidelity`.

# Examples

    fidelity(state, ground_state)
    measurements = "data" => Fidelity(ground_state)
"""
fidelity(a::State{Pure}, b::State{Pure}) =
    abs2(inner(a, b)) / (norm(a)^2 * norm(b)^2)

fidelity(p::State{Pure}, r::State{Mixed}) =
    real(inner(mix(p), r)) / (norm(p)^2 * real(trace(r)))

fidelity(r::State{Mixed}, p::State{Pure}) = fidelity(p, r)

fidelity(::State{Mixed}, ::State{Mixed}) =
    error("no fidelity between two mixed representations, it needs the spectrum of a " *
          "density operator. Use hs_fidelity")

"""
    hs_fidelity(a::State{Mixed}, b::State{Mixed})

the normalised Hilbert-Schmidt overlap of two mixed states of the same system,
``\\mathrm{tr}(ab)/\\sqrt{\\mathrm{tr}(a^2)\\mathrm{tr}(b^2)}``, 1 exactly when they are
proportional. This is **not** the Uhlmann fidelity, see `fidelity`.
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
    # `conj`: the directions stay, only the charges are relabelled
    relab = relabeller(i -> adjoint_index(i, strong_names(s)))
    return State(state, MPS([ relabel(conj(state.state[i]), relab) * m
                              for (i, m) in enumerate(adj_maps(s, relab)) ]))
end


"""
    hermitianize(state [; limits])

the state whose density matrix is the Hermitian part ``(\\rho + \\rho^\\dagger)/2`` of that
of `state`, the sum being truncated according to `limits`, a `Limits` (default `Limits()`). A
pure representation is returned as it is.
"""
hermitianize(state::State{Pure}; limits::Limits=Limits()) =
    state
hermitianize(state::State{Mixed}; limits::Limits=Limits()) =
    State(state, 0.5*(+(state.state, dag(state).state;
                        limits.cutoff, limits.maxdim, limits.mindim)))


"""
    hermiticity(::State)

how Hermitian the density matrix is, from 0 when it is anti-Hermitian to 1 when it is
Hermitian, as it should be:
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
representation; given positions, or a cut standing for the sites `1:cut` (from 0 to the
number of sites), that of the state reduced to those sites, 0 when there is none. On a pure
state, positions mix it first, which is much more expensive, while a cut is as cheap as
`entanglement_entropy`.
"""
renyi2(::State{Pure}) = 0.
renyi2(state::State{Mixed}) = -log(trace2(state))

# the weak form, which `partial_trace` does not take on its own, is fine for a number
function renyi2(state::State{Mixed}, a::AbstractVector{Int})
    if isempty(a)
        return 0.0
    end
    return renyi2(partial_trace(weak_form(state), a; keep = true))
end

function renyi2(state::State{Pure}, a::AbstractVector{Int})
    if isempty(a)
        return 0.0
    end
    return renyi2(mix(state), a)
end

# positions in a vector of another element type, `[]` or `Any[1, 3]`
renyi2(state::State, a::AbstractVector) = renyi2(state, Vector{Int}(a))

"""
    unroll(x)

turn an array of results over the sites, or pairs of sites, each holding one value per
operator, into one such array per operator, arranged as the operators; an array of numbers is
returned as it is
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
the tensor of the state on `pos` included, its site left open for an operator.
`Expector()` is the empty one, at `pos = 0`.
"""
struct Expector
    pos::Int
    t::ITensor
end

Expector() =
    Expector(0, ITensor())


"""
    zipto(state, a, i)

the expector `a` carried to site `i`, at or after `a.pos`, see `zip_between`; from the empty
expector, `get_left(state, i)`
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

the contraction `t` of an expector at `from`, its operator placed there, carried to `to`, the
sites in between traced out and the site of `to` left open
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

the expector `a` completed with the environment on the right of `a.pos`, see `get_right`
"""
zipend(state::State, a::Expector) =
    Expector(a.pos, a.t * get_right(state, a.pos))


"""
    expectfactor(state, a, o)
    expectfactor(state, a, t, i)

the expector `a` carried to the site of the one site operator `o` and multiplied by the tensor
measuring it, or to site `i` and multiplied by `t`; a `Multi_F` places `F` on each of its sites
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
`obs_at`; on a pure state the identity only renames the index of the site
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

`expect` without the simplification nor the normalization, on an operator already in the
form `simplify` gives, or an array of them; given `coef` and `factors`, or `coef` and a com,
that of their product. A com is contracted with one expector per channel.
"""
function expect_norm(state::State, coef::Number, subs::Vector{<:IndexedOp{Pure}})
    if coef == 0.
        return 0.
    end
    state = weak_form(state)
    # every expectation value reaches this leaf, `measure` included
    foreach(o -> check_indices(state.system, o), subs)
    foreach(o -> check_one_site(o, "expect"), subs)
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
                # out of place: a piece may share its storage with the cache of `get_left`
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
array of them, divided by the real part of the trace of the state (its squared norm on a pure
representation). With `normalize = false` it is not divided: ``\\mathrm{tr}(A\\rho)``, or
``\\langle\\psi|A|\\psi\\rangle``, as needed by a non Hermitian ``B\\rho`` made with `Left` or
`Right` in a correlation at two times, whose trace may be complex or zero.

`obs` is simplified first: its factors may come in any order and the Jordan-Wigner strings
of fermionic operators are inserted for you. An operator not placed on sites, `X` rather
than `X(1)`, is refused (see `expect1`), and so is a superoperator.

A representation of one's own measures operators through `expect`, see
[Representations of one's own](@ref).

# Examples

    expect(state, X(1)*Y(2) + Y(1)*Z(3))
    expect(state, [X(1)*Y(2), X(3), Z(1)*X(2)])
    expect(state, C(3)*dag(C)(1))
"""
function expect(state::State, op::Op; normalize::Bool = true)
    # before simplify, which moves an identity to the first site and adds Jordan-Wigner
    # strings
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

# before the method below, which would iterate the operator
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
array of operators, `t` being the state contracted on every site with that of `i` left open;
a fermionic operator is refused, see `expect1`
"""
expect1_one(state::State, op::SimpleOp, i::Int, t::ITensor) =
    if isfermionic(op)
        # a refusal, not a gap to fill: its only correct answer is zero, and implementing it
        # would silently open the branch `expect2` picks on the parity of its first operator
        error("expect1 does not take a fermionic operator: its expectation value is odd, " *
              "so it vanishes on any state of definite fermion parity. If you really " *
              "want it on a state that superposes parities, ask for it site by site " *
              "with expect(state, op(i))")
    else
        scalar(t * obs_at(state, op, i))
    end

expect1_one(::State, a::IndexedOp, ::Int, ::ITensor) =
    error("expect1 measures an operator on every site, X rather than X(1): measure $a with expect")

expect1_one(state::State, ops, i::Int, t::ITensor) =
    map(ops) do o
        expect1_one(state, o, i, t)
    end


"""
    expect1(state, op)

the expectation values of the one site operator `op` on every site, as a vector indexed by
site, or, for an array of operators, an array of such vectors arranged as the operators,
normalised as by `expect`. A fermionic operator is refused, its expectation value vanishing
on any state of definite parity: on a state superposing parities, use `expect(state, op(i))`.

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


"""
    check_pair_parities(ops)

refuse a pair of operators of different fermionic parities, whose product is odd: `expect2`
picks the branch and the swap sign of a pair on the parity of its first operator alone
"""
function check_pair_parities(ops)
    for (o1, o2) in ops
        if isfermionic(o1) ≠ isfermionic(o2)
            error("cannot correlate $o1 and $o2: one is fermionic and the other is not, " *
                  "so their product is odd and has no expectation value. Both operators " *
                  "of a pair must have the same fermionic parity")
        end
    end
end

"""
    mirror(state, o1, o2)

the function giving the entry `(j, i)`, `j > i`, of the pair `(o1, o2)` from its entry `(i, j)`,
or `nothing`: for `o1 == o2`, the same value, of opposite sign for a fermionic operator, the
two sites being different; on a pure state, for `o2` the adjoint of `o1`, its conjugate,
`(o1(i) o2(j))†` being `o1(j) o2(i)`. On a mixed state, which truncation leaves only nearly
hermitian, the second is computed.
"""
function mirror(state::State, o1::SimpleOp, o2::SimpleOp)
    if o1 == o2
        return isfermionic(o1) ? (-) : identity
    elseif state isa State{Pure} && simplify(dag(o1)) == simplify(o2)
        return conj
    end
    return nothing
end

"""
    closed_row(state, e, a, i)

the operator `a` placed on site `i` of the expector `e`, with its string when it is fermionic,
and the result completed on each later site `j`, that site left open for a second operator
"""
function closed_row(state::State, e::Expector, a::SimpleOp, i::Int)
    ferm = isfermionic(a)
    env = expectfactor(state, e, obs_at(state, ferm ? a * F : a, i), i)
    n = length(state)
    row = Vector{ITensor}(undef, n)
    for j in i+1:n
        row[j] = zipend(state, zipto(state, env, j)).t
        if j < n
            next = ferm ? tensor_obs(state, F(j)) : identity_at(state, j)
            env = zipto(state, expectfactor(state, env, next, j), j + 1)
        end
    end
    return row
end

# on a pure state, the operator of site i is placed at once and the environment carried to the
# right is closed: left open, it carries the ket and the bra of site i, a factor of the dimension
# squared at every step, for one environment per operator. On a mixed state the environment is
# a vector, which an open site costs little, and the method below shares it between all pairs
function expect2(state::State{Pure}, ops::Vector{<:Tuple{SimpleOp, SimpleOp}})
    check_pair_parities(ops)
    state = weak_form(state)
    n = length(state)
    scale = real(trace(state))
    mirrors = [ mirror(state, o1, o2) for (o1, o2) in ops ]
    # the operators placed first, the second of a pair only for its entries (j, i)
    firsts = unique([first.(ops); [ o2 for ((_, o2), m) in zip(ops, mirrors) if isnothing(m) ]])
    r = Matrix{Any}(undef, n, n)
    for i in 1:n
        e = zipto(state, Expector(), i)
        e = Expector(e.pos, e.t / scale)
        t = zipend(state, e).t
        r[i, i] = map(((o1, o2),) -> expect1_one(state, o1 * o2, i, t), ops)
        rows = Dict(a => closed_row(state, e, a, i) for a in firsts)
        for j in i+1:n
            r[i, j] = map(((o1, o2),) -> scalar(rows[o1][j] * obs_at(state, o2, j)), ops)
            # swapping two fermionic operators costs a sign
            r[j, i] = map(zip(ops, mirrors, r[i, j])) do ((o1, o2), m, v)
                if !isnothing(m)
                    return m(v)
                end
                w = scalar(rows[o2][j] * obs_at(state, o1, j))
                return isfermionic(o1) ? -w : w
            end
        end
    end
    return unroll(r)
end

function expect2(state::State{Mixed}, ops::Vector{<:Tuple{SimpleOp, SimpleOp}})
    check_pair_parities(ops)
    state = weak_form(state)
    oplist = [first.(ops) ; last.(ops)]
    need_fermionic = any(isfermionic, oplist)
    need_non_fermionic = any(x->!isfermionic(x), oplist)
    n = length(state)
    mirrors = [ mirror(state, o1, o2) for (o1, o2) in ops ]
    # the element type follows what `scalar` returns; `unroll` gives a concretely typed result
    r = Matrix{Any}(undef, n, n)
    scale = real(trace(state))
    for i in 1:n
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
            # a and b of the same parity: a fermionic a takes its string through the F of
            # its site
            correlation(a, b) =
                isfermionic(a) ?
                    scalar(ef.t * tensor_obs(state, (a * F)(i)) * obs_at(state, b, j)) :
                    scalar(enf.t * obs_at(state, a, i) * obs_at(state, b, j))
            r[i, j] = map(((o1, o2),) -> correlation(o1, o2), ops)
            # swapping two fermionic operators costs a sign
            r[j, i] = map(zip(ops, mirrors, r[i, j])) do ((o1, o2), m, v)
                if !isnothing(m)
                    return m(v)
                end
                return isfermionic(o1) ? -correlation(o2, o1) : correlation(o2, o1)
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

the matrix whose entry `(i, j)` is the expectation value of `o1(i) * o2(j)`, normalised as by
`expect`, the diagonal being that of `(o1 * o2)(i)`; for a vector of pairs, a vector of such
matrices. The two operators of a pair must both be fermionic or both not; the Jordan-Wigner
strings are inserted for you.

# Examples

    expect2(state, (X, X))
    expect2(state, [(X, Y), (X, Z), (Y, Z)])
"""
expect2(state::State, ops::Tuple{SimpleOp, SimpleOp}) =
    expect2(state, [ops])[1]

expect2(::State, ops) =
    error("expect2 correlates two operators on every pair of sites, (X, Y) rather than " *
          "(X(1), Y(2)): measure X(1) * Y(2) with expect")

"""
    variance(state::State{Pure}, hamiltonian)
    variance(state::State{Pure}, ::MPO)

the variance of the energy, ``\\langle H^2 \\rangle - \\langle H \\rangle^2``, zero exactly
when the state is an eigenstate of the Hamiltonian: the convergence check of a ground state
search, one stuck in a metastable state also ceasing to progress, but with a large variance.
The state need not be normalised, and a mixed representation is refused. The Hamiltonian may
be given as its `make_mpo(state, hamiltonian)`, to build it once.

It costs about a `dmrg` sweep at the same bond dimension: it belongs in
`final_measurements`, or under a large `measurements_period`.

# Examples

    energy, gs = dmrg(hamiltonian, state; nsweeps = 10)
    variance(gs, hamiltonian)
    measurements = "data" => Variance(hamiltonian)
"""
function variance(state::State{Pure}, mpo::MPO)
    st = state.state
    n2 = real(dot(st, st))
    # the bra primed: unprimed, `inner` takes a deprecated fallback
    e = real(inner(prime(st), mpo, st)) / n2
    e2 = real(inner(mpo, st, mpo, st)) / n2
    return e2 - e^2
end

variance(state::State{Pure}, h) = variance(state, make_mpo(state, h))

variance(::State{Mixed}, _) =
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

the entropy ``-\\sum_i p_i \\log p_i`` of the probabilities `p`, a zero probability, which a
`mindim` above the Schmidt rank gives, adding nothing rather than NaN
"""
shannon_entropy(p) = -sum(x * log(x) for x in p if x > 0)

"""
    entropy_spectrum(S)

the entanglement entropy and spectrum read off the singular values `S` of a cut, as
`entanglement_entropy` returns them
"""
function entropy_spectrum(S::ITensor)
    # sorted: with conserved quantities the diagonal runs block by block
    sp = sort([ S[i,i]^2 for i in 1:dim(S, 1) ]; rev = true)
    sp /= sum(sp)
    return (shannon_entropy(sp), sp)
end

"""
    entanglement_entropy(state, cut::Int)

the entanglement entropy across the cut between sites `cut` and `cut + 1`, and the spectrum
it is computed from, the eigenvalues of the reduced density matrix of the sites up to `cut`
in decreasing order. `cut` runs from 1 to the number of sites, the last one giving 0. On a
mixed representation, this is the operator space entanglement entropy (OSEE).

# Examples

    ee, spectrum = entanglement_entropy(state, 3)   # cut between sites 3 and 4
"""
entanglement_entropy(state::State, cut::Int) = entropy_spectrum(cut_svd(state, cut)[2])

"""
    entanglement_by_sector(state::State{Pure}, cut::Int)

the entanglement across the cut between sites `cut` and `cut + 1`, resolved by the charge of
the sites up to `cut`: a `Dict` from each charge, a `QN` (`using ITensors: QN`), to
`(weight, entropy, spectrum)`, the probability of that charge and the entropy and decreasing
spectrum of the reduced density matrix restricted to it and normalised. The entanglement
entropy is `Σ weight * entropy - Σ weight * log(weight)`, the second sum being the number
entropy. A state conserving nothing has the single sector `QN()`; a mixed representation is
refused.

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
    # U carrying no flux, what flows through `u` is what the sites up to `cut` hold
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
        # a sector kept with nothing in it
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
fermionic signs included: ``U \\rho U^\\dagger``, ``U = (-1)^{\\sum n_l n_k}``, `l` running
over the fermionic sites traced out, `k` over the fermionic sites kept on their right, and `n`
the parity of a site. It is the product by an MPO of bond dimension 2 carrying the parity
traced out on the left, projected on the ket side alone for a site traced out. A state
needing no sign is given back as it is.
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
        a = site_array(Float64, idx, l, r)
        id = Matrix{Float64}(I, dim(idx), dim(idx))
        for p in 1:2
            if !odd[i]
                add_block!(a, p, p, id)
            elseif i in keep
                add_block!(a, p, p, p == 1 ? id : site_matrix(tensor(sys, Gate(F)(i)), idx))
            else
                for q in 1:2
                    add_block!(a, p, xor(p - 1, q - 1) + 1, site_matrix(tensor(sys, projs[q](i)), idx))
                end
            end
        end
        return noprime(site_tensor(a, idx, l, r) * mps[i])
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
    partial_trace(state, positions::AbstractVector{Int} [; keep = false])

the state with the sites at `positions` traced out, or, with `keep = true`, all the others: a
mixed state on a new system made of the sites kept, in their order, with the trace of
`state`. An operator of the sites kept has on it the expectation value it has on `state`,
fermionic signs included; tracing out a fermionic site with fermionic sites kept on its
right may double the bond dimension.

A pure representation is refused (use `mix(state)` first), so is a state conserving
something strongly (weaken it first), and so is tracing out every site.

# Examples

    partial_trace(state, [1, 2])
    partial_trace(state, 3:5; keep = true)
"""
function partial_trace(state::State{Mixed}, pos::AbstractVector{Int}; keep::Bool = false)
    if !isempty(strong_names(state.system))
        error("cannot trace out part of a state conserving something strongly, what is left " *
              "spreads over several sectors: weaken it first, giving up the strong symmetry")
    end
    n = length(state)
    # the filter below would silently ignore a position the state does not have
    check_positions(state, pos, "partial_trace")
    if keep
        kept_sites = sort(unique(pos))
    else
        kept_sites = filter(e->e ∉ pos, 1:n)
    end
    kn = length(kept_sites)
    if kn == 0
        error("partial_trace cannot trace all sites of a state")
    end
    state = trace_signs(state, kept_sites)
    mps = state.state
    sys = state.system
    j = 0
    t = Vector{ITensor}(undef, kn)
    for (i, k) in enumerate(kept_sites)
        if i == 1
            # not divided by the trace, which the result keeps and which may be zero
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
    t[kn] *= get_right(state, kept_sites[kn])
    for i in 1:kn-1
        idx = commonind(t[i], t[i+1])
        jdx = settags(idx, "Link, l=$i")
        replaceind!(t[i], idx, jdx)
        replaceind!(t[i+1], idx, jdx)
    end
    kept = System(AbstractSite[ sys[k] for k in kept_sites ],
                  Index[ SysIndex{Pure}(sys, k) for k in kept_sites ],
                  Index[ SysIndex{Mixed}(sys, k) for k in kept_sites ])
    return State{Mixed}(kept, MPS(t))
end

# not mixed behind the caller's back as in `renyi2`, the result being a state
partial_trace(::State{Pure}, ::AbstractVector{Int}; kwargs...) =
    error("partial_trace needs a mixed representation, use mix(state) first")

"""
    check_cut(state, cut, what)

refuse a cut the state does not have, `what` naming the function given it, `cut` counting the
sites on its left, from 0 to the number of sites
"""
function check_cut(state::State, cut::Int, what)
    n = length(state)
    if !(0 ≤ cut ≤ n)
        error("$what was given the cut $cut, which the state does not have: " *
              "a cut counts the sites on its left, from 0 to $n")
    end
    return nothing
end

function renyi2(state::State, cut::Int)
    check_cut(state, cut, "renyi2")
    if cut == 0
        return 0.0
    elseif cut == length(state)
        return renyi2(state)
    end
    return renyi2(state, collect(1:cut))
end

# read off the Schmidt spectrum rather than from a partial trace of a mix
function renyi2(state::State{Pure}, cut::Int)
    check_cut(state, cut, "renyi2")
    if cut == 0 || cut == length(state)
        return 0.0
    end
    return -log(sum(abs2, last(entanglement_entropy(state, cut))))
end

"""
    parts_renyi2(state, a)

``S_2(A) + S_2(B) - S_2(A \\cup B)`` for part A at the positions `a` and B the rest, neither
of them empty; ``2 S_2(A)`` on a pure state
"""
function parts_renyi2(state::State{Mixed}, a::AbstractVector{Int})
    w = weak_form(state)
    return renyi2(partial_trace(w, a; keep = true)) +
           renyi2(partial_trace(w, a; keep = false)) -
           renyi2(w)
end

parts_renyi2(state::State{Pure}, a::AbstractVector{Int}) = 2 * renyi2(state, a)

"""
    mutual_info_renyi2(state::State, cut::Int)
    mutual_info_renyi2(state::State, a::AbstractVector{Int})

the Rényi-2 analogue of the mutual information between two parts of the state,
``S_2(A) + S_2(B) - S_2(A \\cup B)``, which can be negative on a mixed state. Part A is given
by its positions, or by a cut, sites `1:cut`, from 0 to the number of sites; part B is the
rest. An empty part gives 0.

On a pure state it is twice the `renyi2` of part A, with the same costs, the result holding
for fermions on a state of definite parity, see `RandomState`.
"""
function mutual_info_renyi2(state::State, a::AbstractVector{Int})
    n = length(state)
    check_positions(state, a, "mutual_info_renyi2")
    # partial_trace would have no site to keep
    k = length(unique(a))
    if k == 0 || k == n
        return 0.0
    end
    return parts_renyi2(state, a)
end

function mutual_info_renyi2(state::State, cut::Int)
    check_cut(state, cut, "mutual_info_renyi2")
    return mutual_info_renyi2(state, collect(1:cut))
end

# the cut checked here for the message to name this function
function mutual_info_renyi2(state::State{Pure}, cut::Int)
    check_cut(state, cut, "mutual_info_renyi2")
    return 2 * renyi2(state, cut)
end

# positions in a vector of another element type
mutual_info_renyi2(state::State, a::AbstractVector) = mutual_info_renyi2(state, Vector{Int}(a))



"""
    draw(p, d, rnd)

the first outcome among `0:d-1` whose cumulated probability exceeds `rnd`, `p(x)` being the
probability of `x`, the last one taking what the others leave
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

a random outcome, numbered from 0, of measuring the state in the computational basis: a vector
with one outcome per site, or, given `pos`, the outcome of that site alone. `rng` is the random
number generator (default the global one).

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
    # the conjugate leg, the direction the tensor takes on a charged site
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
        # rescaled, the probabilities to come being ratios: it would underflow on long chains
        l /= norm(l)
    end
    return result
end