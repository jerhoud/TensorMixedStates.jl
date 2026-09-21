export trace, trace2, norm, normalize, hermitianize, hermiticity, renyi2
export inner, dot, fidelity, hs_fidelity
export expect, expect1, expect2
export entanglement_entropy, partial_trace, mutual_info_renyi2, sample, variance

function tensor_trace(state::State{Mixed}, i::Int)
    s = state.system
    j = SysIndex{Pure}(s, i)
    k = SysIndex{Mixed}(s, i)
    return dense(delta(j', j)) * combinerto(k, j', j)
end

tensor_obs(state::State{Pure}, ind::AtIndex{Pure, 1}) =
    tensor(state.system, ind)

function tensor_obs(state::State{Mixed}, ind::AtIndex{Pure, 1})
    s = state.system
    t = tensor(s, ind)
    j = SysIndex{Pure}(s, ind.index...)
    k = SysIndex{Mixed}(s, ind.index...)
    return t * combinerto(k, j, j')
end

# `(c * A)(i)` keeps its coefficient outside the AtIndex, so it has to be taken off here:
# everything downstream of `tensor_obs` works on the one site tensor alone
tensor_obs(state::State, a::ScalarOp{Pure, Indexed, 1}) =
    a.coef * tensor_obs(state, a.arg)

tensor_obs(state::State{Mixed}, i::Int) =
    tensor_trace(state, i) * state.state[i]

function tensor_dag(state::State, i::Int)
    s = state.system
    j = SysIndex{Pure}(s, i)
    k = SysIndex{Mixed}(s, i)
    return dag(state.state[i]) * combinerto(k, j', j) * combinerto(k, j, j')
end

function get_loc(state::State, i::Int)
    l = state.preobs.loc
    if isempty(l)
        create_loc!(l, state)
    end
    return l[i]
end

function get_right(state::State, i::Int)
    r = state.preobs.right
    if isempty(r)
        create_right!(r, state)
    end
    return r[i]
end

function get_left(state::State, i::Int)
    l = state.preobs.left
    if length(l) < i
        create_left!(l, state, i)
    end
    return l[i]
end

"""
    trace(::State)

Return the trace of the system, mostly useful for mixed representations.
This should be one.
"""
function trace(state::State)
    t = state.preobs.trace
    if isempty(t)
        create_trace!(t, state)
    end
    return t[1]
end

create_loc!(_, ::State{Pure}) = error("bug: get_loc on pure states")
function create_loc!(l, state::State{Mixed})
    n = length(state)
    resize!(l, n)
    for i in 1:n
        l[i] = tensor_obs(state, i)
    end
    return l
end

function create_right!(r, state::State{Pure})
    st = state.state
    s = state.system
    n = length(state)
    rl = ITensorMPS.rightlim(st) - 1
    resize!(r, n)
    r[n] = dag(st[n]')
    for i in n-1:-1:1
        rlink = commonind(st[i], st[i+1])
        v = if i >= rl
            delta(rlink, rlink')
        else
            k = SysIndex{Pure}(s, i+1)
            r[i+1] * delta(k, k') * st[i+1]
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

function create_trace!(t, state::State{Pure})
    resize!(t, 1)
    k = SysIndex{Pure}(state.system, 1)
    t[1] = scalar(get_right(state, 1) * delta(k, k') * state.state[1])
    return t
end

function create_trace!(t, state::State{Mixed})
    resize!(t, 1)
    t[1] = scalar(get_loc(state, 1) * get_right(state, 1))
    return t
end

function create_left!(l, state::State{Pure}, i::Int)
    st = state.state
    s = state.system
    ll = ITensorMPS.leftlim(st) + 1
    j = length(l)
    resize!(l, i)
    if j == 0
        l[1] = st[1] / real(trace(state))
        j = 1
    end
    for k in j+1:i
        llink = commonind(st[k-1], st[k])
        v = if k <= ll
            delta(llink, llink') / real(trace(state))
        else
            idx = SysIndex{Pure}(s, k-1)
            l[k-1] * delta(idx, idx') * dag(st[k-1]')
        end
        l[k] = v * st[k]
    end
    return l
end

function create_left!(l, state::State{Mixed}, i::Int)
    st = state.state
    j = length(l)
    resize!(l, i)
    if j == 0
        l[1] = st[1] / real(trace(state))
        j = 1
    end
    for k in j+1:i
        l[k] = l[k-1] * tensor_trace(state, k-1) * st[k]
    end
    return l
end

"""
    trace2(::State)

Return the trace of the square density matrix, mostly useful for mixed representations.
This is one for pure representations.
"""
trace2(::State{Pure}) = 1.
trace2(state::State{Mixed}) = (norm(state.state) / real(trace(state))) ^ 2

"""
    norm(::State)

Return the norm of the state, mostly useful for pure representations.
This should be one for pure representation.
"""
norm(state::State) = norm(state.state)


"""
    normalize(::State)

normalize state so that norm = 1 for pure state and trace = 1 for mixed state
"""
normalize(state::State{Pure}) =
    State(state, normalize(state.state))
normalize(state::State{Mixed}) =
    State(state, state.state / real(trace(state)))

"""
    check_same_system(a, b)

refuse two states that cannot be contracted together. The indices of a `System` are drawn
afresh, so two states of two systems have nothing in common even when their sites match,
and `ITensorMPS` would contract them anyway, on a deprecated fallback that matches by
position and warns.
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

the inner product of two states of the same system, `dot` being an alias of `inner`.

On pure representations this is the overlap ``\\langle a | b \\rangle``, conjugating the
first argument. On mixed ones it is the Hilbert-Schmidt product ``\\mathrm{tr}(a^\\dagger b)``,
the two density matrices being contracted as the vectors they are stored as. Neither is
normalised: divide by the norms, or use `fidelity`.

The two states must share their `System`, see `State(::System, ::State)`.

# Examples

    inner(state, ground_state)
    abs2(inner(a, b))              # the Loschmidt echo of a pure state
"""
function inner(a::State{R}, b::State{R}) where R
    check_same_system(a, b)
    return dot(a.state, b.state)
end

# A pure state is a vector of the Hilbert space and a mixed one a vector of the space of
# operators on it, so there is no product of the two to take. Said here rather than left to
# the method above, which would only report that none matches. Both orders are spelled out
# on purpose: a catch-all `inner(::State, ::State)` is what Julia picks over the diagonal
# method above, so it would capture the matching pairs as well.
different_representations() =
    error("no inner product between a pure and a mixed representation. Mix the pure " *
          "one, or use fidelity")

inner(::State{Pure}, ::State{Mixed}) = different_representations()
inner(::State{Mixed}, ::State{Pure}) = different_representations()

dot(a::State, b::State) = inner(a, b)

"""
    fidelity(a, b)

the fidelity of two states of the same system, a number between 0 and 1, normalised so
that the norm and the trace of its arguments do not matter.

On two pure representations this is ``|\\langle a | b \\rangle|^2``. On a pure and a mixed
one, in either order, it is ``\\langle \\psi | \\rho | \\psi \\rangle``, which is the
fidelity of a mixed state with a pure target and costs no more than an overlap.

Two mixed representations are refused, with a message saying what to reach for instead:
the Uhlmann fidelity ``(\\mathrm{tr}\\sqrt{\\sqrt{\\rho}\\sigma\\sqrt{\\rho}})^2`` needs the
square root of a density operator, hence its spectrum, which is out of reach for a matrix
product state. See `hs_fidelity` for what can be computed.

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

the normalised Hilbert-Schmidt overlap of two mixed states of the same system, that is
``\\mathrm{tr}(ab)/\\sqrt{\\mathrm{tr}(a^2)\\mathrm{tr}(b^2)}``. It is 1 exactly when the two
density matrices are proportional, and the normalisation by the traces cancels out, which
leaves it the cosine between the two states seen as vectors.

This is **not** the Uhlmann fidelity, which is out of reach for a matrix product state, see
`fidelity`. It is a cheaper indicator of how close two mixed states are, and it is what to
reach for when comparing an evolution with a reference density matrix.
"""
hs_fidelity(a::State{Mixed}, b::State{Mixed}) =
    real(inner(a, b)) / (norm(a) * norm(b))

"""
    dag(::State)

adjoint of density matrix for mixed representation
"""
dag(state::State{Pure}) =
    error("dag is meaningless on pure representations")
function dag(state::State{Mixed})
    n = length(state)
    st = MPS(n)
    for i in 1:n
        st[i] = tensor_dag(state, i)
    end
    State(state, st)
end


"""
    hermitianize(::State)

modify the state so that it is Hermitian (only useful for mixed state)
"""
hermitianize(state::State{Pure}; kwargs...) =
    state
hermitianize(state::State{Mixed}; limits::Limits=Limits()) =
    State(state, 0.5*(+(state.state, dag(state).state; limits.cutoff, limits.maxdim)))


"""
    hermiticity(::State)

hermiticity measure whether density matrix for mixed state is Hermitian as it should.

return a value from 0 (anti Hermitian) to 1 (Hermitian)

return 1 for pure state
"""
hermiticity(::State{Pure}) = 1.
hermiticity(state::State{Mixed}) =
    0.5 + 0.5 * real(dot(state.state, dag(state).state)) / norm(state.state)^2

"""
    renyi2(::State)
    renyi2(::State, ::AbstractVector{<:Integer})

renyi2 returns the Renyi entropy of order 2 of the state. It is 0. for a pure
representation, whose density matrix has a single non zero eigenvalue.

When given an array of positions it returns the Renyi entropy of the corresponding
substate, that is `-log(tr(ρ²))` of the state reduced to those sites. On a pure state this
is a measure of how much those sites are entangled with the rest, and is not 0.: the state
is first turned into its mixed representation, which is much more expensive, because a
partial trace needs a density matrix.
"""
renyi2(::State{Pure}) = 0.
renyi2(state::State{Mixed}) = -log(trace2(state))

renyi2(state::State{Mixed}, a::AbstractVector{<:Integer}) =
    renyi2(partial_trace(state, a; keepers = true))

# a subsystem of a pure state is not pure, so this is an entanglement measure rather than
# 0. There is no cheap route for an arbitrary subset, the same way mutual_info_renyi2 has
# none: a partial trace needs a density matrix.
renyi2(state::State{Pure}, a::AbstractVector{<:Integer}) =
    renyi2(mix(state), a)

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


struct Expector
    pos::Int
    t::ITensor
end

Expector() =
    Expector(0, ITensor())


zipto(state::State{Pure}, a::Expector, i::Int) = 
    if a.pos == 0
        Expector(i, get_left(state, i))
    elseif a.pos == i
        a
    else
        st = state.state
        t = a.t * dag(st[a.pos]')
        for k in a.pos+1:i-1
            idk = SysIndex{Pure}(state.system, k)
            t *= st[k] * delta(idk, idk')
            t *= dag(st[k]')
        end
        t *= st[i]
        Expector(i, t)
    end

zipto(state::State{Mixed}, a::Expector, i::Int) =
    if a.pos == 0
        Expector(i, get_left(state, i))
    elseif a.pos == i
        a
    else
        t = a.t
        for k in a.pos+1:i-1
            t *= get_loc(state, k)
        end
        t *= state.state[i]
        Expector(i, t)
    end

zipend(state::State, a::Expector) =
    Expector(a.pos, a.t * get_right(state, a.pos))


function expectfactor(state::State, a::Expector, o::AtIndex)
    a = zipto(state, a, o.index...)
    Expector(a.pos, a.t * tensor_obs(state, o))
end

function expectfactor(state::State, a::Expector, o::Multi_F)
    for i in o.start:o.stop
        a = expectfactor(state, a, F(i))
    end
    return a
end


"""
    expect_norm(::State, obs)

expectation values of an operator already in the form `simplify` produces. This is what
`measure` uses, its operators having been normalised once and for all by `make_obs`.
"""
function expect_norm(state::State, coef::Number, subs::Vector{<:IndexedOp{Pure}})
    if coef == 0.
        return 0.
    end
    # every expectation value reaches this leaf, `measure` included, which calls
    # `expect_norm` rather than `expect`. Checking here covers them all at the cost of a
    # few integer comparisons per term, nothing next to the contractions below
    foreach(o -> check_indices(state.system, o), subs)
    e = Expector()
    for o in subs
        e = expectfactor(state, e, o)
    end
    e = zipend(state, e)
    return coef * scalar(e.t)
end

"""
    expect(::State, obs)

Compute expectation values of `obs` on the given state.

`obs` is simplified first, so its factors may be given in any order and the Jordan-Wigner
strings of fermionic operators are inserted for you.

# Examples
    expect(state, X(1)*Y(2) + Y(1)*Z(3))
    expect(state, [X(1)*Y(2), X(3), Z(1)*X(2)])
    expect(state, C(3)*dag(C)(1))

"""
expect(state::State, op) = expect_norm(state, simplify(op))

expect_norm(state::State, p::IndexedOp{Pure}) =
    expect_norm(state, scalarcoef(p), prodsubs(p))

expect_norm(state::State, op::SumOp{Pure, Indexed}) =
    sum(op.subs) do p
        expect_norm(state, p)
    end

expect_norm(state::State, op) =
    map(op) do o
        expect_norm(state, o)
    end

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
        scalar(t * tensor_obs(state, op(i)))
    end

expect1_one(state::State, ops, i::Int, t::ITensor) =
    map(ops) do o
        expect1_one(state, o, i, t)
    end


"""
    expect1(::State, op)

Compute the expectation values of the given operators on all sites.

# Examples
    expect1(state, X)
    expect1(state, [X, Y, Z])

"""
function expect1(state::State, op)
    n = length(state)
    r = [ expect1_one(state, op, i, zipend(state, zipto(state, Expector(), i)).t) for i in 1:n ]
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
    oplist = [first.(ops) ; last.(ops)]
    need_fermionic = any(isfermionic, oplist)
    need_non_fermionic = any(x->!isfermionic(x), oplist)
    n = length(state)
    # `Matrix{Any}` is said out loud rather than left to the bare `Matrix(undef, ...)`,
    # which means the same thing without showing it. The element type is not known before a
    # cell is computed, since it follows what `scalar` returns for this state, and `unroll`
    # rebuilds a concretely typed result at the end, so the untyped container stays internal
    r = Matrix{Any}(undef, n, n)
    for i in 1:n
        t = zipend(state, zipto(state, Expector(), i)).t
        r[i, i] = map(ops) do (o1, o2)
            expect1_one(state, o1 * o2, i, t)
        end
        lnf = zipto(state, Expector(), i)
        lf = lnf
        for j in i+1:n
            enf = lnf
            ef = lf
            if need_non_fermionic
                enf = zipend(state, zipto(state, lnf, j))
            end
            if need_fermionic
                ef = zipend(state, zipto(state, lf, j))
            end
            r[i, j] = map(ops) do (o1, o2)
                if isfermionic(o1)
                    scalar(ef.t * tensor_obs(state, (o1 * F)(i)) * tensor_obs(state, o2(j)))
                else
                    scalar(enf.t * tensor_obs(state, o1(i)) * tensor_obs(state, o2(j)))
                end
            end
            r[j, i] = map(ops) do (o1, o2)
                if isfermionic(o1)
                    # swapping the two fermionic operators costs a sign
                    -scalar(ef.t * tensor_obs(state, (o2 * F)(i)) * tensor_obs(state, o1(j)))
                else
                    scalar(enf.t * tensor_obs(state, o2(i)) * tensor_obs(state, o1(j)))
                end
            end
            if j < n
                if need_non_fermionic
                    lnf = zipto(state, expectfactor(state, lnf, Id(j)), j+1)
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
    expect2(::State, op_pairs)

Compute the 2-point correlations of the given pairs of operators on all sites.

# Examples
    expect2(state, (X, X))
    expect2(state, [(X, Y), (X, Z), (Y, Z)])
    ...
"""
expect2(state::State, ops::Tuple{SimpleOp, SimpleOp}) =
    expect2(state, [ops])[1]

"""
    variance(hamiltonian, ::State{Pure})
    variance(::MPO, ::State{Pure})

the variance of the energy, ``\\langle H^2 \\rangle - \\langle H \\rangle^2``, which is
zero exactly when the state is an eigenstate of the hamiltonian and nowhere else.

This is the convergence check of a ground state search. The `tolerance` of `GroundState`
stops on the progress of the energy between two sweeps, which says that the optimisation
has stopped moving, not that it has arrived: a search stuck in a metastable state has no
progress left and a large variance. The variance also gives the error bar, by running
several bond dimensions and extrapolating the energy to zero variance, the two being
asymptotically linear in one another.

``H^2`` is never formed. TMS gives every term of a sum a channel of its own and does not
compress, so squaring a hamiltonian squares its number of terms and the bond dimension of
its MPO with them: on a transverse field Ising chain of forty sites, 3 becomes 1603. What
is computed instead is ``\\langle H\\psi | H\\psi \\rangle``, one contraction with the
MPO of `H` on either side, which is linear in the number of sites and does not depend on
the number of terms. It was measured 195 times faster than `expect(state, H * H)` there,
for the same value.

Pass an `MPO` to reuse one already built. The cost is that of a `dmrg` sweep at the same
bond dimension, so this belongs in `final_measures`, or under a large `measures_period`,
rather than at every sweep.

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
    entanglement_entropy(::State, ::Int)

Return the entanglement entropy of the given state across the cut on the right of the
given site, that is between `pos` and `pos + 1`, together with the spectrum it is computed
from.

`pos` runs from 1 to the number of sites. The last one cuts the whole state from nothing,
so it is always 0 and the cuts that say something are 1 to n-1.

The spectrum returned is the squared singular values of that cut, normalized to sum to
one, which are the eigenvalues of the reduced density matrix of the sites up to `pos`.

On a mixed representation the same quantity is computed on the vectorized density matrix,
which makes it the operator space entanglement entropy (OSEE) rather than an entanglement.

# Examples
    ee, spectrum = entanglement_entropy(state, 3)   # cut between sites 3 and 4
"""
function entanglement_entropy(state::State, pos::Int)
    n = length(state)
    if !(1 ≤ pos ≤ n)
        error("cannot compute the entanglement entropy at site $pos of a $n site state, " *
              "the cut is on the right of a site so pos must be between 1 and $n")
    end
    s = orthogonalize(state.state, pos)
    _, S = svd(s[pos], (linkinds(s, pos-1)..., siteinds(s, pos)...))
    sp = [ S[i,i]^2 for i in 1:dim(S, 1) ]
    sp /= sum(sp)
    ee = -sum(p * log(p) for p in sp)
    return (ee, sp)
end

"""
    partial_trace(::State, ::AbstractVector{<:Integer} [; keepers = false])

return the state partially traced at the given positions
alternatively one can give the positions to keep by setting `keepers = true`
"""
function partial_trace(state::State{Mixed}, pos::AbstractVector{<:Integer}; keepers::Bool = false)
    n = length(state)
    if keepers
        keep = sort(unique(pos))
    else
        keep = filter(e->e ∉ pos, 1:n)
    end
    kn = length(keep)
    if kn == 0
        error("partial_trace cannot trace all sites of a state")
    end
    mps = state.state
    sys = state.system
    j = 0
    t = Vector{ITensor}(undef, kn)
    s = Vector{AbstractSite}(undef, kn)
    sp = Vector{Index}(undef, kn)
    sm = Vector{Index}(undef, kn)
    for (i, k) in enumerate(keep)
        s[i] = sys[k]
        sp[i] = SysIndex{Pure}(sys, k)
        sm[i] = SysIndex{Mixed}(sys, k)
        if i == 1
            t[1] = copy(get_left(state, k))
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
    return State{Mixed}(System(s, sp, sm), MPS(t))
end

# a partial trace needs a density matrix: the reduced state of a subsystem is mixed in
# general, so there is nothing to hand back in pure representation. `renyi2` and
# `mutual_info_renyi2` mix on their own because they return a number; this one returns a
# state, and changing its representation behind the caller's back would be a surprise
partial_trace(::State{Pure}, ::AbstractVector{<:Integer}; kwargs...) =
    error("partial_trace needs a mixed representation, use mix(state) first")

"""
    mutual_info_renyi2(state::State, cut::Int)
    mutual_info_renyi2(state::State, a::AbstractVector{<:Integer})

return an approximation of the mutual information using renyi2 entropy.
You define the two parts either by giving the position of the cut between the left and right parts or by giving the list of positions for one of the parts.

On a pure state and for a cut, this is read directly from the entanglement spectrum and
costs nothing more than the entanglement entropy. For a list of positions the pure state
is first turned into its mixed representation, which is much more expensive.
"""
mutual_info_renyi2(state::State, a::AbstractVector{<:Integer}) =
    renyi2(partial_trace(state, a; keepers = true)) +
    renyi2(partial_trace(state, a; keepers = false)) -
    renyi2(state)

mutual_info_renyi2(state::State, cut::Int) =
    mutual_info_renyi2(state, collect(1:cut))

# a pure state has no entropy of its own, and the two sides of a cut share their Schmidt
# spectrum, so the mutual information is just twice the renyi2 entropy of either side
mutual_info_renyi2(state::State{Pure}, cut::Int) =
    -2 * log(sum(abs2, last(entanglement_entropy(state, cut))))

# partial_trace needs a density matrix, there is no cheap route for an arbitrary subset
mutual_info_renyi2(state::State{Pure}, a::AbstractVector{<:Integer}) =
    mutual_info_renyi2(mix(state), a)



"""
    sample(::State [; rng])
    sample(::State, pos::Int [; rng])

randomly sample the state in the computational basis, the outcomes are numbered from 0

Without a position, a complete configuration is sampled (one outcome per site), taking
the correlations between sites into account. With a position, only that site is sampled,
the other sites being traced out.

`rng` is the random number generator to use, it defaults to the global one.

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
    ptot = 0.
    ind = siteind(st, pos)
    d = dim(ind)
    for i in 1:d - 1
        ti = t * onehot(ind => i)
        ptot += real(scalar(ti * dag(ti)))
        if r < ptot
            return i - 1
        end
    end
    return d - 1
end

function sample(state::State{Mixed}, pos::Int; rng = Random.default_rng())
    sys = state.system
    l = get_left(state, pos)
    r = get_right(state, pos)
    d = dim(SysIndex{Pure}(sys, pos))
    rnd = rand(rng)
    ptot = 0.
    for x in 0:d - 2
        ptot += real(scalar(l * tensor_obs(state, Proj(x)(pos)) * r))
        if rnd < ptot
            return x
        end
    end
    return d - 1
end

function sample(state::State{Mixed}; rng = Random.default_rng())
    sys = state.system
    n = length(state)
    result = Vector{Int}(undef, n)
    l = ITensor(1.)
    for pos in 1:n
        a = l * state.state[pos]
        r = get_right(state, pos)
        d = dim(SysIndex{Pure}(sys, pos))
        tot = real(scalar(a * tensor_trace(state, pos) * r))
        rnd = rand(rng) * tot
        ptot = 0.
        x = d - 1
        al = a * tensor_obs(state, Proj(x)(pos))
        for i in 0:d - 2
            candidate = a * tensor_obs(state, Proj(i)(pos))
            ptot += real(scalar(candidate * r))
            if rnd < ptot
                x = i
                al = candidate
                break
            end
        end
        result[pos] = x
        l = al
    end
    return result
end