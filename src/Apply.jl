# apply on pure states: a product of operators applied factor by factor, right to left, never
# developed, each piece as a dense block on its sites, as one-site gates, or through its MPO. A
# channel of a single Kraus operator is applied as that operator. On mixed states, the factors
# of a channel likewise: a gate as its operator on the kets then on the bras, through its MPO
# made to act on the mixed indices, any other channel or a generator as a dense block on the
# mixed indices or through its MPO. A gate never brings sites together: on sites apart, it goes
# through its MPO on the sites between.

export apply

"""
    dense_limit

the largest product of the dimensions of the sites of a piece of more than two sites that
`apply` takes as a dense block; beyond, the piece goes through its MPO
"""
const dense_limit = 1024

"""
    check_apply_algo(apply_algo)

refuse an algorithm of the product of an MPO by a state other than `"densitymatrix"` and
`"naive"`: `"fit"` needs sweeps of its own, and `"zipup"` a later ITensorMPS than TMS requires
"""
function check_apply_algo(apply_algo::String)
    if apply_algo ∉ ("densitymatrix", "naive")
        error("apply_algo is \"densitymatrix\" or \"naive\", not $(repr(apply_algo))")
    end
end

"""
    piece_matrix(system, op, sites)

the matrix of the placed operator `op` on `sites`, in that order, the Jordan-Wigner strings
between them included, `op` acting on no other site: its products and sums taken as matrices,
nothing developed
"""
piece_matrix(system::System, a::AtIndex, sites) = matrix_on(system, a, sites)
piece_matrix(system::System, a::ProdOp, sites) = prod(piece_matrix(system, x, sites) for x in a.subs)
piece_matrix(system::System, a::SumOp, sites) = sum(piece_matrix(system, x, sites) for x in a.subs)
piece_matrix(system::System, a::ScalarOp, sites) = a.coef * piece_matrix(system, a.arg, sites)

function piece_matrix(system::System, ::IdentityOp, sites)
    d = prod(dim(system[k]) for k in sites)
    return Matrix{Float64}(I, d, d)
end

"""
    stays_on(system, m, sites)

whether the operator of matrix `m` on `sites`, increasing, lays no Jordan-Wigner string on any
other site: it is even, or no fermionic site is on its left, and it commutes with the parity of
its sites after each fermionic site between them that is not one of them
"""
function stays_on(system::System, m::AbstractMatrix, sites)
    fs = [ matrix(F, system[k]) for k in sites ]
    if next_fermion(system, 1) < first(sites) && parity_of(m, foldl(kron, fs)) ≠ 0
        return false
    end
    t = next_fermion(system, first(sites) + 1)
    while t < last(sites)
        if !(t in sites)
            p = foldl(kron, [ k > t ? f : one(f) for (k, f) in zip(sites, fs) ])
            if !nearly(p * m * p, m)
                return false
            end
        end
        t = next_fermion(system, t + 1)
    end
    return true
end

"""
    gate_on(m, idx)

the tensor of the gate of matrix `m` acting on the indices `idx`, the last varying fastest
"""
function gate_on(m::AbstractMatrix, idx)
    a, inds = on_legs(Matrix(m), [ i' for i in idx ], [ dag(i) for i in idx ])
    return charged_itensor(a, inds)
end

"""
    term_blocks(t)

the blocks of the term `t`, its coefficient left out, as pairs of a matrix and of its sites
"""
term_blocks(t::Term) = [ (m, [k]) for (k, m) in t.blocks ]

"""
    pieces(system, op)

the factors of the product `op`, gathered, see `gathered`, cut into the pieces `apply` applies
one by one, and the number they are multiplied by: the factors of one site of definite parity
are brought together where they commute with the pieces they cross, a sign for two odd ones,
so that their strings compensate, every other factor being a piece of its own
"""
function pieces(system::System, a::IndexedOp{IsOperator})
    g = gathered(system, a)
    coef = scalarcoef(g)
    ps = Vector{IndexedOp{IsOperator}}[]
    runs = Bool[]
    # the last run of factors of one site, and the pieces after it, each of one factor: their
    # sites, how many are odd, and whether one has no definite parity
    run = 0
    after, odd, none = Set{Int}(), 0, false
    for x in prodsubs(scalararg(g))
        if x isa IdentityOp
            continue
        end
        sx, px = sites_of(x), parity_on(system, x)
        joins = x isa AtIndex && length(x.index) == 1 && !isnothing(px)
        # it reaches the last run when it (anti)commutes with every piece after it, see
        # `sort_factors!`
        if joins && run > 0 && isdisjoint(sx, after) && (px == 0 || !none)
            if px == 1 && isodd(odd)
                coef = -coef
            end
            push!(ps[run], x)
            continue
        end
        push!(ps, [x])
        push!(runs, joins)
        if joins
            run = length(ps)
            after, odd, none = Set{Int}(), 0, false
        else
            union!(after, sx)
            odd += px == 1
            none |= isnothing(px)
        end
    end
    return coef, [ (ProdOp(p), r) for (p, r) in zip(ps, runs) ]
end

"""
    piece_action(system, p, run)

how `apply` applies the piece `p`: `(blocks, number)`, dense blocks to apply as gates, pairs of
a matrix and of its sites, and the number they bring, or `(terms, 1)`, its terms to apply
through their MPO. A run of factors of one site is a single term, applied as gates of one site.
Any other piece on at most two sites, or on sites of dimensions of product at most
`dense_limit`, laying no string elsewhere, is a dense block; otherwise it is developed, and
goes through its MPO, or through gates of one site when it is a single term.
"""
function piece_action(system::System, p::IndexedOp{IsOperator}, run::Bool)
    sites = sort(sites_of(p))
    if !run && (length(sites) ≤ 2 || prod(dim(system[k]) for k in sites) ≤ dense_limit)
        m = piece_matrix(system, p, sites)
        if stays_on(system, m, sites)
            return [ (m, sites) ], 1
        end
    end
    terms = lower(system, p)
    if isempty(terms)
        return Tuple{Matrix, Vector{Int}}[], 0
    elseif length(terms) == 1
        t = laid_strings(Pure, system, only(terms))
        return term_blocks(t), t.coef
    end
    return terms, 1
end

# A workaround of ITensors, which cannot factorise a tensor of charged indices that is exactly
# zero, having no block: when the result of applying gates or an MPO vanishes, the SVD, the QR
# or the eigen of ITensorMPS fails, and so does any later step on that null state. The functions
# below apply them in place of ITensorMPS, giving the null state in that case, so that apply is
# written as it would be without the bug.

"""
    vanished(st)

whether the MPS `st` is the null state for having a tensor of zero, as a gate of one site or a
null result of the functions below leaves it
"""
vanished(st::MPS) = any(t -> iszero(norm(t)), st)

"""
    guarded(f, vanishes, st)

the MPS `f(st)`, gates or an MPO applied to `st` by ITensorMPS, or the null state when the
result vanishes: when `f` fails on charged indices, `st` itself if it is already null, which
ITensors fails on, and `0 * st` if `vanishes(st)` tells exactly that the result is zero; any
other failure is raised again. `f` is given a copy, which a failure may leave half changed.
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
    vanishing_gate(g, st)

whether the gate `g` takes the MPS `st` exactly to zero: the norm of `g` on the block of the
sites it acts on, `st` orthogonalised on the first of them
"""
function vanishing_gate(g::ITensor, st::MPS)
    a, b = extrema(findsites(st, g))
    o = orthogonalize(st, a)
    return iszero(norm(g * prod(o[k] for k in a:b)))
end

"""
    gates_applied(gates, on, st, tr)

the MPS `st` with each gate of `gates` applied in turn on its consecutive sites `on`, given
rather than looked for as ITensorMPS does at a cost of order `n` per gate: a gate of one site
is contracted with the tensor of its site, any other applied by ITensorMPS, truncated by `tr`
"""
function gates_applied(gates::Vector{ITensor}, on::Vector{Vector{Int}}, st::MPS, tr)
    for (g, ks) in zip(gates, on)
        if length(ks) == 1
            st[only(ks)] = apply(g, st[only(ks)])
        else
            st = apply(g, st, ks; tr...)
        end
    end
    return st
end

"""
    applied_gates(gates, on, st, tr)

the MPS `st` with `gates` applied on their sites `on`, see `gates_applied`, truncated by `tr`,
or the null state when the result vanishes, see `guarded`: on a failure on charged indices, the
gates are applied again one by one, so that the one taking the state to zero is found
"""
function applied_gates(gates::Vector{ITensor}, on::Vector{Vector{Int}}, st::MPS, tr)
    try
        return gates_applied(gates, on, copy(st), tr)
    catch
        if !any(hasqns, st)
            rethrow()
        end
    end
    for g in gates
        st = guarded(s -> apply(g, s; tr...), s -> vanishing_gate(g, s), st)
    end
    return st
end

"""
    applied_mpo(mpo, st, limits; alg = "naive")

the MPS `st` with `mpo` applied by the algorithm `alg`, truncated to `limits`, or the null state
when the result vanishes, see `guarded`, which ``\\langle A\\psi | A\\psi \\rangle`` tells
exactly. The cutoff is `eps()` by default for `"densitymatrix"`, whose density matrices give the
squares of the singular values only to `eps()` of the largest: below, it would keep rounding.
"""
function applied_mpo(mpo::MPO, st::MPS, limits::Limits; alg::String = "naive")
    tr = alg == "densitymatrix" ? truncation(limits; cutoff = eps()) : truncation(limits)
    return guarded(s -> apply(mpo, s; alg, tr...), s -> iszero(inner(mpo, s, mpo, s)), st)
end

"""
    window_contracted(mpo, lo, st)

the MPS `st` with the tensors of its sites from `lo` on contracted with those of `mpo`, their
links joined, the others left as they are, `st` orthogonalised on `lo` first: the exact product
of `mpo` on these sites, whose orthogonality limits are left around them
"""
function window_contracted(mpo::MPO, lo::Int, st::MPS)
    ψ = orthogonalize(st, lo)
    ts = [ noprime(mpo[j] * ψ[lo + j - 1]) for j in eachindex(mpo) ]
    for j in 1:length(mpo) - 1
        k = lo + j - 1
        c = combiner(commonind(ψ[k], ψ[k + 1]), commonind(mpo[j], mpo[j + 1]); tags = "Link,l=$k")
        ts[j] *= c
        ts[j + 1] *= dag(c)
    end
    for (j, t) in enumerate(ts)
        ψ[lo + j - 1] = t
    end
    set_ortho_lims!(ψ, lo:lo + length(mpo) - 1)
    return ψ
end

"""
    applied_window(mpo, lo, st, tr)

the MPS `st` with `mpo` applied to its sites from `lo` on alone, see `window_contracted`, these
sites truncated by `tr` in one sweep, the exact product truncated as with `"naive"`, at a cost
that does not grow with the other sites; the null state when the result vanishes, see
`guarded`, which the norm of the exact product tells exactly
"""
applied_window(mpo::MPO, lo::Int, st::MPS, tr) =
    guarded(s -> truncate!(window_contracted(mpo, lo, s); site_range = lo:lo + length(mpo) - 1, tr...),
            s -> iszero(norm(window_contracted(mpo, lo, s))), st)

"""
    gate_mpo(g, idx)

the MPO of the gate `g` on the consecutive sites of indices `idx`, from the first site it acts
on to the last: `g` decomposed on its sites by singular value decompositions, which leave the
least links, and on the sites between the identity of the link alone, which leaves the index
of the site as it is
"""
function gate_mpo(g::ITensor, idx::Vector{<:Index})
    on = [ s for s in idx if hasind(g, s) ]
    m = MPO(g, on; cutoff = rounding_tol^2)
    ts = ITensor[]
    link = nothing
    j = 0
    for s in idx
        if hasind(g, s)
            j += 1
            t = isnothing(link) ? m[j] : replaceind(m[j], commonind(m[j - 1], m[j]) => link)
            push!(ts, t)
            link = j < length(m) ? commonind(m[j], m[j + 1]) : nothing
        else
            next = sim(link)
            push!(ts, denseblocks(delta(dag(link), next)))
            link = next
        end
    end
    return MPO(ts)
end

"""
    channel_factors(a)

the factors of the channel or the generator `a`, its coefficient left out, those of the
identity dropped, and each run of gates as the product of their operators, a gate of a product
being the product of the gates: an operator stands for its gate
"""
function channel_factors(a::IndexedOp)
    fs = IndexedOp[]
    for x in prodsubs(a)
        k = x isa Gate ? x.arg :
            x isa AtIndex && x.op isa Gate ? AtIndex(x.op.arg, x.index) : nothing
        if x isa IdentityOp
            continue
        elseif isnothing(k)
            push!(fs, x)
        elseif !isempty(fs) && last(fs) isa IndexedOp{IsOperator}
            fs[end] = last(fs) * k
        else
            push!(fs, k)
        end
    end
    return fs
end

"""
    mixed_block(system, terms, sites)

the matrix of the superoperator of terms `terms`, acting on `sites` alone, on the mixed indices
of those sites, see `mixed_elements`
"""
function mixed_block(system::System, terms::Vector{Term}, sites)
    orders = Dict(k => mixed_elements(system, k) for k in sites)
    return sum(terms) do t
        bs = Dict(t.blocks)
        t.coef * foldl(kron, [ haskey(bs, k) ? bs[k][orders[k], orders[k]] :
                                   Matrix{Float64}(I, dim(system[k])^2, dim(system[k])^2)
                               for k in sites ])
    end
end

"""
    GateBatch(gates, on)

a step of `apply`: gates applied together, each on its consecutive sites `on`
"""
struct GateBatch
    gates::Vector{ITensor}
    on::Vector{Vector{Int}}
end

"""
    WindowMPO

a step of `apply`: an MPO applied on the sites from `lo` on alone, see `applied_window`, that of
the terms of a piece, or of a gate on sites apart, see `gate_mpo`, which costs less than
bringing the sites together
"""
struct WindowMPO
    mpo::MPO
    lo::Int
end

"""
    add_gates!(steps, gates, on)

append `gates`, each on its consecutive sites `on`, to the steps of `apply`, joining them to the
last step when it is a batch of gates
"""
function add_gates!(steps::Vector, gates::Vector{ITensor}, on::Vector{Vector{Int}})
    if !isempty(steps) && last(steps) isa GateBatch
        append!(last(steps).gates, gates)
        append!(last(steps).on, on)
    else
        push!(steps, GateBatch(gates, on))
    end
    return steps
end

"""
    mpo_step!(steps, R, system, terms)

append to the steps of `apply` the MPO of `terms`, on the indices of representation `R`, on the
sites it acts on, see `mpo_window`, and return the number that brings: 0 when it acts nowhere,
the operator being zero, 1 otherwise
"""
function mpo_step!(steps::Vector, ::Type{R}, system::System, terms::Vector{Term}) where R
    pre = lay_terms!(PreMPO{R}(system), terms, 1)
    w = mpo_window(pre)
    if isnothing(w)
        return 0
    end
    push!(steps, WindowMPO(window_mpo(pre, [1.], w), first(w)))
    return 1
end

"""
    spread(sites)

whether the increasing `sites` leave others between them
"""
spread(ks) = last(ks) - first(ks) ≥ length(ks)

"""
    block_step!(steps, R, system, m, sites)

append to the steps of `apply` the dense block `m` on the increasing `sites`, on the indices of
representation `R`: a gate on consecutive sites, the MPO of the gate on sites apart, see
`WindowMPO`
"""
function block_step!(steps::Vector, ::Type{R}, system::System, m::AbstractMatrix, ks) where R
    if spread(ks)
        idx = [ siteind(R, system, k) for k in first(ks):last(ks) ]
        push!(steps, WindowMPO(gate_mpo(gate_on(m, idx[ks .- first(ks) .+ 1]), idx), first(ks)))
    else
        add_gates!(steps, [ gate_on(m, [ siteind(R, system, k) for k in ks ]) ], [ collect(ks) ])
    end
    return steps
end

"""
    lifted(system, w, k, left)

the tensor `w` of an MPO on the ket index of site `k` of `system`, made to act on its mixed
index: ``\\rho \\mapsto W\\rho`` when `left`, ``\\rho \\mapsto \\rho W^\\dagger`` otherwise,
whose links carry the opposite charges; a tensor without the site, the identity there, stays so
"""
function lifted(system::System, w::ITensor, k::Int, left::Bool)
    s, x = siteind(Pure, system, k), siteind(Mixed, system, k)
    if !hasind(w, s)
        # the identity on the site, see `gate_mpo`
        return left ? w : dag(w)
    end
    links = [ i for i in inds(w) if i ∉ (s, s') ]
    a = Array(w, s', dag(s), links...)
    d = dim(s)
    one = Matrix{eltype(a)}(I, d, d)
    o = mixed_elements(system, k)
    b = zeros(eltype(a), d^2, d^2, size(a)[3:end]...)
    for c in CartesianIndices(size(a)[3:end])
        m = a[:, :, c]
        b[:, :, c] = (left ? kron(one, m) : kron(conj(m), one))[o, o]
    end
    return charged_itensor(b, [ x', dag(x), (left ? links : dag.(links))... ])
end

"""
    gate_steps!(steps, system, m, sites)

append to the steps of `apply` on a density matrix the gate of matrix `m` on the increasing
`sites`: on one site its superoperator, ``\\rho \\mapsto m\\rho m^\\dagger``, as a gate; on
several, ``\\rho \\mapsto m\\rho`` then ``\\rho \\mapsto \\rho m^\\dagger``, the MPO of `m`,
see `gate_mpo`, made to act on the mixed indices, see `lifted`, each of the links of `m`: on
consecutive sites of dimensions small enough, see `dense_limit`, contracted into one gate
"""
function gate_steps!(steps::Vector, system::System, m::AbstractMatrix, ks)
    if length(ks) == 1
        o = mixed_elements(system, only(ks))
        return block_step!(steps, Mixed, system, kron(conj(m), m)[o, o], ks)
    end
    lo = first(ks)
    idx = [ siteind(Pure, system, k) for k in lo:last(ks) ]
    w = gate_mpo(gate_on(m, idx[ks .- lo .+ 1]), idx)
    for left in (true, false)
        ts = [ lifted(system, w[j], lo + j - 1, left) for j in eachindex(w) ]
        if !spread(ks) && prod(dim(system[k])^2 for k in ks) ≤ dense_limit
            add_gates!(steps, [ prod(ts) ], [ collect(ks) ])
        else
            push!(steps, WindowMPO(MPO(ts), lo))
        end
    end
    return steps
end

"""
    pure_steps(system, op)

the number a pure state is multiplied by and the steps of `apply` for the operator `op`, from
its rightmost piece to its leftmost, see `pieces`: each piece a dense block, see `block_step!`,
or its terms through their MPO on the sites they act on, see `mpo_step!`
"""
function pure_steps(system::System, a::IndexedOp{IsOperator})
    coef, ps = pieces(system, a)
    steps = Union{GateBatch, WindowMPO}[]
    for (p, run) in reverse(ps)
        action, c = piece_action(system, p, run)
        coef *= c
        if coef == 0
            break
        elseif action isa Vector{Term}
            coef *= mpo_step!(steps, Pure, system, action)
        else
            foreach(((m, ks),) -> block_step!(steps, Pure, system, m, ks), action)
        end
    end
    return coef, steps
end

"""
    mixed_steps(system, op)

the number a density matrix is multiplied by and the steps of `apply` for the channel or the
generator `op`, from its rightmost factor to its leftmost, see `channel_factors`: a gate piece
by piece, see `gate_steps!`, any other factor as a dense block on the mixed indices, see
`block_step!`, or the terms of either through their MPO, see `mpo_step!`
"""
function mixed_steps(system::System, a::IndexedOp)
    coef = scalarcoef(a)
    steps = Union{GateBatch, WindowMPO}[]
    for f in reverse(channel_factors(scalararg(a)))
        if f isa IndexedOp{IsOperator}
            c, ps = pieces(system, f)
            coef *= abs2(c)
            for (p, run) in reverse(ps)
                action, c = piece_action(system, p, run)
                coef *= abs2(c)
                if action isa Vector{Term}
                    coef *= mpo_step!(steps, Mixed, system, gate_terms(action))
                else
                    foreach(((m, ks),) -> gate_steps!(steps, system, m, ks), action)
                end
            end
        else
            terms = lower(system, f)
            laid = [ laid_strings(Mixed, system, t) for t in terms ]
            ks = sort!(unique(k for t in laid for (k, _) in t.blocks))
            if isempty(ks)
                coef *= sum(t.coef for t in terms; init = 0.)
            elseif prod(dim(system[k])^2 for k in ks) ≤ dense_limit
                block_step!(steps, Mixed, system, mixed_block(system, laid, ks), ks)
            else
                coef *= mpo_step!(steps, Mixed, system, terms)
            end
        end
        if coef == 0
            break
        end
    end
    return coef, steps
end

"""
    applied_steps(st, steps, limits)

the MPS `st` with the steps of `apply` applied in order, each truncated to `limits`: batches of
gates, see `applied_gates`, and MPOs on the sites they act on, see `applied_window`
"""
function applied_steps(st::MPS, steps, limits::Limits)
    tr = truncation(limits)
    for step in steps
        st = step isa WindowMPO ? applied_window(step.mpo, step.lo, st, tr) :
                                  applied_gates(step.gates, step.on, st, tr)
    end
    return st
end

"""
    apply(op, state; limits = Limits())
    apply(mpo::MPO, state; limits = Limits(), apply_algo = "naive")

the state with the operator `op`, placed on sites, applied, or the MPO `mpo` from
`make_mpo`. A product acts as the operator it denotes, its rightmost factor first, factor by
factor and never developed; each factor is applied as a dense block on its sites, or through
its MPO when it is too large or lays Jordan-Wigner strings on other sites. On a pure state, a
channel of a single Kraus operator `K`, as `Gate(K)`, is applied as `K`. On a mixed state, a
channel is applied factor by factor in the same way, `Gate(K)` as `K` on the kets and its
conjugate on the bras, an operator `K` being `Gate(K)`; a generator `L` gives `L(ρ)`. The
state is not normalised.

- `limits`: the truncation of the bonds a gate of several sites crosses, or of the whole state
  for an MPO, `cutoff = 1e-26` by default, `eps()` with `"densitymatrix"`, whose density
  matrices give the squares of the singular values only to `eps()` of the largest
- `apply_algo`: for an MPO, the algorithm of its product with the state: `"naive"`, the exact
  product truncated, or `"densitymatrix"`, faster when it truncates much, but whose error may
  stay near `1e-8`, see Camaño, Epperly and Tropp, arXiv:2504.06475

# Examples

    apply(controlled(Z)(1, 3) * H(2) * controlled(X)(3, 4), state)
    apply(exp(-0.1im * (X ⊗ X + Y ⊗ Y))(2, 5), state)
    apply(depolarizing_gate(0.01, 2)(1, 2) * controlled(X)(1, 2), ρ)
"""
function apply(a::IndexedOp{IsOperator}, state::State{Pure}; limits::Limits = Limits())
    sys = state.system
    check_indices(sys, a)
    coef, steps = pure_steps(sys, a)
    if coef == 0
        return State(state, 0 * state.state)
    end
    st = applied_steps(state.state, steps, limits)
    return State(state, coef == 1 ? st : coef * st)
end

function apply(a::IndexedOp{IsChannel}, state::State{Pure}; limits::Limits = Limits())
    if !is_pure(a)
        error("$a has more than one Kraus operator, which a pure state cannot take: apply it to " *
              "mix(state)")
    end
    return apply(kraus(a), state; limits)
end

apply(a::IndexedOp{IsGenerator}, ::State{Pure}; limits::Limits = Limits()) =
    error("apply takes an operator or a channel, and $a is a generator")

apply(a::IndexedOp{IsOperator}, state::State{Mixed}; limits::Limits = Limits()) =
    apply(Gate(a), state; limits)

function apply(a::IndexedOp{<:Union{IsChannel, IsGenerator}}, state::State{Mixed};
               limits::Limits = Limits())
    sys = state.system
    check_indices(sys, a)
    coef, steps = mixed_steps(sys, a)
    if coef == 0
        return State(state, 0 * state.state)
    end
    st = applied_steps(state.state, steps, limits)
    return State(state, coef == 1 ? st : coef * st)
end

function apply(mpo::MPO, state::State; limits::Limits = Limits(), apply_algo::String = "naive")
    check_apply_algo(apply_algo)
    return State(state, applied_mpo(mpo, state.state, limits; alg = apply_algo))
end
