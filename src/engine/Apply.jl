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

the largest product of the squares of the dimensions of the sites for which `apply` applies to a
density matrix as one block on their mixed indices a gate of two neighbouring sites, see
`gate_steps!`, or another channel or a generator of one site or two neighbouring ones, see
`mixed_steps`; beyond, a gate goes through the MPOs of its sides, any other through its MPO
"""
const dense_limit = 1024

"""
    dense_enough(R, system, sites)

whether a block on `sites` of `system`, on the indices of representation `R`, is small enough
to be dense, see `dense_limit`
"""
function dense_enough(::Type{R}, system::System, ks) where R
    # a Float64, which does not wrap around as an Int does on many sites
    return prod(float(block_dim(R, system[k])) for k in ks) ≤ dense_limit
end

"""
    check_apply_algo(apply_algo)

refuse an algorithm of the product of an MPO by a state other than `"naive"`,
`"densitymatrix"` and `"zipup"`: `"fit"` needs sweeps of its own
"""
function check_apply_algo(apply_algo::String)
    if apply_algo ∉ ("naive", "densitymatrix", "zipup")
        error("apply_algo is \"naive\", \"densitymatrix\" or \"zipup\", not $(repr(apply_algo))")
    end
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
    neighbouring(sites)

whether the increasing `sites` are one site or two neighbouring ones, where `apply` takes a
piece as a dense block
"""
neighbouring(sites) = length(sites) == 1 || (length(sites) == 2 && sites[2] == sites[1] + 1)

"""
    grouped(system, fs)

the factors `fs` of a product, in their order, with each run of consecutive factors on the same
two neighbouring sites holding a factor of these two sites taken as one, see `neighbouring`,
the factors of one site among them included, those at its ends only when even, so that an odd
one keeps its strings in a run, see `pieces`
"""
function grouped(system::System, fs::Vector{IndexedOp{IsOperator}})
    out = IndexedOp{IsOperator}[]
    seg = IndexedOp{IsOperator}[]
    sites = Int[]
    function close!()
        pairs = findall(x -> length(sites_of(x)) == 2, seg)
        if isempty(pairs)
            append!(out, seg)
        else
            lo, hi = first(pairs), last(pairs)
            even(x) = parity_on(system, x) == 0
            while lo > 1 && even(seg[lo - 1])
                lo -= 1
            end
            while hi < length(seg) && even(seg[hi + 1])
                hi += 1
            end
            append!(out, seg[1:lo-1])
            push!(out, ProdOp(seg[lo:hi]))
            append!(out, seg[hi+1:end])
        end
        empty!(seg)
        empty!(sites)
    end
    for x in fs
        u = sort(union(sites, sites_of(x)))
        if !neighbouring(u)
            close!()
            u = sort(sites_of(x))
        end
        if neighbouring(u)
            push!(seg, x)
            append!(empty!(sites), u)
        else
            push!(out, x)
        end
    end
    close!()
    return out
end

"""
    pieces(system, op)

the factors of the product `op`, gathered, see `gathered`, cut into the pieces `apply` applies
one by one, and the number they are multiplied by: the factors on two neighbouring sites taken
together, see `grouped`, then the factors of one site of definite parity brought together where
they commute with the pieces they cross, a sign for two odd ones, so that their strings
compensate, every other factor being a piece of its own, each piece with whether it is such a
run. On a charged system, a product with a group of factors whose terms change the charges by
different amounts is one piece, see `flux_groups`.
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
    fs = IndexedOp{IsOperator}[ x for x in prodsubs(scalararg(g)) if !(x isa IdentityOp) ]
    groups = grouped(system, fs)
    # a piece needs a definite flux on charged sites, which a group with none does not give
    function fluxes(f)
        gs = flux_groups(system, f)
        return isnothing(gs) ? 1 : length(gs)
    end
    if is_charged(system) && any(x -> fluxes(x) > 1, groups)
        return coef, [ (ProdOp(fs), false) ]
    end
    for x in groups
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

# The workaround of ITensors, which cannot factorise a tensor of charged indices that is exactly
# zero, see `guarded`: the functions below apply gates or an MPO in place of ITensorMPS, giving
# the null state when the result vanishes, so that apply is written as it would be without the
# bug.

"""
    vanishing_gate(g, sites, st)

whether the gate `g` on the consecutive `sites` takes the MPS `st` exactly to zero: the norm of
`g` on the block of these sites, `st` orthogonalised on the first of them
"""
function vanishing_gate(g::ITensor, ks::Vector{Int}, st::MPS)
    a, b = extrema(ks)
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
    for (g, ks) in zip(gates, on)
        st = guarded(s -> gates_applied([g], [ks], s, tr), s -> vanishing_gate(g, ks, s), st)
    end
    return st
end

"""
    twice_maxdim(maxdim)

twice the bond dimension `maxdim`, `typemax(Int)` from `typemax(Int) ÷ 2` on, where the double
would wrap around
"""
twice_maxdim(maxdim::Int) = maxdim ≥ typemax(Int) ÷ 2 ? typemax(Int) : 2maxdim

"""
    zipup_truncation(tr)

the keywords of zipup for the truncation `tr` of the product: at each site, as it is formed,
a tenth of the cutoff and twice the bond dimension, as ITensorMPS recommends after Paeckel et
al., these truncations not being optimal, then `tr` once the product is whole
"""
zipup_truncation(tr::NamedTuple) =
    (cutoff = tr.cutoff / 10, maxdim = twice_maxdim(tr.maxdim), mindim = tr.mindim,
     truncate_kwargs = tr)

"""
    applied_mpo(mpo, st, limits; alg = "naive")

the MPS `st` with `mpo` applied by the algorithm `alg`, truncated to `limits`, or the null state
when the result vanishes, see `guarded`, which ``\\langle A\\psi | A\\psi \\rangle`` tells
exactly. The cutoff is `eps()` by default for `"densitymatrix"`, whose density matrices give the
squares of the singular values only to `eps()` of the largest: below, it would keep rounding.
`"zipup"` truncates more loosely as it forms the product, see `zipup_truncation`.
"""
function applied_mpo(mpo::MPO, st::MPS, limits::Limits; alg::String = "naive")
    tr = alg == "densitymatrix" ? truncation(limits; cutoff = eps()) : truncation(limits)
    kw = alg == "zipup" ? zipup_truncation(tr) : tr
    return guarded(s -> apply(mpo, s; alg, kw...), s -> iszero(inner(mpo, s, mpo, s)), st)
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
function applied_window(mpo::MPO, lo::Int, st::MPS, tr)
    sites = lo:lo + length(mpo) - 1
    return guarded(s -> truncate!(window_contracted(mpo, lo, s); site_range = sites, tr...),
                   s -> iszero(norm(window_contracted(mpo, lo, s))), st)
end

"""
    gate_mpo(R, system, m, sites)

the MPO of the gate of matrix `m` on the increasing `sites` of `system`, on the indices of
representation `R`, from the first of them to the last: the gate decomposed on its sites by
singular value decompositions, which leave the least links, and on the sites between the
identity of the link alone, which leaves the index of the site as it is
"""
function gate_mpo(::Type{R}, system::System, m::AbstractMatrix, ks) where R
    idx = [ siteind(R, system, k) for k in ks ]
    # a dense gate of many sites has as many indices as ITensors warns about, from 14
    w = ITensors.@disable_warn_order MPO(gate_on(m, idx), idx; cutoff = rounding_tol^2)
    ts = ITensor[]
    link = nothing
    j = 0
    for k in first(ks):last(ks)
        if k in ks
            j += 1
            t = isnothing(link) ? w[j] : replaceind(w[j], commonind(w[j - 1], w[j]) => link)
            push!(ts, t)
            link = j < length(w) ? commonind(w[j], w[j + 1]) : nothing
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
    add_gate!(steps, g, sites)

append the gate `g` on its consecutive `sites` to the steps of `apply`, joining it to the last
step when that is a batch of gates
"""
function add_gate!(steps::Vector, g::ITensor, ks)
    if !isempty(steps) && last(steps) isa GateBatch
        push!(last(steps).gates, g)
        push!(last(steps).on, collect(ks))
    else
        push!(steps, GateBatch([g], [collect(ks)]))
    end
    return steps
end

"""
    pre_step!(steps, pre)

append to the steps of `apply` the MPO of `pre` on the sites it acts on, see `mpo_window`, and
return the number that brings: 0 when it acts nowhere, the operator being zero, 1 otherwise
"""
function pre_step!(steps::Vector, pre::PreMPO)
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
        push!(steps, WindowMPO(gate_mpo(R, system, m, ks), first(ks)))
    else
        add_gate!(steps, gate_on(m, [ siteind(R, system, k) for k in ks ]), ks)
    end
    return steps
end

"""
    on_mixed_index(system, w, k, left, bra)

the tensor `w` of an MPO on the ket index of site `k` of `system`, made to act on its mixed
index: ``\\rho \\mapsto W\\rho`` when `left`, ``\\rho \\mapsto \\rho W^\\dagger`` otherwise,
whose links carry the opposite charges, relabelled by `bra`, a `relabeller` the tensors of one
MPO share, to the starred names the bra takes under a strong symmetry; a tensor without the
site, the identity there, stays so
"""
function on_mixed_index(system::System, w::ITensor, k::Int, left::Bool, bra)
    s, x = siteind(Pure, system, k), siteind(Mixed, system, k)
    if !hasind(w, s)
        # the identity on the site, see `gate_mpo`
        return left ? w : relabel(dag(w), bra)
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
    return charged_itensor(b, [ x', dag(x), (left ? links : dag.(bra.(links)))... ])
end

"""
    gate_steps!(steps, R, system, m, sites)

append to the steps of `apply` the gate of matrix `m` on the increasing `sites`: on a pure
state, `R = Pure`, as a dense block, see `block_step!`. On a density matrix, on one site, its
superoperator ``\\rho \\mapsto m\\rho m^\\dagger`` as a block; on several,
``\\rho \\mapsto m\\rho`` then ``\\rho \\mapsto \\rho m^\\dagger``, each the MPO of `m`, see
`gate_mpo`, made to act on the mixed indices, see `on_mixed_index`: one gate on consecutive
sites small enough, see `dense_enough`, an MPO on its sites otherwise
"""
gate_steps!(steps::Vector, ::Type{Pure}, system::System, m::AbstractMatrix, ks) =
    block_step!(steps, Pure, system, m, ks)

function gate_steps!(steps::Vector, ::Type{Mixed}, system::System, m::AbstractMatrix, ks)
    if length(ks) == 1
        o = mixed_elements(system, only(ks))
        return block_step!(steps, Mixed, system, kron(conj(m), m)[o, o], ks)
    end
    lo = first(ks)
    w = gate_mpo(Pure, system, m, ks)
    one_gate = !spread(ks) && dense_enough(Mixed, system, ks)
    bra = relabeller(i -> star(i, strong_names(system)))
    for left in (true, false)
        ts = [ on_mixed_index(system, w[j], lo + j - 1, left, bra) for j in eachindex(w) ]
        if one_gate
            add_gate!(steps, prod(ts), ks)
        else
            push!(steps, WindowMPO(MPO(ts), lo))
        end
    end
    return steps
end

"""
    one_site_gates!(steps, R, system, t)

append to the steps of `apply` the blocks of the term `t` of an operator, its strings laid, as
gates of one site, see `gate_steps!`, and return its coefficient
"""
function one_site_gates!(steps::Vector, ::Type{R}, system::System, t::Term) where R
    t = laid_strings(Pure, system, t)
    foreach(b -> gate_steps!(steps, R, system, b.m, [b.site]), t.blocks)
    return t.coef
end

"""
    piece_steps!(steps, R, system, p, run)

append to the steps of `apply` on a state of representation `R` those of the piece `p` of an
operator, `run` when it is a run of factors of one site, see `pieces`, and return the number
they bring, 0 for a dense piece that is zero. A piece on one site or two neighbouring ones, see `neighbouring`, laying no
string elsewhere and changing the charges by a definite amount, see `definite_on`, is a dense
block, see `gate_steps!`, unless it is a run; an operator given by its matrix on more sites goes
through its MPO as any other piece. Otherwise it goes through gates of one site when it has a
single term, as a run usually does, and through its MPO otherwise, see `compacted_op`, on a
density matrix that of ``\\rho \\mapsto K\\rho`` then that of ``\\rho \\mapsto \\rho K^\\dagger``
"""
function piece_steps!(steps::Vector, ::Type{R}, system::System, p::IndexedOp{IsOperator},
                      run::Bool) where R
    sites = sort(sites_of(p))
    if !run && neighbouring(sites)
        m = matrix_on(system, p, sites)
        # the null state, which ITensors cannot factorise into an MPO on charged sites
        if iszero(m)
            return 0
        elseif stays_on(system, m, sites) && definite_on(system, m, sites)
            gate_steps!(steps, R, system, m, sites)
            return 1
        end
    end
    c = compacted_op(system, p)
    if !isnothing(c.terms) && length(c.terms) == 1
        return one_site_gates!(steps, R, system, only(c.terms))
    elseif R === Pure
        return pre_step!(steps, PreMPO(system, [c]))
    end
    # the sides one after the other, not their product, whose reduction can lose the parts of
    # small relative size of the operator, nor `Gate(p)`, made of the gates of the factors, of
    # no definite flux when they have none
    sides = [ compacted_op(system, p; side = s) for s in (:left, :right) ]
    return prod(pre_step!(steps, PreMPO(system, [c])) for c in sides)
end

"""
    operator_steps!(steps, R, system, op)

append to the steps of `apply` on a state of representation `R` those of the operator `op`,
piece by piece from its rightmost to its leftmost, see `pieces` and `piece_steps!`, and return
the number they bring
"""
function operator_steps!(steps::Vector, ::Type{R}, system::System, a::IndexedOp{IsOperator}) where R
    coef, ps = pieces(system, a)
    for (p, run) in reverse(ps)
        coef *= piece_steps!(steps, R, system, p, run)
        if coef == 0
            break
        end
    end
    return coef
end

"""
    pure_steps(system, op)

the number a pure state is multiplied by and the steps of `apply` for the operator `op`, see
`operator_steps!`
"""
function pure_steps(system::System, a::IndexedOp{IsOperator})
    steps = Union{GateBatch, WindowMPO}[]
    return operator_steps!(steps, Pure, system, a), steps
end

"""
    mixed_steps(system, op)

the number a density matrix is multiplied by and the steps of `apply` for the channel or the
generator `op`, from its rightmost factor to its leftmost, see `channel_factors`: a gate as
the pieces of its operator, see `piece_steps!`; any other factor on one site or two neighbouring
ones, see `neighbouring`, small enough, see `dense_enough`, as a dense block on the mixed
indices, see `block_step!`, and through its MPO otherwise, see `compacted_super`
"""
function mixed_steps(system::System, a::IndexedOp)
    coef = scalarcoef(a)
    steps = Union{GateBatch, WindowMPO}[]
    for f in reverse(channel_factors(a))
        if f isa IndexedOp{IsOperator}
            coef *= abs2(operator_steps!(steps, Mixed, system, f))
        else
            terms = [ laid_strings(Mixed, system, t) for t in lower(system, f) ]
            ks = sort!(unique(b.site for t in terms for b in t.blocks))
            if isempty(ks)
                coef *= sum(t.coef for t in terms; init = 0.)
            elseif neighbouring(ks) && dense_enough(Mixed, system, ks)
                orders = [ mixed_elements(system, k) for k in ks ]
                m = system_matrix(Mixed, system, terms, ks, orders)
                block_step!(steps, Mixed, system, m, ks)
            else
                coef *= pre_step!(steps, PreMPO(system, [compacted_super(system, f)]))
            end
        end
        if coef == 0
            break
        end
    end
    return coef, steps
end

"""
    step_applied(st, coef, steps, k, tr)

the MPS `st` after step `k` of `apply`, `k` running over `1:max(1, length(steps))`: step `k`
applied, truncated by `tr`, see `applied_step`, and the MPS multiplied by `coef` at the last
step, or made null at the first when `coef` is zero; a step past them leaves it
"""
function step_applied(st::MPS, coef::Number, steps, k::Int, tr)
    if coef == 0
        return k == 1 ? 0 * st : st
    elseif k ≤ length(steps)
        st = applied_step(st, steps[k], tr)
    end
    return k == max(1, length(steps)) && coef ≠ 1 ? coef * st : st
end

"""
    applied(state, coef, steps, limits)

the state with the steps of `apply` applied in order, each truncated to `limits`, batches of
gates, see `applied_gates`, and MPOs on the sites they act on, see `applied_window`, and
multiplied by `coef`; the null state when `coef` is zero
"""
function applied(state::State, coef::Number, steps, limits::Limits)
    st, tr = state.state, truncation(limits)
    for k in 1:max(1, length(steps))
        st = step_applied(st, coef, steps, k, tr)
    end
    return State(state, st)
end

"""
    applied_step(st, step, tr)

the MPS `st` with the step `step` of `apply` applied, truncated by `tr`, see `truncation`: a batch
of gates, see `applied_gates`, or an MPO on the sites it acts on, see `applied_window`
"""
applied_step(st::MPS, step, tr) =
    step isa WindowMPO ? applied_window(step.mpo, step.lo, st, tr) :
                         applied_gates(step.gates, step.on, st, tr)

"""
    apply(op, state; limits = Limits())
    apply(mpo::MPO, state; limits = Limits(), apply_algo = "naive")

the state with the operator `op`, placed on sites, applied, or the MPO `mpo` from
`make_mpo`. A product acts as the operator it denotes, its rightmost factor first, factor by
factor and never developed: the factors on one site or two neighbouring sites are applied as a
dense block there, together when they share these sites, any other through its MPO, as a
factor laying Jordan-Wigner strings on other sites; on conserved charges, a product with a
factor of no definite flux outside such a block goes whole through its MPO. On a pure state, a
channel of a single Kraus operator `K`, as `Gate(K)`, is applied as `K`. On a mixed state, a
channel is applied factor by factor in the same way, `Gate(K)` as `K` on the kets and its
conjugate on the bras, an operator `K` being `Gate(K)`; a generator `L` gives `L(ρ)`, a
dephasing of more than two sites through the MPO of the gates of its projectors when they have
few terms. The state is not normalised.

- `limits`: the truncation of the bonds a gate of several sites crosses, or of the whole state
  for an MPO, `cutoff = 1e-26` by default, `eps()` with `"densitymatrix"`, whose density
  matrices give the squares of the singular values only to `eps()` of the largest
- `apply_algo`: for an MPO, the algorithm of its product with the state: `"naive"`, the exact
  product truncated, `"densitymatrix"`, faster when it truncates much, but whose error may stay
  near `1e-8`, or `"zipup"`, which truncates at each site as it multiplies, to a tenth of the
  cutoff and twice the bond dimension, these truncations not being optimal, then the result to
  `limits`: faster still when it truncates much, but less accurate; see Camaño, Epperly and
  Tropp, arXiv:2504.06475

# Examples

    apply(controlled(Z)(1, 3) * H(2) * controlled(X)(3, 4), state)
    apply(exp(-0.1im * (X ⊗ X + Y ⊗ Y))(2, 5), state)
    apply(depolarizing_channel(0.01, 2)(1, 2) * controlled(X)(1, 2), ρ)
"""
apply(a::IndexedOp, state::State; limits::Limits = Limits()) =
    applied(state, apply_steps(a, state)..., limits)

"""
    apply_steps(op, state)

the number `state` is multiplied by and the steps of `apply` for `op`, see `pure_steps` and
`mixed_steps`, which depend on the system of the state alone
"""
function apply_steps(a::IndexedOp{IsOperator}, state::State{Pure})
    check_indices(state.system, a)
    return pure_steps(state.system, a)
end

function apply_steps(a::IndexedOp{IsChannel}, state::State{Pure})
    if !is_pure(a)
        error("$(obs_name(a)) has more than one Kraus operator, which a pure state cannot " *
              "take: apply it to mix(state)")
    end
    return apply_steps(kraus(a), state)
end

apply_steps(a::IndexedOp{IsGenerator}, ::State{Pure}) =
    error("apply takes an operator or a channel, and $(obs_name(a)) is a generator")

apply_steps(a::IndexedOp{IsOperator}, state::State{Mixed}) = apply_steps(Gate(a), state)

function apply_steps(a::IndexedOp{<:Union{IsChannel, IsGenerator}}, state::State{Mixed})
    check_indices(state.system, a)
    return mixed_steps(state.system, a)
end

function apply(mpo::MPO, state::State; limits::Limits = Limits(), apply_algo::String = "naive")
    check_apply_algo(apply_algo)
    return State(state, applied_mpo(mpo, state.state, limits; alg = apply_algo))
end
