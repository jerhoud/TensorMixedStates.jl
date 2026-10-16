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
`apply` takes as a dense block, and of their squares for a block on the mixed indices of a
density matrix; beyond, it goes through its MPO
"""
const dense_limit = 1024

"""
    dense_enough(R, system, sites)

whether a block on `sites` of `system`, on the indices of representation `R`, is small enough
to be dense, see `dense_limit`
"""
dense_enough(::Type{R}, system::System, ks) where R = prod(block_dim(R, system[k]) for k in ks) ≤ dense_limit

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
    for x in prodsubs(g)
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
    gate_mpo(R, system, m, sites)

the MPO of the gate of matrix `m` on the increasing `sites` of `system`, on the indices of
representation `R`, from the first of them to the last: the gate decomposed on its sites by
singular value decompositions, which leave the least links, and on the sites between the
identity of the link alone, which leaves the index of the site as it is
"""
function gate_mpo(::Type{R}, system::System, m::AbstractMatrix, ks) where R
    idx = [ siteind(R, system, k) for k in ks ]
    w = MPO(gate_on(m, idx), idx; cutoff = rounding_tol^2)
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
    mpo_step!(steps, R, system, terms)

append to the steps of `apply` the MPO of `terms`, on the indices of representation `R`, on the
sites it acts on, see `mpo_window`, and return the number that brings: 0 when it acts nowhere,
the operator being zero, 1 otherwise
"""
mpo_step!(steps::Vector, ::Type{R}, system::System, terms::Vector{Term}) where R =
    pre_step!(steps, laid(system, [from_terms(R, system, "a piece", terms)]))

"""
    pre_step!(steps, pre)

append to the steps of `apply` the MPO of `pre` on the sites it acts on, see `mpo_step!`, and
return the number that brings
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
    lifted(system, w, k, left, bra)

the tensor `w` of an MPO on the ket index of site `k` of `system`, made to act on its mixed
index: ``\\rho \\mapsto W\\rho`` when `left`, ``\\rho \\mapsto \\rho W^\\dagger`` otherwise,
whose links carry the opposite charges, relabelled by `bra`, a `relabeller` the tensors of one
MPO share, to the starred names the bra takes under a strong symmetry; a tensor without the
site, the identity there, stays so
"""
function lifted(system::System, w::ITensor, k::Int, left::Bool, bra)
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
`gate_mpo`, made to act on the mixed indices, see `lifted`: one gate on consecutive sites small
enough, see `dense_enough`, an MPO on its sites otherwise
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
        ts = [ lifted(system, w[j], lo + j - 1, left, bra) for j in eachindex(w) ]
        if one_gate
            add_gate!(steps, prod(ts), ks)
        else
            push!(steps, WindowMPO(MPO(ts), lo))
        end
    end
    return steps
end

"""
    piece_steps!(steps, R, system, p, run)

append to the steps of `apply` on a state of representation `R` those of the piece `p` of an
operator, `run` when it is a run of factors of one site, see `pieces`, and return the number
they bring, 0 for a dense piece that is zero. A piece on at most two sites, or small enough, see `dense_enough`, laying no string
elsewhere and changing the charges by a definite amount, see `definite_on`, is a dense block,
see `gate_steps!`, unless it is a run; otherwise it is developed, the descent refusing it when
the charges change by no definite amount, and goes through gates of one site when it gives a
single term, as a run usually does, through
the MPO of its terms otherwise, see `mpo_step!`, on a density matrix those of
``\\rho \\mapsto K\\rho K^\\dagger``, see `gate_terms`
"""
function piece_steps!(steps::Vector, ::Type{R}, system::System, p::IndexedOp{IsOperator},
                      run::Bool) where R
    sites = sort(sites_of(p))
    if !run && (length(sites) ≤ 2 || dense_enough(Pure, system, sites))
        m = matrix_on(system, p, sites)
        # the null state, which ITensors cannot factorise into an MPO on charged sites
        if iszero(m)
            return 0
        elseif stays_on(system, m, sites) && definite_on(system, m, sites)
            gate_steps!(steps, R, system, m, sites)
            return 1
        end
    end
    check_indices(system, p)
    g = gathered(system, p)
    # a product of sums that would develop into more terms than its factors have, or, on a
    # density matrix, Kρ K† of many terms
    if R === Pure && any(by_blocks, gathered_summands(g))
        return pre_step!(steps, laid(system, [compacted_op(system, p)]))
    elseif R === Mixed && many_terms(g)
        return pre_step!(steps, laid(system, [product(sided(system, p, true), sided(system, p, false))]))
    end
    terms = lower(system, p, g)
    if length(terms) == 1
        t = laid_strings(Pure, system, only(terms))
        foreach(((k, m),) -> gate_steps!(steps, R, system, m, [k]), t.blocks)
        return t.coef
    end
    return mpo_step!(steps, R, system, R === Mixed ? gate_terms(terms) : terms)
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
generator `op`, from its rightmost factor to its leftmost, see `channel_factors`: a gate piece
by piece, see `gate_steps!`, any other factor as a dense block on the mixed indices, see
`block_step!`, or the terms of either through their MPO, see `mpo_step!`
"""
function mixed_steps(system::System, a::IndexedOp)
    coef = scalarcoef(a)
    steps = Union{GateBatch, WindowMPO}[]
    for f in reverse(channel_factors(a))
        if f isa IndexedOp{IsOperator}
            coef *= abs2(operator_steps!(steps, Mixed, system, f))
        else
            terms = lower(system, f)
            laid = [ laid_strings(Mixed, system, t) for t in terms ]
            ks = sort!(unique(k for t in laid for (k, _) in t.blocks))
            if isempty(ks)
                coef *= sum(t.coef for t in terms; init = 0.)
            elseif dense_enough(Mixed, system, ks)
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
    applied(state, coef, steps, limits)

the state with the steps of `apply` applied in order, each truncated to `limits`, batches of
gates, see `applied_gates`, and MPOs on the sites they act on, see `applied_window`, and
multiplied by `coef`; the null state when `coef` is zero
"""
function applied(state::State, coef::Number, steps, limits::Limits)
    if coef == 0
        return State(state, 0 * state.state)
    end
    st, tr = state.state, truncation(limits)
    for step in steps
        st = step isa WindowMPO ? applied_window(step.mpo, step.lo, st, tr) :
                                  applied_gates(step.gates, step.on, st, tr)
    end
    return State(state, coef == 1 ? st : coef * st)
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
    check_indices(state.system, a)
    return applied(state, pure_steps(state.system, a)..., limits)
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
    check_indices(state.system, a)
    return applied(state, mixed_steps(state.system, a)..., limits)
end

function apply(mpo::MPO, state::State; limits::Limits = Limits(), apply_algo::String = "naive")
    check_apply_algo(apply_algo)
    return State(state, applied_mpo(mpo, state.state, limits; alg = apply_algo))
end
