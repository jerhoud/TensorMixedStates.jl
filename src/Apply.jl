# apply on pure states: a product of operators applied factor by factor, right to left, never
# developed, each piece as a dense block on its sites, as one-site gates, or through its MPO. A
# channel of a single Kraus operator is applied as that operator. On mixed states, the factors
# of a channel likewise: a gate as its operator on the kets and its conjugate on the bras, any
# other channel or a generator as a dense block on the mixed indices or through its MPO.

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
    fermionic(k) = matrix(F, system[k]) ≠ I
    if any(fermionic, 1:first(sites)-1) && parity_of(m, foldl(kron, fs)) ≠ 0
        return false
    end
    for t in first(sites)+1:last(sites)-1
        if !(t in sites) && fermionic(t)
            p = foldl(kron, [ k > t ? f : one(f) for (k, f) in zip(sites, fs) ])
            if !nearly(p * m * p, m)
                return false
            end
        end
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
    info = Tuple{Vector{Int}, Union{Nothing, Int}}[]
    runs = Bool[]
    for x in prodsubs(scalararg(g))
        if x isa IdentityOp
            continue
        end
        sx, px = sites_of(x), parity_on(system, x)
        joins = x isa AtIndex && length(x.index) == 1 && !isnothing(px)
        # the last run of factors of one site this one can reach
        j = joins ? findlast(runs) : nothing
        if !isnothing(j) && all(k -> !dependent(info[k], (sx, px)), j+1:length(ps))
            for k in j+1:length(ps)
                if info[k][2] == 1 && px == 1
                    coef = -coef
                end
            end
            push!(ps[j], x)
            info[j] = (union(info[j][1], sx), mod(info[j][2] + px, 2))
        else
            push!(ps, [x])
            push!(info, (sx, px))
            push!(runs, joins)
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
result vanishes: `st` itself when it is already null, which ITensors would fail on, and when
`f` fails on charged indices, `0 * st` if `vanishes(st)` tells exactly that the result is
zero; any other failure is raised again. `f` is given a copy, which a failure may leave half
changed.
"""
function guarded(f, vanishes, st::MPS)
    if vanished(st)
        return st
    end
    try
        return f(copy(st))
    catch
        if any(hasqns, st) && vanishes(st)
            return 0 * st
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
    sites_back(ψ, st, tr)

the MPS `st`, made from `ψ` by gates that moved its sites, with each site back where it is in
`ψ`, the swaps truncated by `tr`: only the window of the positions that moved is reordered, the
sites outside it staying in place, where ITensorMPS reorders all of them at a cost of order `n²`
"""
function sites_back(ψ::MPS, st::MPS, tr)
    place = Dict(i => k for k in 1:length(ψ) for i in siteinds(ψ, k))
    now = [ place[first(siteinds(st, n))] for n in 1:length(st) ]
    moved = findall(n -> now[n] ≠ n, eachindex(now))
    if isempty(moved)
        return st
    end
    lo, hi = extrema(moved)
    return movesites(st, [ n => now[n] for n in lo:hi ]; tr...)
end

"""
    gates_applied(gates, on, st, tr)

the MPS `st` with each gate of `gates` applied in turn on its sites `on`, given rather than looked
for as ITensorMPS does at a cost of order `n` per gate: a gate of one site is contracted with
the tensor of its site, any other applied by ITensorMPS on the positions its sites have then,
truncated by `tr`, the sites it brings next to each other left there, see `sites_back`
"""
function gates_applied(gates::Vector{ITensor}, on::Vector{Vector{Int}}, st::MPS, tr)
    # the site at each position and the position of each site
    at = collect(1:length(st))
    pos = collect(1:length(st))
    for (g, ks) in zip(gates, on)
        ns = sort([ pos[k] for k in ks ])
        if length(ns) == 1
            st[only(ns)] = apply(g, st[only(ns)])
            continue
        end
        st = apply(g, st, ns; move_sites_back = false, tr...)
        # the sites of the gate now next to each other from the first, the others of the window
        # after them, each group in its order
        window = ns[1]:ns[end]
        at[window] = [ at[ns] ; [ at[p] for p in window if !(p in ns) ] ]
        pos[at[window]] = window
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
    if vanished(st)
        return st
    end
    try
        return sites_back(st, gates_applied(gates, on, copy(st), tr), tr)
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
    applied_mpo(mpo, st, tr; alg = "densitymatrix")

the MPS `st` with `mpo` applied by the algorithm `alg`, truncated by `tr`, or the null state
when the result vanishes, see `guarded`, which ``\\langle A\\psi | A\\psi \\rangle`` tells
exactly
"""
applied_mpo(mpo::MPO, st::MPS, tr; alg::String = "densitymatrix") =
    guarded(s -> apply(mpo, s; alg, tr...), s -> iszero(inner(mpo, s, mpo, s)), st)

"""
    gate_terms(terms)

the terms of the superoperator ``\\rho \\mapsto K\\rho K^\\dagger`` of the operator `K` of
terms `terms`
"""
gate_terms(ts::Vector{Term}) = Term[ sided_term(t, u) for t in ts for u in ts ]

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
    bra_leg!(bras, system, i)

the index standing for the bra of site `i` of `system` while `separated` keeps it apart from
the ket, of prime level zero so that the gates of `apply` act on it, made once and kept in
`bras`
"""
bra_leg!(bras::Dict{Int, Index}, system::System, i::Int) =
    get!(bras, i) do
        b, _ = mixer(siteind(Pure, system, i), siteind(Mixed, system, i), system[i])
        noprime(sim(dag(b')))
    end

"""
    separated(system, st, bras, sites)
    joined(system, st, bras, sites)

the MPS `st` of a density matrix with the mixed index of each of `sites` split into its ket and
its bra, the bra on `bras`, and the reverse, the other sites left as they are. A combiner
contracted with the site of a tensor changes the basis of that site only, so the limits of
orthogonality of `st` stay true.
"""
function separated(system::System, st::MPS, bras::Dict{Int, Index}, sites)
    ts = [ st[i] for i in 1:length(st) ]
    for i in sites
        b, c = mixer(siteind(Pure, system, i), siteind(Mixed, system, i), system[i])
        ts[i] = replaceind(ts[i] * dag(c), dag(b'), bra_leg!(bras, system, i))
    end
    return MPS(ts; ortho_lims = ortho_lims(st))
end

function joined(system::System, st::MPS, bras::Dict{Int, Index}, sites)
    ts = [ st[i] for i in 1:length(st) ]
    for i in sites
        b, c = mixer(siteind(Pure, system, i), siteind(Mixed, system, i), system[i])
        ts[i] = replaceind(ts[i], bras[i], dag(b')) * c
    end
    return MPS(ts; ortho_lims = ortho_lims(st))
end

"""
    GateBatch(gates, on, apart, sites)

a step of `apply`: gates applied together, each on its sites `on`, on the indices of the state
or, when `apart`, on the kets and the bras of a density matrix split apart on all their sites
`sites`, see `separated`
"""
struct GateBatch
    gates::Vector{ITensor}
    on::Vector{Vector{Int}}
    apart::Bool
    sites::Vector{Int}
end

"""
    add_gates!(steps, gates, on, apart)

append `gates`, each on its sites `on`, to the steps of `apply`, joining them to the last step
when it is a batch of the same kind
"""
function add_gates!(steps::Vector, gates::Vector{ITensor}, on::Vector{Vector{Int}}, apart::Bool)
    if !isempty(steps) && last(steps) isa GateBatch && last(steps).apart == apart
        append!(last(steps).gates, gates)
        append!(last(steps).on, on)
        union!(last(steps).sites, (k for ks in on for k in ks))
    else
        push!(steps, GateBatch(gates, on, apart, unique(k for ks in on for k in ks)))
    end
    return steps
end

"""
    pure_steps(system, op)

the number a pure state is multiplied by and the steps of `apply` for the operator `op`, from
its rightmost piece to its leftmost, see `pieces`: batches of gates, and the terms of a piece
to apply through its MPO
"""
function pure_steps(system::System, a::IndexedOp{IsOperator})
    coef, ps = pieces(system, a)
    steps = Union{GateBatch, Vector{Term}}[]
    for (p, run) in reverse(ps)
        action, c = piece_action(system, p, run)
        coef *= c
        if coef == 0
            break
        elseif action isa Vector{Term}
            push!(steps, action)
        else
            gates = [ gate_on(m, [ siteind(Pure, system, k) for k in ks ]) for (m, ks) in action ]
            add_gates!(steps, gates, [ collect(ks) for (_, ks) in action ], false)
        end
    end
    return coef, steps
end

"""
    mixed_steps(system, op, bras)

the number a density matrix is multiplied by and the steps of `apply` for the channel or the
generator `op`, from its rightmost factor to its leftmost, see `channel_factors`: a gate as
its operator on the kets and its conjugate on the bras, made in `bras`, see `bra_leg!`, piece by
piece, any other factor as a dense block on the mixed indices, or the terms of either to apply
through their MPO
"""
function mixed_steps(system::System, a::IndexedOp, bras::Dict{Int, Index})
    coef = scalarcoef(a)
    steps = Union{GateBatch, Vector{Term}}[]
    for f in reverse(channel_factors(scalararg(a)))
        if f isa IndexedOp{IsOperator}
            c, ps = pieces(system, f)
            coef *= abs2(c)
            for (p, run) in reverse(ps)
                action, c = piece_action(system, p, run)
                coef *= abs2(c)
                if action isa Vector{Term}
                    push!(steps, gate_terms(action))
                else
                    # K on the kets, its conjugate on the bras
                    gates = [ g for (m, ks) in action
                              for g in (gate_on(m, [ siteind(Pure, system, k) for k in ks ]),
                                        gate_on(conj(m), [ bra_leg!(bras, system, k)
                                                           for k in ks ])) ]
                    on = [ collect(ks) for (_, ks) in action for _ in 1:2 ]
                    add_gates!(steps, gates, on, true)
                end
            end
        else
            terms = lower(system, f)
            laid = [ laid_strings(Mixed, system, t) for t in terms ]
            ks = sort!(unique(k for t in laid for (k, _) in t.blocks))
            if isempty(ks)
                coef *= sum(t.coef for t in terms; init = 0.)
            elseif prod(dim(system[k])^2 for k in ks) ≤ dense_limit
                block = gate_on(mixed_block(system, laid, ks), [ siteind(Mixed, system, k) for k in ks ])
                add_gates!(steps, [block], [ks], false)
            else
                push!(steps, terms)
            end
        end
        if coef == 0
            break
        end
    end
    return coef, steps
end

"""
    applied_steps(R, system, st, steps, tr, bras = nothing)

the MPS `st` of a state of representation `R` with the steps of `apply` applied in order, each
truncated by `tr`, the MPOs built as they come, the gates of a density matrix split apart on
their sites, the bras on `bras`
"""
function applied_steps(::Type{R}, system::System, st::MPS, steps, tr,
                       bras = Dict{Int, Index}()) where R
    for step in steps
        st = if step isa Vector{Term}
            applied_mpo(make_mpo(lay_terms!(PreMPO{R}(system), step, 1)), st, tr)
        elseif step.apart
            kets_bras = separated(system, st, bras, step.sites)
            joined(system, applied_gates(step.gates, step.on, kets_bras, tr), bras, step.sites)
        else
            applied_gates(step.gates, step.on, st, tr)
        end
    end
    return st
end

"""
    apply(op, state; limits = Limits())
    apply(mpo::MPO, state; limits = Limits(), apply_algo = "densitymatrix")

the state with the operator `op`, placed on sites, applied, or the MPO `mpo` from
`make_mpo`. A product acts as the operator it denotes, its rightmost factor first, factor by
factor and never developed; each factor is applied as a dense block on its sites, or through
its MPO when it is too large or lays Jordan-Wigner strings on other sites. On a pure state, a
channel of a single Kraus operator `K`, as `Gate(K)`, is applied as `K`. On a mixed state, a
channel is applied factor by factor in the same way, `Gate(K)` as `K` on the kets and its
conjugate on the bras, an operator `K` being `Gate(K)`; a generator `L` gives `L(ρ)`. The
state is not normalised.

- `limits`: the truncation of the bonds a gate of several sites crosses, or of the whole
  state for an MPO, `cutoff = eps()` by default
- `apply_algo`: for an MPO, the algorithm of its product with the state, `"densitymatrix"` or
  `"naive"`

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
    st = applied_steps(Pure, sys, state.state, steps, truncation(limits))
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
    bras = Dict{Int, Index}()
    coef, steps = mixed_steps(sys, a, bras)
    if coef == 0
        return State(state, 0 * state.state)
    end
    st = applied_steps(Mixed, sys, state.state, steps, truncation(limits), bras)
    return State(state, coef == 1 ? st : coef * st)
end

function apply(mpo::MPO, state::State; limits::Limits = Limits(), apply_algo::String = "densitymatrix")
    check_apply_algo(apply_algo)
    return State(state, applied_mpo(mpo, state.state, truncation(limits); alg = apply_algo))
end
