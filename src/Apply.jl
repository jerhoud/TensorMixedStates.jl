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
        t = only(terms)
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
    applied_gates(gates, st, tr)

the MPS `st` with `gates` applied, truncated by `tr`, or the null state when the result
vanishes, see `guarded`: on a failure on charged indices, the gates are applied again one by
one, so that the one taking the state to zero is found
"""
function applied_gates(gates::Vector{ITensor}, st::MPS, tr)
    if vanished(st)
        return st
    end
    try
        return apply(gates, copy(st); move_sites_back_between_gates = false, tr...)
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
    bra_legs(system)

an index for the bra of each site, standing in for the one its mixed index combines while
`separated` keeps them apart: of prime level zero, so that the gates of `apply` act on it
"""
function bra_legs(system::System)
    return map(1:length(system)) do i
        b, _ = mixer(siteind(Pure, system, i), siteind(Mixed, system, i), system[i])
        noprime(sim(dag(b')))
    end
end

"""
    separated(system, st, bras)
    joined(system, st, bras)

the MPS `st` of a density matrix with the mixed index of each site split into its ket and its
bra, the bra on `bras`, and the reverse
"""
function separated(system::System, st::MPS, bras)
    return MPS([ begin
                     b, c = mixer(siteind(Pure, system, i), siteind(Mixed, system, i), system[i])
                     replaceind(st[i] * dag(c), dag(b'), bras[i])
                 end for i in 1:length(system) ])
end

function joined(system::System, st::MPS, bras)
    return MPS([ begin
                     b, c = mixer(siteind(Pure, system, i), siteind(Mixed, system, i), system[i])
                     replaceind(st[i], bras[i], dag(b')) * c
                 end for i in 1:length(system) ])
end

"""
    GateBatch(gates, apart)

a step of `apply`: gates applied together, on the indices of the state or, when `apart`, on
the kets and the bras of a density matrix split apart, see `separated`
"""
struct GateBatch
    gates::Vector{ITensor}
    apart::Bool
end

"""
    add_gates!(steps, gates, apart)

append `gates` to the steps of `apply`, joining them to the last step when it is a batch of the
same kind
"""
function add_gates!(steps::Vector, gates::Vector{ITensor}, apart::Bool)
    if !isempty(steps) && last(steps) isa GateBatch && last(steps).apart == apart
        append!(last(steps).gates, gates)
    else
        push!(steps, GateBatch(gates, apart))
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
            add_gates!(steps, gates, false)
        end
    end
    return coef, steps
end

"""
    mixed_steps(system, op, bras)

the number a density matrix is multiplied by and the steps of `apply` for the channel or the
generator `op`, from its rightmost factor to its leftmost, see `channel_factors`: a gate as
its operator on the kets and its conjugate on the bras `bras`, piece by piece, any other factor
as a dense block on the mixed indices, or the terms of either to apply through their MPO
"""
function mixed_steps(system::System, a::IndexedOp, bras)
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
                                        gate_on(conj(m), bras[ks])) ]
                    add_gates!(steps, gates, true)
                end
            end
        else
            terms = lower(system, f)
            ks = sort!(unique(k for t in terms for (k, _) in t.blocks))
            if isempty(ks)
                coef *= sum(t.coef for t in terms; init = 0.)
            elseif prod(dim(system[k])^2 for k in ks) ≤ dense_limit
                add_gates!(steps, [ gate_on(mixed_block(system, terms, ks),
                                            [ siteind(Mixed, system, k) for k in ks ]) ], false)
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
the bras `bras`
"""
function applied_steps(::Type{R}, system::System, st::MPS, steps, tr, bras = nothing) where R
    for step in steps
        st = if step isa Vector{Term}
            applied_mpo(make_mpo(lay_terms!(PreMPO{R}(system), step, 1)), st, tr)
        elseif step.apart
            joined(system, applied_gates(step.gates, separated(system, st, bras), tr), bras)
        else
            applied_gates(step.gates, st, tr)
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
    bras = bra_legs(sys)
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
