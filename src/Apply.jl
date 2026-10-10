# apply on pure states: a product of operators applied factor by factor, right to left, never
# developed, each piece as a dense block on its sites, as one-site gates, or through its MPO.

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
    block_gate(system, m, sites)

the tensor of the gate of matrix `m` on `sites`
"""
function block_gate(system::System, m::AbstractMatrix, sites)
    idx = [ siteind(Pure, system, k) for k in sites ]
    a, inds = on_legs(Matrix(m), [ i' for i in idx ], [ dag(i) for i in idx ])
    return charged_itensor(a, inds)
end

"""
    term_gates(system, t)

the gates of one site of the term `t`, its coefficient left out
"""
term_gates(system::System, t::Term) = [ block_gate(system, m, [k]) for (k, m) in t.blocks ]

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

how `apply` applies the piece `p`: `(gates, number)`, the gates and the number they bring, or
`(mpo, number)`. A run of factors of one site is a single term, applied as gates of one site.
Any other piece on at most two sites, or on sites of dimensions of product at most
`dense_limit`, laying no string elsewhere, is a dense block; otherwise it is developed, and
goes through its MPO, or through gates of one site when it is a single term.
"""
function piece_action(system::System, p::IndexedOp{IsOperator}, run::Bool)
    sites = sort(sites_of(p))
    if !run && (length(sites) ≤ 2 || prod(dim(system[k]) for k in sites) ≤ dense_limit)
        m = piece_matrix(system, p, sites)
        if stays_on(system, m, sites)
            return [ block_gate(system, m, sites) ], 1
        end
    end
    terms = lower(system, p)
    if isempty(terms)
        return ITensor[], 0
    elseif length(terms) == 1
        t = only(terms)
        return term_gates(system, t), t.coef
    end
    pre = lay_terms!(PreMPO{Pure}(system), terms, 1)
    return make_mpo(pre), 1
end

"""
    apply(op, state::State{Pure}; limits = Limits())
    apply(mpo::MPO, state; limits = Limits(), apply_algo = "densitymatrix")

the pure state with the operator `op`, placed on sites, applied, or the MPO `mpo` from
`make_mpo`. A product acts as the operator it denotes, its rightmost factor first, factor by
factor and never developed; each factor is applied as a dense block on its sites, or through
its MPO when it is too large or lays Jordan-Wigner strings on other sites. The state is not
normalised.

- `limits`: the truncation of the bonds a gate of several sites crosses, or of the whole
  state for an MPO, `cutoff = eps()` by default
- `apply_algo`: for an MPO, the algorithm of its product with the state, `"densitymatrix"` or
  `"naive"`

# Examples

    apply(controlled(Z)(1, 3) * H(2) * controlled(X)(3, 4), state)
    apply(exp(-0.1im * (X ⊗ X + Y ⊗ Y))(2, 5), state)
"""
function apply(a::IndexedOp{IsOperator}, state::State{Pure}; limits::Limits = Limits())
    sys = state.system
    check_indices(sys, a)
    tr = truncation(limits)
    coef, ps = pieces(sys, a)
    st = state.state
    gates = ITensor[]
    flush!() =
        if !isempty(gates)
            st = apply(gates, st; move_sites_back_between_gates = false, tr...)
            empty!(gates)
        end
    for (p, run) in reverse(ps)
        action, c = piece_action(sys, p, run)
        coef *= c
        if coef == 0
            return State(state, 0 * state.state)
        elseif action isa MPO
            flush!()
            st = apply(action, st; alg = "densitymatrix", tr...)
        else
            append!(gates, action)
        end
    end
    flush!()
    return State(state, coef == 1 ? st : coef * st)
end

function apply(mpo::MPO, state::State; limits::Limits = Limits(), apply_algo::String = "densitymatrix")
    check_apply_algo(apply_algo)
    return State(state, apply(mpo, state.state; alg = apply_algo, truncation(limits)...))
end
