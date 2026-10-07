# apply, which applies gates, or an MPO, to a state or a simulation.

export apply

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
    apply(op, ::State; limits::Limits)
    apply(mpo, ::State; limits::Limits, apply_algo)
    apply(op, ::Simulation; limits::Limits)

the state, or the simulation, with the gates `op`, or the MPO `mpo`, from `make_mpo`, applied.
A product of gates acts as the operator it denotes, its rightmost factor first; applying all
the gates in a single call is much more efficient. A pure gate `A` applied to a mixed state
acts as ``\\rho \\mapsto A \\rho A^\\dagger``, see `Gate`. A sum is refused: on a pure state,
apply its MPO, `make_mpo(state, op)`, instead; on a mixed state that MPO is not a gate.

- `limits`: the truncation of each gate of several sites, on the bonds it spans or crosses,
  the other bonds and the gates of one site being left untouched, or of the whole result for
  an MPO (default `Limits()`)
- `apply_algo`: for an MPO, the algorithm of its product with the state, `"densitymatrix"`
  (default) or `"naive"`

# Examples

    apply(controlled(Z)(1, 3) * H(2) * controlled(X)(3, 4), state)
"""
function apply(a::IndexedOp{Pure}, state::State{Mixed}; kwargs...)
    # before any rewriting, so that a site out of the system is named as the caller wrote it
    check_indices(state.system, a)
    # prepared first: build_gate distributes over the product removeMulti leaves, down to the
    # one site factors it knows how to lift
    return apply(build_gate(prepare_gate(a)), state; kwargs...)
end

apply(a::IndexedOp{Mixed}, ::State{Pure}; kwargs...) =
    error("$a acts on a density matrix, which a pure state is not: apply it to mix(state)")

function apply(a::IndexedOp{R}, state::State{R}; limits::Limits=Limits()) where R
    check_indices(state.system, a)
    coef, ops = make_ops(state.system, prepare_gate(a))
    # ITensorMPS applies a list of gates first to last, and the factors of a product act
    # right to left: A*B is B applied first
    st = apply(reverse(ops), state.state; move_sites_back_between_gates=false,
            limits.cutoff, limits.maxdim, limits.mindim)
    # the coefficient is carried apart: a gate made only of identities places no tensor
    return State(state, coef == 1 ? st : coef * st)
end

function apply(mpo::MPO, state::State; limits::Limits=Limits(),
               apply_algo::String = "densitymatrix")
    check_apply_algo(apply_algo)
    return State(state, apply(mpo, state.state; alg = apply_algo, limits.cutoff, limits.maxdim,
                              limits.mindim))
end
    
"""
    parity_sign

the gate ``(-1)^{p_1 p_2}`` of two sites, `p` being the fermionic parity of a site, which gives
a function of fermionic operators the strings of the sites it skips, see `strung_function`
"""
const parity_sign = Operator{2}("ParitySign", Id ⊗ Id - 2 * (((Id - F) / 2) ⊗ ((Id - F) / 2)),
                                involution_op)

"""
    order_signs(sites, index)

the diagonal, as a vector, of the matrix that takes an operator of the `sites`, placed on the
sites `index`, from the order of its factors to the order of the sites: the sign of each pair
of odd basis states whose sites it swaps, a basis written with the last site varying fastest
"""
function order_signs(sites, index)
    ps = [ (1 .- real(diag(matrix(F, s)))) ./ 2 for s in sites ]
    d = ones(prod(length, ps))
    for a in eachindex(index), b in a+1:length(index)
        if index[a] > index[b]
            d .*= 1 .- 2 .* foldl(kron, [ k == a || k == b ? ps[k] : ones(length(ps[k]))
                                          for k in eachindex(ps) ])
        end
    end
    return d
end

"""
    laid_with_strings(a, index)

the operator `a` of several sites as an operator of the sites `index` in their order, see
`order_signs`, and the product `K` of the parity signs that gives it the strings of the sites
in between, see `strung_function`
"""
function laid_with_strings(a, index)
    function laid(sites...)
        d = order_signs(sites, index)
        return d .* matrix(a, sites...) .* transpose(d)
    end
    g = Operator{length(index)}(repr(a), laid, plain_op)
    k = ProdOp(IndexedOp{Pure}[ parity_sign(s, j) for j in min(index...)+1:max(index...)-1
                                if j ∉ index for s in index if s < j ])
    return g, k
end

"""
    strung_function(a, index)

the gate of the function `a` of an even fermionic operator, placed on the sites `index`: its
matrix laid in the order of the sites, see `order_signs`, conjugated by
``K = \\prod (-1)^{p_s p_k}`` over the sites `s` of the gate and the sites `k` in between on
their right, which gives each odd factor the strings of the sites in between on its right,
equivalent for an even operator to those on its left. A function of an operator of odd or no
definite parity is left whole, for `prepare_gate` to refuse.
"""
function strung_function(a, index)
    if fermion_parity(a, false) ≠ 0
        return a(index...)
    end
    g, k = laid_with_strings(a, index)
    return k * g(index...) * k
end

"""
    expand_gate(op)

the gate `op` with each factor of several sites holding a fermionic operator replaced, where
its definition allows, by a product of factors of one site in the same order, `(A ⊗ B)(i, j)`
being `A(i) * B(j)`, and a function of an even fermionic operator by its gate with its
strings, see `strung_function`. Any other factor of several sites is kept whole.
"""
expand_gate(a::ProdOp{R, Indexed, 1}) where R = ProdOp(IndexedOp{R}[ expand_gate(x) for x in a.subs ])
expand_gate(a::SumOp{R, Indexed, 1}) where R = SumOp(IndexedOp{R}[ expand_gate(x) for x in a.subs ])
expand_gate(a::ScalarOp{R, Indexed, 1}) where R = a.coef * expand_gate(a.arg)
expand_gate(a::AtIndex{R, N}) where {R, N} =
    N > 1 && hasfermionic(a.op) ? expand_placed(a.op, a.index) : a
expand_gate(a::IndexedOp) = a

expand_placed(a::TensorOp{R}, index) where R =
    ProdOp(IndexedOp{R}[ expand_gate(o(index[p]...)) for (o, p) in zip(a.subs, factor_sites(a)) ])
expand_placed(a::ProdOp{R}, index) where R =
    ProdOp(IndexedOp{R}[ expand_gate(o(index...)) for o in a.subs ])
expand_placed(a::SumOp{R}, index) where R =
    SumOp(IndexedOp{R}[ expand_gate(o(index...)) for o in a.subs ])
expand_placed(a::IntPowOp, index) = ProdOp(fill(expand_gate(a.arg(index...)), a.expo))
expand_placed(a::Union{ExpOp, GenPowOp{Pure}, ModOp, Proj}, index) = strung_function(a, index)
expand_placed(a::Operator, index) = a.expr isa Op ? expand_gate(a.expr(index...)) : a(index...)
expand_placed(a::DagOp, index) = placed_dag(expand_gate(a.arg(index...)))
expand_placed(a::Left, index) = sided(Left, expand_gate(a.arg(index...)))
expand_placed(a::Right, index) = sided(Right, expand_gate(a.arg(index...)))
function expand_placed(a::Gate, index)
    e = expand_gate(a.arg(index...))
    return ProdOp([sided(Left, e), sided(Right, e)])
end
# Σ Gate(K Π K) = Gate(K) Σ Gate(Π) Gate(K), the same K giving every projector its strings
function expand_placed(a::Dephase, index)
    g, k = laid_with_strings(a.arg, index)
    gk = build_gate(k)
    return gk * Dephase(g)(index...) * gk
end
expand_placed(a, index) = a(index...)

"""
    placed_dag(op)

the adjoint of an operator on pure states placed on sites, the factors of a product taken in
the reverse order, each one keeping its sites
"""
placed_dag(a::ProdOp{Pure, Indexed, 1}) = ProdOp(reverse(placed_dag.(a.subs)))
placed_dag(a::SumOp{Pure, Indexed, 1}) = SumOp(placed_dag.(a.subs))
placed_dag(a::ScalarOp{Pure, Indexed, 1}) = conj(a.coef) * placed_dag(a.arg)
placed_dag(a::AtIndex{Pure}) = dag(a.op)(a.index...)
placed_dag(a::IdentityOp) = a

"""
    spans_sites(op)

whether an operator placed on sites holds a factor acting on several sites at once
"""
spans_sites(a::AtIndex) = length(a.index) > 1
spans_sites(a::Union{ProdOp, SumOp}) = any(spans_sites, a.subs)
spans_sites(a::ScalarOp) = spans_sites(a.arg)
spans_sites(::Op) = false

"""
    prepare_gate(op)

the gate `op` with the Jordan-Wigner strings of its fermionic factors inserted and spelled out
as one factor per site, see `removeMulti`. A gate with a fermionic factor is expanded, see
`expand_gate`, then simplified piece by piece, each run of factors of one site on its own and
the factors of several sites left whole, which `simplify` would rewrite into sums or move
strings past. The gate is refused if a piece becomes a sum, or if a function of an odd
fermionic operator remains, which mixes the two parities.
"""
function prepare_gate(a::IndexedOp{R}) where R
    if !hasfermionic(a)
        return a
    end
    e = expand_gate(a)
    pieces = IndexedOp{R}[]
    run = IndexedOp{R}[]
    function close_run()
        p = simplify(ProdOp(run))
        if scalararg(p) isa SumOp
            error("cannot apply $a as a gate: inserting its Jordan-Wigner strings makes " *
                  "it a sum, which apply cannot place. Use make_mpo to build an MPO instead")
        end
        push!(pieces, p)
        empty!(run)
    end
    for x in prodsubs(e)
        if !spans_sites(x)
            push!(run, x)
        elseif x isa AtIndex
            close_run()
            push!(pieces, x)
        else
            error("cannot apply sums as gates ($a)")
        end
    end
    close_run()
    b = removeMulti(scalarcoef(e) * ProdOp(pieces))
    if hasfermionic(b)
        error("cannot apply $a as a gate: it holds a function of an odd fermionic operator, " *
              "which mixes the two parities")
    end
    return b
end

"""
    make_ops(::System, op)

the coefficient of the gate `op` and the tensors to place for it, one per factor, a sum being
refused
"""
make_ops(::System, a::SumOp) =
    error("cannot apply sums as gates ($a)")

make_ops(::System, a::ComOp) =
    error("cannot apply sums as gates ($a is a sum gathered by compact)")

make_ops(::System, a::Evolver) =
    error("cannot apply sums as gates ($a is Left + Right of its argument)")

# a null gate, `0Id`, places no tensor and makes the state null through its coefficient
function make_ops(s::System, a::ScalarOp)
    coef, ops = make_ops(s, a.arg)
    return (coef * a.coef, ops)
end

function make_ops(s::System, a::ProdOp)
    coef = 1
    ops = ITensor[]
    for x in a.subs
        c, o = make_ops(s, x)
        coef *= c
        append!(ops, o)
    end
    return (coef, ops)
end

make_ops(s::System, a::AtIndex) = (1, [ tensor(s, a) ])

# a factor that contributes no tensor must not leave the gate list untyped: ITensorMPS.product
# has no method for a Vector{Any}
make_ops(::System, ::IdentityOp) = (1, ITensor[])
