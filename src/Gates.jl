# apply, which applies gates, or an MPO, to a state or a simulation.

export apply

"""
    check_apply_algo(apply_algo)

refuse an algorithm of the product of an MPO by a state that is not one of those offered,
`"densitymatrix"` and `"naive"`. ITensorMPS has `"fit"` too, which needs a number of sweeps of
its own, and `"zipup"`, which it only has for a state from version 0.3.45 on, above the lowest
version TMS accepts.
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

the state, or the simulation, with the gates `op`, or the MPO `mpo`, applied.

A product of gates is the operator it denotes, its rightmost factor acting first, and applying
all the gates in a single call is much more efficient. A pure gate `A` applied to a mixed
state acts as ``\\rho \\mapsto A \\rho A^\\dagger``, see `Gate`. A sum is refused: on a pure
state, apply its MPO, `make_mpo(state, op)`, instead. On a mixed state that MPO is not a gate,
see `make_mpo`.

`limits` (default `Limits()`, no truncation) constrains the truncations made while a gate of
several sites is applied, on the bond it spans and on those crossed to bring its sites
together; a gate of one site, and the other bonds, are not truncated. An MPO truncates the
whole result.

`apply_algo`, for an MPO only, is the algorithm of its product with the state, as
`ITensorMPS.apply` takes it: `"densitymatrix"` (default) or `"naive"`.

# Examples

    apply(controlled(Z)(1, 3) * H(2) * controlled(X)(3, 4), state)
"""
function apply(a::IndexedOp{Pure}, state::State{Mixed}; kwargs...)
    # on the operator as it was written, as the pure path does, so that a site out of the
    # system is named as the caller wrote it and not as a factor of its string
    check_indices(state.system, a)
    # prepared before the gate is built: build_gate distributes over the product removeMulti
    # leaves behind, down to the one site factors it knows how to lift
    return apply(build_gate(prepare_gate(a)), state; kwargs...)
end

function apply(a::IndexedOp{R}, state::State{R}; limits::Limits=Limits()) where R
    check_indices(state.system, a)
    coef, ops = make_ops(state.system, prepare_gate(a))
    # ITensorMPS applies a list of gates first to last, and the factors of a product act
    # right to left: A*B is B applied first
    st = apply(reverse(ops), state.state; move_sites_back_between_gates=false,
            limits.cutoff, limits.maxdim, limits.mindim)
    # the coefficient is carried here rather than laid on the first tensor, because a gate
    # whose factors are all identities places no tensor at all and there would be nothing
    # to lay it on: `2Id(1)` then went through leaving the state unscaled
    return State(state, coef == 1 ? st : coef * st)
end

function apply(mpo::MPO, state::State; limits::Limits=Limits(),
               apply_algo::String = "densitymatrix")
    check_apply_algo(apply_algo)
    return State(state, apply(mpo, state.state; alg = apply_algo, limits.cutoff, limits.maxdim,
                              limits.mindim))
end
    
"""
    expand_gate(op)

the gate `op` with each factor of several sites holding a fermionic operator replaced by a
product of factors of one site, by the definitions of the tensor product, `(A ⊗ B)(i, j)` being
`A(i) * B(j)`, of the product, of the integer power, of the adjoint, of an operator defined by
an expression, and of `Left`, `Right` and `Gate`, the order of the factors being kept. Any
other factor of several sites is kept whole: a function of a fermionic operator, which
`prepare_gate` then refuses, or an operator with no fermionic factor, which `apply` places as
it is.
"""
expand_gate(a::ProdOp{R, Indexed, 1}) where R = ProdOp(IndexedOp{R}[ expand_gate(x) for x in a.subs ])
expand_gate(a::SumOp{R, Indexed, 1}) where R = SumOp(IndexedOp{R}[ expand_gate(x) for x in a.subs ])
expand_gate(a::ScalarOp{R, Indexed, 1}) where R = a.coef * expand_gate(a.arg)
expand_gate(a::AtIndex{R, N}) where {R, N} =
    N > 1 && has_fermionic(a.op) ? expand_placed(a.op, a.index) : a
expand_gate(a::IndexedOp) = a

expand_placed(a::TensorOp, index) =
    ProdOp(IndexedOp{Pure}[ expand_gate(o(index[p]...)) for (o, p) in zip(a.subs, factor_sites(a)) ])
expand_placed(a::ProdOp{R}, index) where R =
    ProdOp(IndexedOp{R}[ expand_gate(o(index...)) for o in a.subs ])
expand_placed(a::SumOp{R}, index) where R =
    SumOp(IndexedOp{R}[ expand_gate(o(index...)) for o in a.subs ])
expand_placed(a::IntPowOp, index) = ProdOp(fill(expand_gate(a.arg(index...)), a.expo))
expand_placed(a::Operator, index) = a.expr isa Op ? expand_gate(a.expr(index...)) : a(index...)
expand_placed(a::DagOp, index) = placed_dag(expand_gate(a.arg(index...)))
expand_placed(a::Left, index) = sided(Left, expand_gate(a.arg(index...)))
expand_placed(a::Right, index) = sided(Right, expand_gate(a.arg(index...)))
function expand_placed(a::Gate, index)
    e = expand_gate(a.arg(index...))
    return ProdOp([sided(Left, e), sided(Right, e)])
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
as one factor per site, see `removeMulti`, as `PreMPO` does.

`apply` places one tensor per factor and cannot build a string: only `simplify` inserts them.
Simplifying a whole gate is not an option. It would replace a gate defined by an expression,
such as `Swap`, with that expression, and a product of those becomes a sum `apply` cannot
place. And it sorts the factors by site, moving strings past a factor of several sites, which
is only right when the string holds all of its sites or none of them. So a gate with a
fermionic factor is first expanded into factors of one site where it can be, see
`expand_gate`, and then simplified piece by piece: each run of factors of one site on its own,
the factors of several sites left whole between them, so that `simplify` never sees one. The
gate is refused if a piece becomes a sum, or if a fermionic factor remains inside a function of
several sites.
"""
function prepare_gate(a::IndexedOp{R}) where R
    if !has_fermionic(a)
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
    if has_fermionic(b)
        error("cannot apply $a as a gate: its Jordan-Wigner strings cannot be inserted " *
              "into a function of a fermionic operator acting on several sites")
    end
    return b
end

"""
    make_ops(::System, op)

the coefficient of the gate `op` and the tensors to place for it, one per factor, a sum being
refused. The two are kept apart because an identity places no tensor, so a gate made only of
identities has no tensor to carry the coefficient.
"""
make_ops(::System, a::SumOp) =
    error("cannot apply sums as gates ($a)")

make_ops(::System, a::ComOp) =
    error("cannot apply sums as gates ($a is a sum gathered by compact)")

# Left(H) + Right(H), a sum
make_ops(::System, a::Evolver) =
    error("cannot apply sums as gates ($a is Left + Right of its argument)")

# a null gate is `0Id`, which places no tensor and makes the state null, as a gate that
# annihilates it does
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
    
