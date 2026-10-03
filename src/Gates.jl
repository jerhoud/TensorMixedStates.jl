# apply, which applies gates, or an MPO, to a state or a simulation.

export apply, kraus_operators

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

apply(a::IndexedOp{Mixed}, ::State{Pure}; kwargs...) =
    error("$a acts on a density matrix, which a pure state is not: apply it to mix(state)")

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
    strung_function(a, index)

the gate of the function `a` of an even fermionic operator, placed on the sites `index`. Its
matrix on these sites, laid in the order of the sites with the sign of each pair of odd factors
it swaps, see `order_signs`, is the gate on the sites taken as neighbours. The strings its odd
factors take through the sites in between come from the conjugation by
``K = \\prod (-1)^{p_s p_k}``, over the sites `s` of the gate and the sites `k` in between on
their right: `K` is diagonal and squares to the identity, conjugating an operator by it gives
each odd factor the strings of the sites in between on its right, which an even operator cannot
tell from those on its left, and a function goes through the conjugation. A function of an odd
operator, or of one of no definite parity, mixes the two parities and is left whole, for
`prepare_gate` to refuse.
"""
function strung_function(a, index)
    if fermion_parity(a, false) ≠ 0
        return a(index...)
    end
    n = length(index)
    function laid(sites...)
        d = order_signs(sites, index)
        return d .* matrix(a, sites...) .* transpose(d)
    end
    g = Operator{n}(repr(a), laid, plain_op)
    k = ProdOp(IndexedOp{Pure}[ parity_sign(s, j) for j in min(index...)+1:max(index...)-1
                                if j ∉ index for s in index if s < j ])
    return k * g(index...) * k
end

"""
    expand_gate(op)

the gate `op` with each factor of several sites holding a fermionic operator replaced by a
product of factors of one site, by the definitions of the tensor product, `(A ⊗ B)(i, j)` being
`A(i) * B(j)`, of the product, of the integer power, of the adjoint, of an operator defined by
an expression, and of `Left`, `Right` and `Gate`, the order of the factors being kept. A
function of an even fermionic operator becomes its gate with its strings, see
`strung_function`. Any other factor of several sites is kept whole: a function of an odd
fermionic operator, which `prepare_gate` then refuses, or an operator with no fermionic factor,
which `apply` places as it is.
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
expand_placed(a::Union{ExpOp, GenPowOp{Pure}, ModOp}, index) = strung_function(a, index)
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
gate is refused if a piece becomes a sum, or if a function of an odd fermionic operator
remains whole, on several sites or on the first one, where it has no string to take but mixes
the two parities.
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
        error("cannot apply $a as a gate: it holds a function of an odd fermionic operator, " *
              "which mixes the two parities")
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
    


"""
    reset_kraus(system, a, i)

the Kraus operators of the superoperator `SetState` `a` on site `i` of `system`: one
``|s\\rangle\\langle k|`` for each basis state ``k`` of the site, for every eigenvector ``s`` of
its density matrix, weighted by the square root of its eigenvalue, for a mixed state
"""
function reset_kraus(system::System, a::SetState, i::Int)
    site = system[i]
    if matrix(F, site) != I
        error("kraus_operators does not take SetState on the fermionic site $i: its Kraus " *
              "operators would move fermions, which the reset does not with strings")
    end
    s = state(site, a.state)
    d = dim(site)
    if s isa AbstractVector
        vs = [ (1., s) ]
    else
        e = eigen(Hermitian(s))
        vs = [ (p, v) for (p, v) in zip(e.values, eachcol(e.vectors)) if p > rounding_tol ]
    end
    return IndexedOp{Pure}[ Operator{1}("Reset", sqrt(p) * v * [ j == k ? 1. : 0. for j in 1:d ]',
                                        plain_op)(i)
                            for (p, v) in vs for k in 1:d ]
end

"""
    kraus_channel(system, a)

the Kraus operators of the factor `a` of a product of gates on mixed states, a sum of terms of
weights not negative, each a `SetState`, gates of a pure operator, `Gate(A)`, placed on sites,
or a sum of them placed on the same sites, a noisy gate
"""
function kraus_channel(system::System, a::IndexedOp{Mixed})
    terms = IndexedOp{Pure}[]
    for t in sumsubs(a)
        c, b = scalarcoef(t), scalararg(t)
        if !isreal(c) || real(c) < 0
            error("a gate of weight $c has no Kraus operator: the weight must be real and not negative")
        end
        w = sqrt(real(c))
        if b isa IdentityOp
            push!(terms, w * IdentityOp{Pure, Indexed, 1}())
        elseif b isa AtIndex && b.op isa SetState
            append!(terms, w .* reset_kraus(system, b.op, only(b.index)))
        elseif b isa AtIndex && b.op isa SumOp
            # c₁ Gate(A₁) + c₂ Gate(A₂) placed on the same sites, as a noisy gate
            for u in sumsubs(b.op)
                append!(terms, w .* kraus_channel(system, scalarcoef(u) * AtIndex(scalararg(u), b.index)))
            end
        else
            # gates on several sites at once, as Gate(X)(1) * Gate(X)(2)
            k = w * IdentityOp{Pure, Indexed, 1}()
            for f in prodsubs(b)
                if !(f isa AtIndex && f.op isa Gate)
                    error("kraus_operators reads gates, sums of them and SetState, and $a holds $f")
                end
                k = k * f.op.arg(f.index...)
            end
            push!(terms, k)
        end
    end
    return terms
end

"""
    kraus_operators(system, gates)

the channels of the product of gates `gates` on `system`, in the order they act, the rightmost
first, each the vector of its Kraus operators, operators on pure states placed on sites: the
channel ``\\rho \\mapsto \\sum_k K_k \\rho K_k^\\dagger``. A gate on pure states, as `X(1)`, is a
channel of a single operator, a noisy gate `c₁ Gate(A₁)(i) + c₂ Gate(A₂)(i)` one of
`sqrt(c₁) * A₁(i)` and `sqrt(c₂) * A₂(i)`, and `SetState(s)(i)` one of an operator
``|s\\rangle\\langle k|`` for each basis state ``k`` of site `i`. A coefficient of the whole
product goes to its first channel.

What it reads is what a representation of a state needs to apply the gates of a `Gates`
phase by sampling their Kraus operators, quantum trajectories for instance. A superoperator of
another form, as `Left(X)(1)`, or a gate of negative or complex weight is refused.

# Examples

    kraus_operators(System(2, Qubit()), (0.9Gate(Id) + 0.1Gate(X))(1) * Gate(H)(2))
"""
function kraus_operators(::System, gates::IndexedOp{Pure})
    c = scalarcoef(gates)
    channels = [ IndexedOp{Pure}[f] for f in reverse(prodsubs(scalararg(gates)))
                 if !(f isa IdentityOp) ]
    if isempty(channels)
        return [ IndexedOp{Pure}[ c * IdentityOp{Pure, Indexed, 1}() ] ]
    end
    channels[1] = c .* channels[1]
    return channels
end

function kraus_operators(system::System, gates::IndexedOp{Mixed})
    c = scalarcoef(gates)
    if !isreal(c) || real(c) < 0
        error("gates of weight $c have no Kraus operator: the weight must be real and not negative")
    end
    channels = [ kraus_channel(system, f) for f in reverse(prodsubs(scalararg(gates)))
                 if !(f isa IdentityOp) ]
    if isempty(channels)
        return [ IndexedOp{Pure}[ sqrt(real(c)) * IdentityOp{Pure, Indexed, 1}() ] ]
    end
    channels[1] = sqrt(real(c)) .* channels[1]
    return channels
end
