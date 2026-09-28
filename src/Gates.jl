export apply

"""
    apply(op, ::State; limits::Limits)
    apply(mpo, ::State; limits::Limits)
    apply(op, ::Simulation; limits::Limits)

Apply the given gates to the state. `limits` constrain the truncations made while a gate of
several sites is applied, on the bond it spans and on those crossed to bring its sites
together: a gate of one site, and the other bonds, are not truncated. An MPO truncates the
whole result. It is much more efficient to apply all the gates in a single call to apply.
A product of gates is the operator it denotes, its rightmost factor acting first.

# Examples
    apply(controlled(Z)(1, 3)*H(2)*controlled(X)(3, 4), state)

"""
function apply(a::IndexedOp{Pure}, state::State{Mixed}; kwargs...)
    # on the operator as it was written, as the pure path does, so that a site out of the
    # system is named as the caller wrote it and not as a factor of its string
    check_indices(state.system, a)
    # prepared before the Gate wrapping: Gate distributes over the product removeMulti
    # leaves behind, down to the one site factors it knows how to lift
    return apply(Gate(prepare_gate(a)), state; kwargs...)
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

apply(mpo::MPO, state::State; limits::Limits=Limits()) =
    State(state, apply(mpo, state.state; limits.cutoff, limits.maxdim, limits.mindim))
    
"""
    prepare_gate(op)

`apply` places one local tensor per factor and has no way to build the Jordan-Wigner
string a fermionic operator needs: only `simplify` inserts those. Simplifying every gate
is not an option, since it replaces a gate defined by an expression, such as `Swap`, with
that expression, and a product of those becomes a sum `apply` cannot place. So only an
operator that still has a fermionic factor is simplified, which leaves every other gate
untouched, and the result is refused if it came out as a sum. It is refused as well if a
fermionic factor survived, which happens inside an exponential or a power of several sites:
`simplify` keeps those whole and has no string to put into them.

`removeMulti` then spells the string out as one factor per site, which is what `PreMPO`
does too. Those one site factors are built by the `Multi_F` constructor, which is where
the knowledge of whether the string acts on the left of the density matrix, on its right,
or on both, already lives.
"""
prepare_gate(a) =
    if has_fermionic(a)
        b = removeMulti(simplify(a))
        if scalararg(b) isa SumOp
            error("cannot apply $a as a gate: inserting its Jordan-Wigner strings makes " *
                  "it a sum, which apply cannot place. Use make_mpo to build an MPO instead")
        end
        if has_fermionic(b)
            error("cannot apply $a as a gate: its Jordan-Wigner strings cannot be inserted " *
                  "into a function of a fermionic operator acting on several sites")
        end
        b
    else
        a
    end

"""
    make_ops(::System, op)

the coefficient of a gate and the tensors to place for it, one per factor.

The two are kept apart because a factor may place no tensor: an identity contributes nothing,
and a gate made of nothing else leaves an empty list, which no coefficient can ride.
"""
make_ops(::System, a::SumOp) =
    error("cannot apply sums as gates ($a)")

# Left(H) + Right(H), a sum
make_ops(::System, a::Evolver) =
    error("cannot apply sums as gates ($a is Left + Right of its argument)")

make_ops(s::System, a::ScalarOp) =
    if a.coef == 0
        error("cannot apply null gate")
    else
        coef, ops = make_ops(s, a.arg)
        (coef * a.coef, ops)
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
    
