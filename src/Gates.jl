export apply

"""
    apply(op, ::State; limits::Limits)
    apply(mps, ::State; limits::Limits)
    apply(op, ::Simulation; limits::Limits)

Apply the given gates to the state and truncate the result according to limits.
It is much more efficient to apply all the gates in a single call to apply.

# Examples
    apply(controlled(Z)(1, 3)*H(2)*controlled(X)(3, 4), state)

"""
apply(a::IndexedOp{Pure}, state::State{Mixed}; kwargs...) =
    # prepared before the Gate wrapping: Gate distributes over the product removeMulti
    # leaves behind, down to the one site factors it knows how to lift
    apply(Gate(prepare_gate(a)), state; kwargs...)

function apply(a::IndexedOp{R}, state::State{R}; limits::Limits=Limits()) where R
    check_indices(state.system, a)
    coef, ops = make_ops(state.system, prepare_gate(a))
    st = apply(ops, state.state; move_sites_back_between_gates=false,
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
untouched, and the result is refused if it came out as a sum.

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

make_ops(s::System, a::AtIndex{R, N}) where {R, N} =
    if a == MakeIdentity{R, Indexed, 1}()
        # a factor that contributes no tensor must not leave the gate list untyped:
        # ITensorMPS.product has no method for a Vector{Any}
        (1, ITensor[])
    else
        (1, [ tensor(s, a) ])
    end
    
