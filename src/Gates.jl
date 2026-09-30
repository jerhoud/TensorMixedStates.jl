# apply, which applies gates, or an MPO, to a state or a simulation.

export apply

"""
    apply(op, ::State; limits::Limits)
    apply(mpo, ::State; limits::Limits)
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

# Examples

    apply(controlled(Z)(1, 3) * H(2) * controlled(X)(3, 4), state)
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

the gate `op` with the Jordan-Wigner strings of its fermionic factors inserted and spelled out
as one factor per site, see `removeMulti`, as `PreMPO` does.

`apply` places one tensor per factor and cannot build a string: only `simplify` inserts them.
Simplifying every gate is not an option, since it replaces a gate defined by an expression,
such as `Swap`, with that expression, and a product of those becomes a sum `apply` cannot
place. So only a gate with a fermionic factor is simplified, and it is refused if it becomes a
sum, or if a fermionic factor remains inside a function of several sites, which `simplify`
keeps whole.
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

the coefficient of the gate `op` and the tensors to place for it, one per factor, a sum being
refused. The two are kept apart because an identity places no tensor, so a gate made only of
identities has no tensor to carry the coefficient.
"""
make_ops(::System, a::SumOp) =
    error("cannot apply sums as gates ($a)")

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
    
