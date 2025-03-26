export apply

"""
    apply(op, ::State; limits::Limits)
    apply(op, ::Simulation; limits::Limits)

Apply the given gates to the state and truncate the result according to maxdim and cutoff.
It is much more efficient to apply all the gates in a single call to apply.

# Examples
    apply(CZ(1,3)*H(2)*CNOT(3,4), state; limits)

"""
function apply(a::ExprIndexed, state::State; limits::Limits=Limits())
    ops = make_ops(state.type, state.system, a)
    st = apply(ops, state.state; move_sites_back_between_gates=false, limits.cutoff, limits.maxdim, limits.mindim)
    return State(state, st)
end

make_ops(::PM, ::System, ::SumOp{IndexOp}) = error("cannot apply sums as gates (sums of operator are possible)")

function make_ops(tp::PM, s::System, a::ProdOp{IndexOp})
    r = map(a.subs) do b
        make_ops(tp, s, b)
    end
    r = vcat(r...)
    if !isempty(r)
        r[1] *= a.coef
    end
    return r
end

make_ops(::Mixed, s::System, a::Gate{IndexOp}) = make_ops(Mixed(), s, a.arg)
make_ops(::T, s::System, a::Indexed{T, IndexOp}) where T = [ tensor(s, a) ]
make_ops(::Mixed, s::System, a::Indexed{Pure, IndexOp}) = [ tensor(s, Gate(a)) ]