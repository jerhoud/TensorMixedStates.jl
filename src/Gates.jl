export apply

make_ops(type::PM, system::System, a::Indexed) = tensor(type; system, a)

function make_ops(type::PM, system::System, a::ProdOp{IndexOp})
    r = [ make_ops(type, system, i) for i in a.subs ]
    r[1] *= a.coef
    return r    
end

"""
    apply(op, ::State; maxdim, cutoff)
    apply(op, ::Simulation; maxdim, cutoff)

Apply the given gates to the state and truncate the result according to maxdim and cutoff.
It is much more efficient to apply all the gates in a single call to apply.

# Examples
    apply(CZ(1,3)*H(2)*CNOT(3,4), state; maxdim=10, cutoff=1e-10)

"""
function apply(a::ExprIndexed, state::State; kwargs...)
    ops = make_ops(state.type, state.system, gate(a))
    st = apply(ops, state.state; move_sites_back_between_gates=false, kwargs...)
    return State(state, st)
end
