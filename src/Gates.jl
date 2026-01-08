export apply

"""
    apply(op, ::State; limits::Limits)
    apply(mps, ::State; limits::Limits)
    apply(op, ::Simulation; limits::Limits)

Apply the given gates to the state and truncate the result according to limits.
It is much more efficient to apply all the gates in a single call to apply.

# Examples
    apply(CZ(1,3)*H(2)*CNOT(3,4), state)

"""
apply(a::IndexedOp{Pure}, state::State{Mixed}; kwargs...) =
    apply(Gate(a), state; kwargs...)

function apply(a::IndexedOp{R}, state::State{R}; limits::Limits=Limits()) where R
    ops = make_ops(state.system, a)
    st = apply(ops, state.state; move_sites_back_between_gates=false,
            limits.cutoff, limits.maxdim)
    return State(state, st)
end

apply(mpo::MPO, state::State; limits::Limits=Limits()) =
    State(state, apply(mpo, state.state; limits.cutoff, limist.maxdim))
    
make_ops(::System, a::SumOp) =
    error("cannot apply sums as gates ($a)")

make_ops(s::System, a::ScalarOp) =
    if a.coef == 0
        error("cannot apply null gate")
    else
        ops = make_ops(s, a.arg)
        if !isempty(ops)
            ops[1] *= a.coef
        end
        ops 
    end    

make_ops(s::System, a::ProdOp) =
    reduce(vcat, map(x->make_ops(s, x), a.subs))

make_ops(s::System, a::AtIndex{R, N}) where {R, N} =
    if a == MakeIdentity{R, Indexed, N}
        []
    else
        [ tensor(s, a) ]
    end

function make_ops(s::System, a::Multi_F)
    ops = []
    for i in a.start:a.stop
        f = F_info(s[i])
        if f == Id
            continue
        end
        push!(ops, tensor(f, s[i]))
    end
    ops
end    
