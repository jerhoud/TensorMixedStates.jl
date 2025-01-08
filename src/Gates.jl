export apply_gates

function Lit_to_ops(state::State, a::ProdLit)
    if a.coef == 0
        return []
    end
    s = state.system
    r = map(a.ls) do l
        if l.dissipator
            error("Cannot apply a dissipator as a gate")
        end
        idx = map(i->s.pure_sites[i], l.index)
        oi = op(l.opname, idx...; l.param...)
        if state.type == Pure
            return oi
        else
            jdx = sim(idx)
            oj = replaceinds(oi, (idx..., idx'...), (jdx..., jdx'...))
            kdx = map(i->s.mixed_sites[i], l.index)
            return *(oi, dag(oj), combinerto.(jdx, idx, kdx)..., combinerto.(jdx', idx', kdx')...)
        end
    end
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
function apply(a::ProdLit, state::State; kwargs...)
    ops = Lit_to_ops(state, insertFfactors(a))
    st = apply(ops, state.state; move_sites_back_between_gates=false, kwargs...)
    return State(state, st)
end
