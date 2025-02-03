export apply_gates

function Lit_to_ops(::TPure, system::System, a::ProdLit)
    r = [ l(system.pure_sites) for l in a.ls ]
    r[1] *= a.coef
    return r    
end

function Lit_to_ops(::TMixed, system::System, a::ProdLit)
    r = map(a.ls) do l
        try
            l(system.mixed_sites)
        catch
            make_operator(MixGate, system, l(system.pure_sites), l.index)
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
    if a.coef == 0
        error("Cannot apply a null gate")
    end
    if dissipLit(a)
        error("Cannot apply a dissipator as a gate")
    end
    ops = Lit_to_ops(state.type, state.system, insertFfactors(a))
    st = apply(ops, state.state; move_sites_back_between_gates=false, kwargs...)
    return State(state, st)
end
