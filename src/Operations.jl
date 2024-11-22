
struct MixObservable end
struct MixEvolve end
struct MixEvolve2 end
struct MixGate end
struct MixDissipator end

function make_operator(::Type{MixEvolve}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    return *(tj, delta(k, k'), combinerto(k, j, i), combinerto(k', j', i ))
end

function make_operator(::Type{MixEvolve2}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    return *(tk, delta(j, j'), combinerto(k, j, i), combinerto(k', j', i ))
end

function make_operator(::Type{MixGate}, tj::ITensor, i::Index...)
    j = filter(t->hasplev(t, 0), inds(tj))
    k = sim(j)
    tk = replaceinds(tj, (j..., j'...), (k..., k'...))
    return *(tj, dag(tk), combinerto.(k, j, i)..., combinerto.(k', j', i')...)
end

function make_operator(::Type{MixDissipator}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    atj = swapprime(dag(tj'), 1=>2)
    atk = swapprime(dag(tk'), 1=>2)
    r = tj * dag(tk) -
        0.5 * *(replaceprime(atj * tj, 2 => 1), delta.(k, k')) -
        0.5 * *(replaceprime(atk * tk, 2 => 1), delta.(j, j'))
    return *(r, combinerto(k, j, i), combinerto(k', j', i'))
end




function Lit_to_ops(::Type{Pure}, a::ProdLit, sites)
    if a.coef == 0
        return []
    end
    r = [ op(sites, l.opname, l.index...; l.param...) for l in a.ls ]
    r[1] *= a.coef
    return r    
end
    
function Lit_to_ops(::Type{Mixed}, a::ProdLit, sites)
    if a.coef == 0
        return []
    end
    r = map(a.ls) do l
        idx = map(i->sites[i], l.index) 
        jdx = pure_index.(idx)
        o = make_operator(MixGate, op(l.opname, jdx...; l.param), idx...)
    end
    r[1] *= a.coef
    return r
end

function apply_gate(tp, p::MPS, g::ProdLit; kwargs...)
    ops = Lit_to_ops(tp, g, siteinds(p))
    return apply(ops, p; move_sites_back_between_gates=false, kwargs...)
end

function entanglement_entropy(s::MPS, pos::Int)
    s = orthogonalize(s, pos)
    _,S = svd(s[pos], (linkinds(s, pos-1)..., siteinds(s, pos)...))
    sp = [ S[i,i]^2 for i in 1:dim(S, 1) ]
    ee = -sum(p * log(p) for p in sp)
    return (ee, sp)
end

function signature(p)
    n = length(p)
    t = fill(false, n)
    r = 1
    for i in 1:n
        if t[i] continue end
        j = p[i]
        while j ≠ i
            r = -r
            t[j] = true
            j = p[j]
        end
    end
    return r
end

function reorder(a::ProdLit)
    fs = filter(t->t.fermionic, a.ls)
    p = sortperm(fs; by=l->l.index)
    s = signature(p)
    ProdLit(s * a.coef, sort(a.ls; by=l->l.index))
end

function reorder(a::SumLit)
    ps = sort(map(reorder, a.ps))
    pps = ProdLit[]
    c = 0
    l = Lit[]
    for p in ps
        if p.ls == l
            c += p.coef
        else
            if c ≠ 0
                push!(pps, ProdLit(c, l))
            end
            c = p.coef
            l = p.ls
        end
    end
    if c ≠ 0
        push!(pps, ProdLit(c, l))
    end
    return SumLit(pps)
end

litF(idx) = Lit("F", "F", (idx,), (;), false, false)

function insertFfactors(a::ProdLit)
    nls = Lit[]
    fermion_idx = 0
    for l in reverse(a.ls)
        idx = first(l.index)
        if fermion_idx ≠ 0
            append!(nls, (litF(i) for i in reverse(idx + 1:fermion_idx - 1)))
            if l.fermionic
                push!(nls, litF(idx))
                fermion_idx = 0
            else
                fermion_idx = idx
            end
        elseif l.fermionic
            fermion_idx = idx
        end
        push!(nls, l)
    end
    if fermion_idx ≠ 0
        append!(nls, (litF(i) for i in reverse(1:fermion_idx - 1)))
    end
    return ProdLit(a.coef, reverse(nls))
end

insertFfactors(a::SumLit) =
    SumLit(map(insertFfactors, a.ps))