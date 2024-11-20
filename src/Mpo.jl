struct PreMPO
    linkdims::Vector{Int}
    terms::Vector{Vector{Tuple{Int, Int, ITensor}}}
end

function create_one_mpo_term(tp, p::ProdLit, sites, pre::PreMPO)
    fst = p.ls[1].index[1]
    lst = p.ls[end].index[1]
    for k in fst:lst-1
        pre.linkdims[k] += 1
    end
    i = 0
    t = ITensor()
    for l in p.ls
        j = l.index[1]
        if l.dissipator && tp ≠ MixDissipator
            die("Dissipators cannot be used in pure representation, as gates or multiplied by other operators")
        end
        if i ≠ j
            if i ≠ 0
                if i == fst
                    push!(pre.terms[i], (1, pre.linkdims[i], make_operator(tp, p.coef * t, sites[i])))
                else
                    push!(pre.terms[i], (pre.linkdims[i-1], pre.linkdims[i], make_operator(tp, t, sites[i])))
                end
                for k in i+1:j-1
                    push!(pre.terms[k],(pre.linkdims[k-1], pre.linkdims[k], delta(sites[k], sites[k]')))
                end
            end
            idx = sites[j]
            if tp ≠ Pure
                idx = Index(dim(idx)÷2, tags(idx))
            end
            t = op(l.name, idx; l.param...)
            i = j
        else
            idx = noprime(ind(t, 1))
            t = replaceprime(t' * op(l.name, idx; l.param...), 2 => 1)
        end
    end
    if i == fst
        push!(pre.terms[i], (1, 1, make_operator(tp, p.coef * t, sites[i])))
    else
        push!(pre.terms[i], (pre.linkdims[i-1], 1, make_operator(tp, t, sites[i])))
    end
end

function prepare_mpo(tp, a::SumLit, sites)
    n = length(sites)
    pre = PreMPO(fill(1, n - 1), [ Tuple{Int, Int, ITensor}[] for _ in 1:n ])
    for p in a.ps
        if length(p.ls) == 1 && p.ls[1].dissipator
            create_one_mpo_term(MixDissipator, p, sites, pre)
        else
            create_one_mpo_term(tp, p, sites, pre)
            if tp == MixEvolve
                create_one_mpo_term(MixEvolve2, -p, sites, pre)
            end
        end
    end
    return pre
end


function make_mpo(tp, a::SumLit, sites)
    p = prepare_mpo(tp, a, sites)
    n = length(sites)
    ts = Vector{ITensor}(undef, n)
    rdim = 1
    rlink = Index(2, "Link, l=0")
    for (i, idx) in enumerate(sites)
        ldim = rdim
        if i == n
            rdim = 1
        else
            rdim = p.linkdims[i]
        end
        llink = rlink
        rlink = Index(1 + rdim, "Link, l=$i")
        w = ITensor(ComplexF64, idx, idx', llink, rlink)
        id = delta(idx, idx')
        for j in eachindval(idx, idx')
            w[llink => 1, rlink => 1, j...] = id[j...]
            w[llink => 1 + ldim, rlink => 1 + rdim, j...] = id[j...]
        end
        for (l, r, u) in p.terms[i]
            if r == 1
                r += rdim
            end
            for j in eachindval(idx, idx')
                w[llink=>l, rlink=>r, j...] += u[j...]
            end
        end
        if i == 1
            w *= ITensor([1, 0], llink)
        end
        if i == n
            w *= ITensor([0, 1], rlink)
        end
        ts[i] = w
    end
    return MPO(ts)
end

function make_approx_W1(tp, a::SumLit, tau::Number, sites)
    p = prepare_mpo(tp, a, sites)
    n = length(sites)
    ts = Vector{ITensor}(undef, n)
    rdim = 1
    rlink = Index(1, "Link, l=0")
    for (i, idx) in enumerate(sites)
        ldim = rdim
        if i == n
            rdim = 1
        else
            rdim = p.dimlinks[i]
        end
        llink = rlink
        rlink = Index(rdim, "Link, l=$i")
        w = ITensor(llink, rlink, idx, idx')
        id = delta(idx, idx')
        for (l, r, u) in p.terms
            if r == 1
                u *= tau
                if l == 1
                    u += id
                end
            end
            for j in eachindval(idx, idx')
                w[llink=>l, rlink=>r, j...] += u[j...]
            end
        end
        if i == 1
            w *= ITensor([1], llink)
        end
        if i == n
            w *= ITensor([1], rlink)
        end
        ts[i] = w
    end
    return MPO(ts)
end

function make_approx_W2(tp, a::SumLit, tau::Number, sites)
    p = prepare_mpo(tp, a, sites)
    n = length(sites)
    ts = Vector{ITensor}(undef, n)
    rdim = 1
    rlink = Index(1, "Link, l=0")
    for (i, idx) in enumerate(sites)
        ldim = rdim
        if i == n
            rdim = 1
        else
            rdim = p.dimlinks[i]
        end
        llink = rlink
        rlink = Index(rdim, "Link, l=$i")
        dl = dim(llink)
        dr = dim(rlink)
        v = fill(ITensor(), (dl, dr))
        for (l, r, u) in p.terms
            if r == 1
                u *= tau
            end
            v[l, r] += u
        end
        e = v[1, 1] = exp(v[1, 1])
        for l in 2:dl, r in 2:dr
            v[l, r] += replaceprime(v[l, 1]'' * e' * v[1, r], 3=>1)
        end
        for r in 2:dr
            v[1, r] = replaceprime(e' * v[1, r], 2=>1)
        end
        for l in 2:dl
            v[l, 1] = replaceprime(v[l, 1]' * e, 2=>1)
        end
        
        w = ITensor(llink, rlink, idx, idx')
        for l in 1:dl, r in 1:dr
            u = v[l, r]
            if isempty(u) continue end
            for j in eachindval(idx, idx')
                w[llink=>l, rlink=>r, j...] += u[j...]
            end
        end
        if i == 1
            w *= ITensor([1], llink)
        end
        if i == n
            w *= ITensor([1], rlink)
        end
        ts[i] = w
    end
    return MPO(ts)
end
