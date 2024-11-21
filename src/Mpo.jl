struct PreMPO
    sites
    linkdims::Vector{Int}
    terms::Vector{Vector{Tuple{Int, Int, ITensor, Int}}}
    function PreMPO(sites)
        n = length(sites)
        return new(sites, fill(1, n - 1), [ Tuple{Int, Int, ITensor, Int}[] for _ in 1:n ])
    end
end

function create_one_mpo_term(tp, p::ProdLit, pre::PreMPO, ref::Int)
    sites = pre.sites
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
                    push!(pre.terms[i], (1, pre.linkdims[i], make_operator(tp, p.coef * t, sites[i], ref)))
                else
                    push!(pre.terms[i], (pre.linkdims[i-1], pre.linkdims[i], make_operator(tp, t, sites[i]), ref))
                end
                for k in i+1:j-1
                    push!(pre.terms[k],(pre.linkdims[k-1], pre.linkdims[k], delta(sites[k], sites[k]'), ref))
                end
            end
            idx = sites[j]
            if tp ≠ Pure
                idx = Index(dim(idx)÷2, tags(idx))
            end
            t = op(l.opname, idx; l.param...)
            i = j
        else
            idx = noprime(ind(t, 1))
            t = replaceprime(t' * op(l.opname, idx; l.param...), 2 => 1)
        end
    end
    if i == fst
        push!(pre.terms[i], (1, 1, make_operator(tp, p.coef * t, sites[i]), ref))
    else
        push!(pre.terms[i], (pre.linkdims[i-1], 1, make_operator(tp, t, sites[i]), ref))
    end
end

function prepare_mpo(tp, sites, a::SumLit, pre::PreMPO=PreMPO(sites), ref::Int=1)
    for p in a.ps
        if length(p.ls) == 1 && p.ls[1].dissipator
            create_one_mpo_term(tp, p, pre, ref)
        else
            create_one_mpo_term(tp, p, pre, ref)
            if tp == MixEvolve
                create_one_mpo_term(MixEvolve2, -p, pre, ref)
            end
        end
    end
    return pre
end

function prepare_mpo(tp, sites, as::Vector{SumLit})
    pre = PreMPO(sites)
    for (i, a) in enumerate(as)
        prepare_mpo(tp, a, pre, i)
    end
    return pre
end

function make_mpo(p::PreMPO, coefs=[1.0])
    sites = p.sites
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
        for (l, r, u, ref) in p.terms[i]
            if r == 1
                r += rdim
                c = coefs[ref]
                if c ≠ 1
                    u *= c
                end
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

function make_mpo(tp, sites, a::SumLit)
    make_mpo(prepare_mpo(tp, sites, a))
end

function make_approx_W1(p::PreMPO, tau::Number, coefs=[1.0])
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
        for (l, r, u, ref) in p.terms
            if r == 1
                c = tau * coefs[ref]
                u *= tau
                if l == 1
                    u += id
                end
            end
            c = coefs[ref]
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

function make_approx_W2(p::PreMPO, tau::Number, coefs=[1.0])
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
        for (l, r, u, ref) in p.terms
            if r == 1
                c = tau * coefs[ref]
                u *= c
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
