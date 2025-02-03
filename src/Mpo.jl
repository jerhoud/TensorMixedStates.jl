export PreMPO, make_mpo, make_approx_W1, make_approx_W2

"""
    struct PreMPO
    PreMPO(state, operator)

a data type to hold precomputations for building an MPO
"""
struct PreMPO
    type
    system::System
    linkdims::Vector{Int}
    terms::Vector{Vector{Tuple{Int, Int, ITensor, Int}}}
    function PreMPO(type, state::State)
        if state.type == Pure
            type = Pure
        end
        n = length(state)
        return new(type, state.system, fill(1, n - 1), [ Tuple{Int, Int, ITensor, Int}[] for _ in 1:n ])
    end
end

function PreMPO!(tp, pre::PreMPO, p::ProdLit, ref::Int)
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    fst = p.ls[1].index[1]
    lst = p.ls[end].index[1]
    for k in fst:lst-1
        ld[k] += 1
    end
    i = 0
    t = ITensor()
    for l in p.ls
        j = l.index[1]
        o = l(sys.pure_sites)
        if i ≠ j
            if i ≠ 0
                if i == fst
                    push!(tm[i], (1, ld[i], make_operator(tp, sys, p.coef * t, i), ref))
                else
                    push!(tm[i], (ld[i-1], ld[i], make_operator(tp, sys, t, i), ref))
                end
                for k in i+1:j-1
                    kdx = sites(tp, sys)[k]
                    push!(tm[k],(ld[k-1], ld[k], delta(kdx, kdx'), ref))
                end
            end
            t = o
            i = j
        else
            t = replaceprime(t' * o, 2 => 1)
        end
    end
    if i == fst
        push!(tm[i], (1, 1, make_operator(tp, sys, p.coef * t, i), ref))
    else
        push!(tm[i], (ld[i-1], 1, make_operator(tp, sys, t, i), ref))
    end
    return pre
end

function PreMPO!(p::ProdLit, pre::PreMPO, ref::Int=1)
    if p.coef == 0
        return
    end
    if dissipLit(p) && (length(p.ls) ≠ 1 || pre.type ≠ MixEvolve)
        error("dissipators cannot be used on pure representations or multiplied by other operators")
    end
    if pre.type == MixEvolve
        p1 = p.ls[1]
        if p1.op.dissipator
            if p1.op.fermionic
                PreMPO!(MixGate, pre, ProdLit(p.coef, [[ litF(i) for i in 1:p1.index[1]-1 ]; p.ls]), ref)
                PreMPO!(MixDissipatorF, pre, p, ref)
            else
                PreMPO!(MixDissipator, pre, p, ref)
            end
        else
            PreMPO!(MixEvolve, pre, p, ref)
            PreMPO!(MixEvolve2, pre, p, ref)
        end
    else
        PreMPO!(pre.type, pre, p, ref)
    end
    return pre
end

function PreMPO!(a::SumLit, pre::PreMPO, ref::Int=1)
    for p in a.ps
        PreMPO!(p, pre, ref)
    end
    return pre
end

function PreMPO!(as, pre::PreMPO)
    for (i, a) in enumerate(as)
        PreMPO!(a, pre, i)
    end
    return pre
end

function PreMPO(state::State, a, tp=MixEvolve)
    PreMPO!(insertFfactors(a), PreMPO(tp, state))
end

"""
    make_mpo(PreMPO[, coef])
    make_mpo(state, operator)

build an mpo representing an operator
"""
function make_mpo(pre::PreMPO, coefs=[1.])
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    sit = sites(pre.type, sys)
    n = length(sys)
    ts = Vector{ITensor}(undef, n)
    rdim = 1
    rlink = Index(2, "Link, l=0")
    for (i, idx) in enumerate(sit)
        ldim = rdim
        llink = rlink
        if i == n
            rdim = 1
        else
            rdim = ld[i]
        end
        rlink = Index(1 + rdim, "Link, l=$i")
        w = ITensor(ComplexF64, idx', idx, llink, rlink)
        id = delta(idx, idx')
        for j in eachindval(idx, idx')
            w[llink => 1, rlink => 1, j...] = id[j...]
            w[llink => 1 + ldim, rlink => 1 + rdim, j...] = id[j...]
        end
        for (l, r, u, ref) in tm[i]
            c = coefs[ref]
            if c ≠ 0
                if r == 1
                    r += rdim
                end
                for j in eachindval(idx, idx')
                    w[llink=>l, rlink=>r, j...] += c * u[j...]
                end
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

make_mpo(state, a, tp = MixEvolve) = 
    make_mpo(PreMPO(state, a, tp))

"""
    make_approx_W1(PreMPO, tau[, coefs])
    make_approx_W1(state, operator, tau)

build MPO representing approximation WI of a given operator and time step
"""
function make_approx_W1(pre::PreMPO, tau::Number, coefs=[1.])
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    sit = sites(pre.type, sys)
    n = length(sys)
    ts = Vector{ITensor}(undef, n)
    rdim = 1
    rlink = Index(1, "Link, l=0")
    for (i, idx) in enumerate(sit)
        llink = rlink
        if i == n
            rdim = 1
        else
            rdim = ld[i]
        end
        rlink = Index(rdim, "Link, l=$i")
        w = ITensor(ComplexF64, idx', idx, llink, rlink)
        id = delta(idx, idx')
        for j in eachindval(idx, idx')
            w[llink=>1, rlink=>1, j...] = id[j...]
        end
        for (l, r, u, ref) in tm[i]
            c = coefs[ref]
            if c ≠ 0
                if r == 1
                    c *= tau
                end
                for j in eachindval(idx, idx')
                    w[llink=>l, rlink=>r, j...] += c * u[j...]
                end
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

make_approx_W1(state, a, tau::Number) = 
    make_approx_W1(PreMPO(state, a, MixEvolve), tau)

"""
    make_approx_W2(PreMPO, tau[, coefs])
    make_approx_W2(state, operator, tau)

build MPO representing approximation WII of a given operator and time step
"""
function make_approx_W2(pre::PreMPO, tau::Number, coefs=[1.])
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    sit = sites(pre.type, sys)
    n = length(sys)
    ts = Vector{ITensor}(undef, n)
    rdim = 1
    rlink = Index(1, "Link, l=0")
    for (i, idx) in enumerate(sit)
        ldim = rdim
        if i == n
            rdim = 1
        else
            rdim = ld[i]
        end
        llink = rlink
        rlink = Index(rdim, "Link, l=$i")
        v = fill(ITensor(), (ldim, rdim))
        for (l, r, u, ref) in tm[i]
            c = coefs[ref]
            if c ≠ 0
                if r == 1
                    c *= tau
                end
                v[l, r] += c * u
            end
        end
        d = v[1, 1]
        if isempty(d)
            e = delta(idx, idx')
        else
            e = exp(d)
        end
        v[1, 1] = e
        for l in 2:ldim, r in 2:rdim
            vl = v[l, 1]
            vr = v[1, r]
            if !isempty(vl) && !isempty(vr)
                v[l, r] += replaceprime(vl'' * e' * vr, 3=>1)
            end
        end
        for r in 2:rdim
            vr = v[1, r]
            if !isempty(vr)
                v[1, r] = replaceprime(e' * vr, 2=>1)
            end
        end
        for l in 2:ldim
            vl = v[l, 1]
            if !isempty(vl)
                v[l, 1] = replaceprime(vl' * e, 2=>1)
            end
        end
        
        w = ITensor(ComplexF64, idx', idx, llink, rlink)
        for l in 1:ldim, r in 1:rdim
            u = v[l, r]
            if !isempty(u)
                for j in eachindval(idx, idx')
                    w[llink=>l, rlink=>r, j...] += u[j...]
                end
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

make_approx_W2(state, a, tau::Number) = 
    make_approx_W2(PreMPO(state, a, MixEvolve), tau)
    