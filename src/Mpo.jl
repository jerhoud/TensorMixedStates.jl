export PreMPO, make_mpo, make_approx_W1, make_approx_W2

struct MixEvolve end
struct MixEvolve2 end
struct MixDissipator end

function make_operator(::TPure, ::System, t::ITensor, ::Int)
    return t
end
    
function make_operator(::Type{MixEvolve}, system::System, t::ITensor, i::Int)
    idx = system.pure_index[i]
    jdx = sim(idx)
    kdx = system.mixed_sites[i]
    return *(t, delta(jdx, jdx'), combinerto(jdx, idx, kdx), combinerto(jdx', idx', kdx'))
end

function make_operator(::Type{MixEvolve2}, system::System, t::ITensor, i::Int)
    idx = system.pure_index[i]
    jdx = sim(idx)
    kdx = system.mixed_sites[i]
    return *(dag(t), delta(jdx, jdx'), combinerto(idx, jdx, kdx), combinerto(idx', jdx', kdx'))
end

function make_operator(::Type{MixDissipator}, system::System, ti::ITensor, i::Int)
    idx = system.pure_sites[i]
    jdx = sim(idx)
    kdx = system.mixed_sites[i]
    tj = replaceinds(ti, (idx, idx'), (jdx, jdx'))
    ati = swapprime(dag(ti'), 1=>2)
    atj = swapprime(dag(tj'), 1=>2)
    r = ti * dag(tj) -
        0.5 * *(replaceprime(ati * ti, 2 => 1), delta.(jdx, jdx')) -
        0.5 * *(replaceprime(atj * tj, 2 => 1), delta.(idx, idx'))
    return *(r, combinerto(jdx, idx, kdx), combinerto(jdx', idx', kdx'))
end

struct PreMPO
    type::Union{TPure, TMixed}
    system::System
    linkdims::Vector{Int}
    terms::Vector{Vector{Tuple{Int, Int, ITensor, Int}}}
    function PreMPO(state::State)
        n = length(state)
        return new(state.type, state.system, fill(1, n - 1), [ Tuple{Int, Int, ITensor, Int}[] for _ in 1:n ])
    end
end

function PreMPO(tp, pre::PreMPO, p::ProdLit, ref::Int)
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
        if l.dissipator && tp ≠ MixDissipator
            error("dissipators cannot be used in pure representations or multiplied by other operators")
        end
        j = l.index[1]
        o = op(l.opname, sys.pure_sites[j]; l.param...)
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

PreMPO(::TPure, p::ProdLit, pre::PreMPO, ref::Int=1) =
    PreMPO(Pure, pre, p, ref)

function PreMPO(::TMixed, p::ProdLit, pre::PreMPO, ref::Int=1)
    if length(p.ls) == 1 && p.ls[1].dissipator
        PreMPO(MixDissipator, pre, p, ref)
    else
        PreMPO(MixEvolve, pre, p, ref)
        PreMPO(MixEvolve2, pre, -p, ref)
    end
    return pre
end

function PreMPO(type::Union{TPure, TMixed}, a::SumLit, pre::PreMPO, ref::Int=1)
    for p in a.ps
        PreMPO(type, p, pre, ref)
    end
    return pre
end

function PreMPO(type::Union{TPure, TMixed}, as, pre::PreMPO)
    for (i, a) in enumerate(as)
        PreMPO(type, a, pre, i)
    end
    return pre
end

PreMPO(state::State, a) =
    PreMPO(state.type, insertFfactors(a), PreMPO(state))


function make_mpo(pre::PreMPO, coefs=(1.0,))
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

make_mpo(type::Union{TPure, TMixed}, system::System, a::SumLit) = 
    make_mpo(PreMPO(type, system, a))

function make_approx_W1(pre::PreMPO, tau::Number, coefs=(1.0,))
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
            if r == 1
                c = tau * coefs[ref]
                u *= c
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

function make_approx_W2(pre::PreMPO, tau::Number, coefs=(1.0,))
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
            if r == 1
                c = tau * coefs[ref]
                u *= c
            end
            v[l, r] += u
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
