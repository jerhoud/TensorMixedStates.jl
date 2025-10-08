export PreMPO, make_mpo, make_approx_W1, make_approx_W2

struct PreMPO{R <: PM}
    system::System
    linkdims::Vector{Int}
    terms::Vector{Vector{Tuple{Int, Int, ITensor, Int}}}
    function PreMPO{R}(system::System) where R
        n = length(system)
        return new{R}(system, fill(1, n - 1), [ Tuple{Int, Int, ITensor, Int}[] for _ in 1:n ])
    end
end

function PreMPO!(pre::PreMPO{R}, coef::Number, subs::Vector{<:IndexedOp{R}}, ref::Int=1) where R
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    fst = subs[1].index[1]
    lst = subs[end].index[1]
    if fst == lst
        push!(tm[fst], (1, 1, coef * tensor(sys, subs[1]), ref))
    else
        for k in fst:lst-1
            ld[k] += 1
        end    
        push!(tm[fst], (1, ld[fst], coef * tensor(sys, subs[1]), ref))
        i = fst
        for ind in subs[2:end-1]
            j = ind.index[1]
            for k in i+1:j-1
                kdx = SysIndex{R}(sys, k)
                push!(tm[k],(ld[k-1], ld[k], delta(kdx', kdx), ref))
            end
            push!(tm[j], (ld[j-1], ld[j], tensor(sys, ind), ref))
            i = j
        end
        for k in i+1:lst-1
            kdx = SysIndex{R}(sys, k)
            push!(tm[k],(ld[k-1], ld[k], delta(kdx', kdx), ref))
        end
        push!(tm[lst], (ld[lst-1], 1, tensor(sys, subs[end]), ref))
    end
    return pre
end

PreMPO!(pre::PreMPO{R}, a::IndexedOp{R}, ref::Int = 1) where {R <: PM} =
    PreMPO!(pre, scalarcoef(a), prodsubs(a), ref)

function PreMPO!(pre::PreMPO{R}, s::SumOp{R, Indexed}, ref::Int=1) where {R <: PM}
    for p in s.subs
        PreMPO!(pre, p, ref)
    end
    return pre
end

function PreMPO!(pre::PreMPO, as)
    for (i, a) in enumerate(as)
        PreMPO!(a, pre, i)
    end
    return pre
end

"""
    PreMPO(::State, op)

preprocess an operator (or vector of operators).
The result can be passed wherever an operator that must be turned into an MPO is expected
"""
PreMPO(state::State{R}, a; kwargs...) where R =
    PreMPO!(PreMPO{R}(state.system), removeMulti(simplify(a; kwargs...)))

"""
    make_mpo(::PreMPO[, coefs])
    make_mpo(::State, operator)

build an mpo representing an operator
"""
function make_mpo(pre::PreMPO{R}, coefs=[1.]) where R
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    n = length(sys)
    ts = Vector{ITensor}(undef, n)
    rdim = 1
    rlink = Index(2, "Link, l=0")
    for i in 1:n
        idx = SysIndex{R}(sys, i)
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

make_mpo(state::State, a) = make_mpo(PreMPO(state, a))

"""
    make_approx_W1(::PreMPO, tau[, coefs])
    make_approx_W1(::State, operator, tau)

build MPO representing approximation WI of a given operator and time step
"""
function make_approx_W1(pre::PreMPO{R}, tau::Number, coefs=[1.]) where R
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    n = length(sys)
    ts = Vector{ITensor}(undef, n)
    rdim = 1
    rlink = Index(1, "Link, l=0")
    for i in 1:n
        idx = SysIndex{R}(sys, i)
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

make_approx_W1(state::State, a, tau::Number) = make_approx_W1(PreMPO(state, a), tau)

"""
    make_approx_W2(::PreMPO, tau[, coefs])
    make_approx_W2(::State, operator, tau)

build MPO representing approximation WII of a given operator and time step
"""
function make_approx_W2(pre::PreMPO{R}, tau::Number, coefs=[1.]) where R
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    n = length(sys)
    ts = Vector{ITensor}(undef, n)
    rdim = 1
    rlink = Index(1, "Link, l=0")
    for i in 1:n 
        idx = SysIndex{R}(sys, i)
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

make_approx_W2(state::State, a, tau::Number) = make_approx_W2(PreMPO(state, a), tau)
    