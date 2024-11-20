function Pure2Mixed(st::MPS; maxdim, cutoff)
    s = dense(st) # to remove QNs
    n = length(s)
    v = Vector{ITensor}(undef, n)
    leftlink = Index(0)
    r = ITensor()
    for (i, t) in enumerate(s)
        idx = siteind(s, i)
        midx = mixed_index(idx)
        mt = t * t' * combinerto(idx', idx, midx)
        if (i > 1)
            leftidx = (midx, leftlink)
            mt *= r
        else
            leftidx = (siteidx,)
        end
        if (i < n)
            ut, st, vt = svd(mt, leftidx...; lefttags="Link,l=$i", maxdim, cutoff)
            leftlink = commonind(ut, st)
            mt = ut
            r = st * vt
        end
        v[i] = mt
    end
    return MPS(v)
end

function random_mixed_state(sites, linkdims::Int)
    n = length(sites)
    super = vcat(sites, sites)
    smps = random_mps(ComplexF64, super; linkdims = (linkdims + 1) ÷ 2)
    smps = Pure2Mixed(smps; maxdim = linkdims, cutoff = 0)
    t = ITensor(1)
    for i in 2n:-1:n+1
        t *= smps[i] * obs(siteind(smps, i))
    end
    smps[n] *= t
    return MPS(smps[1:n])
end

trace(p::MPS) =
    scalar(prod(pi * obs(siteind(p, i)) for (i, pi) in enumerate(p)))

trace2(p::MPS) = real(inner(p, p))

struct Preprocess
    loc::Vector{ITensor}
    left::Vector{ITensor}
    right::Vector{ITensor}
end

preprocess(::Type{Pure}, ::MPS) = nothing

function preprocess(::Type{Mixed}, st::MPS)
    n = length(st)
    vloc = [st[i] * obs(siteind(st, i)) for i in 1:n]
    v = ITensor(1)
    vleft = vcat([st[1]], [(v *= vloc[i]; v * st[i+1]) for i in 1:n-1])
    v = ITensor(1)
    vright = reverse(vcat([ITensor(1)], [v *= vloc[i] for i in n:-1:2]))
    return Preprocess(vloc, vleft, vright)
end

unroll(x) =
    if x[1] isa Number
        x
    else
        v = [map(t->t[i], x) for i in 1:length(x[1])]
        if x[1] isa Matrix
            return reshape(v, size(x[1]))
        else
            return v
        end
    end


function expect(::Type{Pure}, st::MPS, p::ProdLit, ::Nothing)
    if p.coef == 0
        return 0
    end
    imin = p.ls[1].index[1]
    imax = p.ls[end].index[1]
    if isortho(st)
        c = orthocenter(st)
        if c < imin
            orthogonalize!(st, imin)
        elseif c > imax
            orthogonalize!(st, imax)
        end
    else
        orthogonalize!(st, imin)
    end
    r = ITensor(1)
    o = ITensor()
    j = 0
    foreach(p.ls) do l
        i = l.index[1]
        if j == i
            idx = noprime(ind(o, 1))
            o = replaceprime(o' * op(l.opname, idx; l.param...))
        else
            if j == 0
                if i ≠ 1
                    idx = commonind(st[i-1], st[i])
                    r = delta(idx, idx')
                end
            else
                r *= st[j]
                r *= o
                r *= dag(st[j]')
                for k in j+1:i-1
                    r *= st[k]
                    r *= dag(st[k]')
                end
            end
            j = i
            o = op(l.opname, siteind(st, i); l.param...)
        end
    end
    r *= st[j]
    if j ≠ length(st)
        idx = commonind(st[j], st[j+1])
        r *= delta(idx, idx')
    end
    r *= o
    r *= dag(st[j])
    return p.coef * scalar(r)
end

function expect(::Type{Mixed}, st::MPS, p::ProdLit, prep::Preprocess)
    r = ITensor(1)
    o = ITensor()
    j = 0
    foreach(p.ls) do l
        i = l.index[1]
        if j == i
            idx = noprime(ind(o, 1))
            o = replaceprime(o' * op(l.opname, idx; l.param...))
        else
            if j == 0
                r = prep.left[i]
            else
                r *= make_operator(MixObservable, o, siteind(st, j))
                for k in j+1:i-1
                    r *= prep.loc[k]
                end
                r *= st[i]
            end
            j = i
            idx = pure_index(siteind(st, i))
            o = op(l.opname, idx; l.param...)
        end
    end
    if j ≠ 0
        r *= make_operator(MixObservable, o, siteind(st, j))
        r *= prep.right[j]
    end
    return p.coef * scalar(r)
end

expect(tp, s, op::SumLit, prep) =
    sum(op.ps; init=0) do p
        expect(tp, s, p, prep)
    end

expect(tp, s, op, prep) =
    map(op) do o
        expect(tp, s, o, prep)
    end

expect1_one(::Type{Pure}, o::String, i::Index, t::ITensor) =
    scalar(op(o, i) * t)

expect1_one(::Type{Mixed}, o::String, i::Index, t::ITensor) =
    scalar(t * make_operator(MixObservable, op(o, pure_index(i)), i))
    
expect1_one(tp, op, i::Index, t::ITensor) =
    map(op) do o
        expect1_one(tp, o, i, t)
    end


function expect1(::Type{Pure}, st::MPS, op, ::Nothing)
    n = length(st)
    r = Vector(undef, n)
    for i in 1:n
        orthogonalize!(st, i)
        t = st[i]
        if i > 1
            idx = commonind(st[i-1], st[i])
            t *= delta(idx, idx')
        end
        if i < n
            idx = commonind(st[i], st[i+1])
            t *= delta(idx, idx')
        end
        t *= dag(st[i]')
        r[i] = expect1_one(Pure, op, siteind(st, i), t)
    end
    return unroll(r)
end
    
function expect1(::Type{Mixed}, st::MPS, op, prep::Preprocess)
    n = length(st)
    r = [expect1_one(Mixed, op, siteind(st, i), prep.left[i] * prep.right[i]) for i in 1:n]
    return unroll(r)
end

function expect2_one(::Type{Pure}, ops::Tuple{String, String}, i1::Index, i2::Index, t::ITensor, rev::Bool)
    if rev
        o2, o1 = ops
    else
        o1, o2 = ops
    end
    if i1 == i2
        o = replaceprime(op(o1, i1') * op(o2, i1), 2=>1)
        return scalar(o * t)
    else
        o1 = op(o1, idx1)
        o2 = op(o2, idx2)
        return scalar(o1 * t * o2)
    end
end

function expect2_one(::Type{Mixed}, ops::Tuple{String, String}, i1::Index, i2::Index, t::ITensor, rev::Bool)
    if rev
        o2, o1 = ops
    else
        o1, o2 = ops
    end
    if i1 == i2
        idx = pure_index(i1)
        o = make_operator(MixObservable, replaceprime(op(o1, idx') * op(o2, idx), 2=>1), i1) 
        return scalar(o * t)
    else
        idx1 = pure_index(i1)
        o1 = make_operator(MixObservable, op(o1, idx1), i1)
        idx2 = pure_index(i2)
        o2 = make_operator(MixObservable, op(o2, idx2), i2)
        return scalar(o1 * t * o2)
    end
end

expect2_one(tp, ops, i1::Index, i2::Index, t::ITensor, rev::Bool) =
    map(ops) do op
        expect2_one(tp, op, i1, i2, t, rev)
    end

function expect2(::Type{Pure}, st::MPS, ops, prep::Nothing)
    n = length(st)
    sites = siteinds(st)
    c = Matrix(undef, n, n)
    for (i, idx) in enumerate(sites)
        orthogonalize!(st, i)
        l = st[i]
        if i > 1
            lidx = commonind(st[i-1], st[i])
            l *= delta(lidx, lidx')
        end
        l *= dag(st[i]')
        t = l
        if i < n
            ridx = commonind(st[i], st[i+1])
            t *= delta(ridx, ridx')
        end
        c[i, i] = expect2_one(Pure, ops, idx, idx, t, false)
        for j in i+1:n
            l *= st[j]
            t = l
            if i < n
                ridx = commonind(st[i], st[i+1])
                t *= delta(ridx, ridx')
            end
            t *= dag(st[j]')
            jdx = sites[j]
            c[i, j] = expect2_one(Mixed, ops, idx, jdx, t, false)
            c[j, i] = expect2_one(Mixed, ops, idx, jdx, t, true)
            l *= delta(jdx, jdx') * dag(st[j]')
        end
    end
    return unroll(c)
end
    
function expect2(::Type{Mixed}, st::MPS, ops, prep::Preprocess)
    n = length(st)
    sites = siteinds(st)
    r = Matrix(undef, n, n)
    for (i, idx) in enumerate(sites)
        l = prep.left[i]
        r[i, i] = expect2_one(Mixed, ops, idx, idx, l * prep.right[i], false)
        for j in i+1:n
            t = l * (st[j] * prep.right[j])
            r[i, j] = expect2_one(Mixed, ops, idx, sites[j], t, false)
            r[j, i] = expect2_one(Mixed, ops, idx, sites[j], t, true)
            l *= prep.loc[j]
        end
    end
    return unroll(r)
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