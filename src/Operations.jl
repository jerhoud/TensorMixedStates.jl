mpos::Dict{Lits, MPO} = IdDict()

clear_mpos() =
    global mpos = IdDict()

function get_mpo(op::Lits, st::MPS)
    get!(mpos, op) do
        MPO(Lit_to_OpSum(op), siteinds(st))
    end
end

function Pure2Mixed(st::MPS; maxdim, cutoff)
    p = dense(st) # to remove QNs
    n = length(p)
    v = Vector{ITensor}(undef, n)
    leftlink = nothing
    r = nothing
    for (i, t) in enumerate(p)
        idx = siteind(p, i)
        midx = mixed_index(idx)
        mt = t * t' * combinerto(idx', idx, midx)
        if (i > 1)
            leftidx = (midx, leftlink)
            mt *= r
        else
            leftidx = (siteidx,)
        end
        if (i < n)
            ut, st, vt = svd(tt, leftidx...; lefttags="Link,l=$i", maxdim, cutoff)
            leftlink = commonind(ut, st)
            tt = ut
            r = st * vt
        end
        v[i] = tt
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
        t *= smps[i] * make_operator(MixObservable, siteind(smps, i))
    end
    smps[n] *= t
    return MPS(smps[1:n])
end

trace(p::MPS) =
    prod(pi * make_opeartor(MixObservable, siteind(p, i)) for (i, pi) in enumerate(p))[1]

trace2(p::MPS) = real(inner(p, p))

struct Preprocess
    loc::Vector{ITensor}
    left::Vector{ITensor}
    right::Vector{ITensor}
end

preprocess(::Type{Pure}, ::MPS) = nothing

function preprocess(::Type{Mixed}, p::MPS)
    n = length(p)
    vloc = [p[i] * op("obs", siteind(p, i)) for i in 1:n]
    v = ITensor(1)
    vleft = vcat([p[1]], [(v *= vloc[i]; v * p[i+1]) for i in 1:n-1])
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


expect(::Type{Pure}, p::MPS, op::Union{String, Vector{String}}, prep::Nothing) =
    ITensorMPS.expect(p, op)

expect(::Type{Pure}, p::MPS, op::ProdLit, prep::Nothing) =
    inner(p', get_mpo(op, p), p)

function expect(::Type{Mixed}, p::MPS, pl::ProdLit, prep::Preprocess)
    sites = siteinds(p)
    n = length(p)
    t::Vector{Any} = fill(nothing, n)
    foreach(pl.ls) do l
        if length(l.index) ≠ 1
            die("Multiple site observables are not supported")
        else
            i = l.index[1]
            obs = op(sites, l.name, i)
            if isnothing(t[i])
                t[i] = obs
            else
                t[i] = replaceprime(prime(t[i]) * obs, 2 => 1)
            end
        end
    end
    r = ITensor(1)
    j = 0
    for i in 1:n
        o = t[i]
        if isnothing(o)
            continue
        end
        o *= op("obs", siteind(p, i)')
        if j == 0
            r = prep.left[i] * o
        else
            for k in j+1:i-1
                r *= prep.loc[k]
            end
            r *= p[i] * o
        end
        j = i
    end
    r *= prep.right[j]
    return r[1] * pl.coef
end

expect(tp, s, op::Union{Vector{ProdLit}}, prep) =
    map(op) do o 
        expect(tp, s, o, prep)
    end



expect_one(o::String, i::Index, t::ITensor) =
    (op("obs", i') * op(o, i) * t)[1]
    
expect_one(op::Union{Vector{String}}, i::Index, t::ITensor) =
    map(op) do o
        expect_one(o, i, t)
    end
    
function expect(::Type{Mixed}, p::MPS, op::Union{String, Vector{String}}, prep::Preprocess)
    n = length(p)
    r = [expect_one(op, siteind(p, i), prep.left[i] * prep.right[i]) for i in 1:n]
    return unroll(r)
end




correlations(::Type{Pure}, p::MPS, ops::Union{Vector{String}, Tuple{Vararg{String}}}, prep::Nothing) =
    correlation_matrix(p, ops[1], ops[2])

correlations(::Type{Pure}, p::MPS, ops, prep::Nothing) =
    map(ops) do op
        correlations(Pure, p, op, nothing)
    end


function correlations_one(ops::Union{Vector{String}, Tuple{Vararg{String}}}, i1::Index, i2::Index, t::ITensor)
    if i1 == i2
        obs = op("obs", i1'') * op(ops[1], i1') * op(ops[2], i1)
        return (obs * t)[1]
    else
        obs1 = op("obs", i1') * op(ops[1], i1)
        obs2 = op("obs", i2') * op(ops[2], i2)
        return (obs1 * t * obs2)[1]
    end
end

correlations_one(ops, i1::Index, i2::Index, t::ITensor) =
    map(ops) do op
        correlations_one(op, i1, i2, t)
    end


function correlations(::Type{Mixed}, p::MPS, ops, prep::Preprocess)
    n = length(p)
    sites = siteinds(p)
    r = Matrix{Any}(undef, n, n)
    for i in 1:n
        idx = sites[i]
        l = prep.left[i]
        r[i, i] = correlations_one(ops, idx, idx, l * prep.right[i])
        for j in i+1:n
            r[i, j] = r[j, i] = correlations_one(ops, idx, sites[j], l * p[j] * prep.right[j])
            l *= prep.loc[j]
        end
    end
    return unroll(r)
end

function apply_gate(p::MPS, g::ProdLit; kwargs...)
    ops = Lit_to_ops(g, siteinds(p))
    return apply(ops, p; move_sites_back_between_gates=false, kwargs...)
end

function entanglement_entropy(s::MPS, pos::Int)
    s = orthogonalize(s, pos)
    _,S = svd(s[pos], (linkinds(s, pos-1)..., siteinds(s, pos)...))
    sp = [ S[i,i]^2 for i in 1:dim(S, 1) ]
    ee = -sum(p * log(p) for p in sp)
    return (ee, sp)
end