mpos::Dict{Tuple{String, Lits}, MPO} = IdDict()

clear_mpos() =
    global mpos = IdDict()

function get_mpo(op::Lits, st::MPS, ext::String="")
    get!(mpos, (ext, op)) do
        MPO(Lit_to_OpSum(op, ext), siteinds(st))
    end
end

function Pure2Mixed(st::MPS; maxdim, cutoff)
    p = dense(st) # to remove eventual QNs
    n = length(p)
    v = Vector{ITensor}(undef, n)
    leftlink = nothing
    r = nothing
    for i in 1:n
        t = p[i]
        idx = siteind(p, i)
        tp = "Mixed" * string(tags(idx)[1]) # hoping that is the tag with the type of the site
        combsite = combiner(prime(idx), idx; tags="$tp,Site,n=$i")
        siteidx = combinedind(combsite)
        tt = t * dag(prime(t)) * combsite
        if (i > 1)
            leftidx = (siteidx, leftlink)
            tt *= r
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



trace(::Type{Pure}, ::MPS) = 1

function trace(::Type{Mixed}, p::MPS)
    n = length(p)
    return prod(p[i] * op("obs", siteind(p, i)) for i in 1:n)[1]
end

trace2(::Type{Pure}, ::MPS) = 1

trace2(::Type{Mixed}, p::MPS) = inner(p, p)

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
            obs = op(sites, "obs" * l.name, i)
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
    (op("obs", i') * op("obs" * o, i) * t)[1]
    
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
        correlations(s, p, op, nothing)
    end


function correlations_one(ops::Union{Vector{String}, Tuple{Vararg{String}}}, i1::Index, i2::Index, t::ITensor)
    if i1 == i2
        obs = op("obs", i1'') * op("obs" * ops[1], i1') * op("obs" * ops[2], i1)
        return (obs * t)[1]
    else
        obs1 = op("obs", i1') * op("obs" * ops[1], i1)
        obs2 = op("obs", i2') * op("obs" * ops[2], i2)
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

function apply_gate(::Type{Pure}, p::MPS, g::ProdLit; kwargs...)
    ops = Lit_to_ops(g, siteinds(p))
    return apply(ops, p; move_sites_back_between_gates=false, kwargs...)
end

function apply_gate(::Type{Mixed}, p::MPS, g::ProdLit; kwargs...)
    ops = Lit_to_ops(g, siteinds(p), "gate")
    return apply(ops, p; move_sites_back_between_gates=false, kwargs...)
end
