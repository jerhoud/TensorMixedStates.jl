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

Ffactor(::Mixed) = Gate(F)
Ffactor(::Pure) = F

function insertFfactors(a::Vector{Indexed{T}}) where T
    f = Ffactor(T())
    r = []
    fermion_idx = 0
    for l in reverse(a)
        idx = first(l.index)
        if fermion_idx ≠ 0
            append!(nls, (litF(i) for i in reverse(idx + 1:fermion_idx - 1)))
            if l.op.fermionic
                push!(nls, litF(idx))
                fermion_idx = 0
            else
                fermion_idx = idx
            end
        elseif l.op.fermionic && !l.op.dissipator
            fermion_idx = idx
        end
        push!(nls, l)
    end
    if fermion_idx ≠ 0
        append!(nls, (litF(i) for i in reverse(1:fermion_idx - 1)))
    end
    return reverse(r)
end

function simplify(a::SumOp{T, IndexOp}) where T
    subs = sort(reduce(vcat, map(t->sumsubs(simplify(t)), a.subs)))
    r = []
    c = 0
    p = []
    for s in subs
        ps = prodsubs(s)
        cs = prodcoef(s)
        if ps == p
            c += cs
        else
            if c ≠ 0
                push!(r, ProdOp{T, IndexOp}(c, p))
            end
            c = cs
            p = ps
        end
    end
    if c ≠ 0
        push!(r, ProdOp{T, IndexOp}(c, p))
    end
    return SumOp{T, IndexOp}(r)
end

distribute(a::Vector{<:Vector}) = a
distribute(a::Vector{<:Vector}, b::Vector, c::Vector...) = 
    if isempty(b)
        [[]]
    else
        r = Vector{Vector}(undef, length(a) * length(b))
        n = 1
        for i in a
            for j in b
                r[n] = vcat(i, [j])
                n += 1
            end
        end
        distribute(r, c...)
    end

distribute(a::Vector...) = distribute([[]], a...)

function simplify(a::ProdOp{T, IndexOp}) where T
    if a.coef == 0
        return a
    end
    args = map(simplify, a.subs)
    c = a.coef * prod(map(prodcoef, args))
    subs = reduce(vcat, map(prodsubs, args))
    if c == 0 || length(subs) <= 1
        return ProdOp{T, IndexOp}(c, subs)
    end
    sums = distribute(map(sumsubs, args)...)
    r = []

    c *= signature(sortperm(filter(fermionic, subs)))
    sort!(subs)



    return SumOp{T, IndexOp}(r)
end