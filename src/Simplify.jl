export expand_index, simplify, insertFfactors

expand_index(a) = apply_expr(expand_index, a)
expand_index(a::Indexed{T, 1}) where T =
    if has_dissipator(a)
        inner_index(a.op, a.index...)
    else
        a
    end
expand_index(a::Indexed{T, N}) where {T, N} = inner_index(a.op, a.index...)

inner_index(a, idx...) = apply_expr(o->inner_index(o, idx...), a)
inner_index(a::Operator{T, 1}, idx...) where T = a(idx...)
inner_index(a::Union{ExpOp{T, 1}, PowOp{T, 1}, SqrtOp{T, 1}}, idx...) where T = a(idx...)
inner_index(a::Operator, idx...) = inner_index(a.expr, idx...)
inner_index(a::TensorOp, idx...) = prod(tensor_apply(inner_index, a, idx...))
inner_index(::Union{ExpOp, PowOp, SqrtOp}, idx...) = error("cannot simplify multiple site functionals")
function inner_index(a::Dissipator{1}, idx...)
    aa = dag(a.arg) * a.arg
    if fermionic(a.arg)
        return Gate(a.arg)(idx...) + (Left(aa) + Right(aa))(idx...)
    else
        return (Gate(a.arg) + Left(aa) + Right(aa))(idx...)
    end
end
function inner_index(a::Dissipator, idx...)
    a = inner_index(a.arg, idx...)
    aa = dag(a) * a
    Gate(a) + Left(aa) + Right(aa)
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

function collect_product(a::Vector{<:ExprIndexed{T}}) where T
    if isempty(a)
        return []
    end
    local o::ExprOp{T, 1}
    idx = 0
    r = []
    for ind in a
        if !(ind isa Indexed{T, 1})
            idx = 0
            continue
        end
        i = ind.index[1]
        if i == idx
            o *= ind.op
        else
            if idx ≠ 0
                push!(r, o(idx))
            end
            idx = i
            o = ind.op
        end
    end
    if idx ≠ 0
        push!(r, o(idx))
    end
    return r
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
    ni = filter(t->!isa(t, Indexed), r)
    f = filter(t->isa(t, Indexed) && fermionic(t.op), r)
    nf = filter(t->isa(t, Indexed) && !fermionic(t.op), r)
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
    args = map(simplify, a.subs)   # get the simplified factors
    c = a.coef * prod(map(prodcoef, args)) # gather all product coefs
    subs = reduce(vcat, map(prodsubs, args)) # gather all product factors
    if c == 0 || length(subs) <= 1           # early bail if it is simple
        return ProdOp{T, IndexOp}(c, subs)
    end
    r = map(distribute(map(sumsubs, subs)...)) do p # develop the inner sums and loop over all resulting products
        cp = c * prod(map(prodcoef, p))             # gather all product coefs for each
        sp = reduce(vcat, map(prodsubs, p))         # gather all product factors for each
        cp *= signature(sortperm(filter(fermionic, sp)))  # get the fermionic order sign
        ProdOp{T, IndexOp}(cp, collect_product(sort(sp)))    # return the final product ordered with operators on the same site collapsed
    end
    return SumOp{T, IndexOp}(r)
end

function simplify(a::Evolve)
    args = sumsubs(simplify(a.arg))
    r = []
    for p in args
        c = prodcoef(p)
        s = prodsubs(p)
        push!(r, ProdOp{Mixed, IndexOp}(c, map(ind->Left(ind.op)(ind.index...), s)))
        push!(r, ProdOp{Mixed, IndexOp}(c, map(ind->Right(ind.op)(ind.index...), s)))
    end
    return SumOp{Mixed, IndexOp}(r)
end

function simplify(a::Gate{IndexOp})
    arg = simplify(a.arg)
    if arg isa SumOp
        error("Gate of sums not implemented: Gate($arg)")
    end
    c = abs2(prodcoef(arg))
    s = map(ind->Gate(ind.op)(ind.index...), prodsubs(arg))
    return ProdOp{Mixed, IndexOp}(c, s)
end

simplify(a::Indexed) = a

simplify(a::DagOp) = simplify_dag(a.arg)

simplify_dag(a::DagOp) = simplify(a.arg)
simplify_dag(a::ProdOp) = simplify(ProdOp(conj(a.coef), reverse(map(dag, a.subs))))
simplify_dag(a) = simplify(apply_expr(dag, a))

Ffactor(::Indexed{Pure, 1}) = F
Ffactor(a::Indexed{Mixed, 1}) = Ffactor(a.arg)
Ffactor(::Left) = Left(F)
Ffactor(::Right) = Right(F)
Ffactor(a) = Gate(F)

insertFfactors(a::SumOp) = apply_expr(insertFfactors, a)
function insertFfactors(a::ProdOp)
    local f::ExrpOp
    r = []
    fermion_idx = 0
    for l in reverse(a)
        idx = first(l.index)
        if fermion_idx ≠ 0
            append!(r, (f(i) for i in reverse(idx + 1:fermion_idx - 1)))
            push!(r, (l.op * f)(idx))
            if fermionic(l)
                fermion_idx = 0
            else
                fermion_idx = idx
            end
        elseif fermionic(l)
            push!(r, l)
            fermion_idx = idx
        end
    end
    if fermion_idx ≠ 0
        append!(r, f(i) for i in reverse(1:fermion_idx - 1))
    end
    return ProdOp(a.coef, reverse(r))
end

