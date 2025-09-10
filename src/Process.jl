export expand_index, simplify, insertFfactors, process

"""
    expand_index(op)

expand multiple site operators into simple site operators if possible. Also expands Dissipator and Evolver.
"""
expand_index(a) = apply_expr(expand_index, a)
expand_index(a::Indexed) = inner_index(a.op, a.index...)
function expand_index(a::Evolver)
    arg = expand_index(a.arg)
    return Left(arg) + Right(arg)
end

inner_index(a, idx...) = apply_expr(o->inner_index(o, idx...), a)
inner_index(a::Operator{T, 1}, idx...) where T = a(idx...)
inner_index(a::Proj, idx...) = a(idx...)
inner_index(a::Union{ExpOp{T, 1}, PowOp{T, 1}, SqrtOp{T, 1}}, idx...) where T = a(idx...)
inner_index(a::Operator, idx...) =
    if isnothing(a.expr)
        error("Cannot simplify matrix multiple site operators ($a)")
    else
        inner_index(a.expr, idx...)
    end
inner_index(a::TensorOp, idx...) = prod(tensor_apply(inner_index, a, idx...))
inner_index(a::Union{ExpOp, PowOp, SqrtOp}, idx...) = error("cannot simplify multiple site functionals ($a)")
inner_index(a::Gate) = error("Gate are not meant to be simplified ($a)")
function inner_index(a::Dissipator, idx...)
    a = inner_index(a.arg, idx...)
    aa = dag(a) * a
    Left(a) * Right(a) - 0.5 * Left(aa) - 0.5 * Right(aa)
end



function collect_sum(a::Vector{<:ExprIndexed{T}}) where T
    subs = sort(a; by=prodsubs)
    r = ExprIndexed{T}[]
    c = 0
    p = ExprIndexed{T}[]
    for s in subs
        ps = prodsubs(s)
        cs = prodcoef(s)
        if ps == p
            c += cs
        elseif c == 0
            c = cs
            p = ps
        elseif length(p) == 1 && length(ps) == 1 && p[1].index == ps[1].index && fermionic(p[1]) == fermionic(ps[1])
            p = [ Indexed{T, 1}(c*p[1].op+cs*ps[1].op, p[1].index) ]
            c = 1
        else
            push!(r, ProdOp{T, IndexOp}(c, p))
            c = cs
            p = ps
        end
    end
    if c ≠ 0
        push!(r, ProdOp{T, IndexOp}(c, p))
    end
    return SumOp{T, IndexOp}(r)
end

"""
    simplify(op)

simplify indexed operator
"""
simplify(a::SumOp{T, IndexOp}) where T =
    collect_sum(reduce(vcat, map(t->sumsubs(simplify(t)), a.subs)))

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
    r = ExprIndexed{T}[]
    for ind in a
        if !(ind isa Indexed{T, 1})
            error("Cannot simplify multiple site tensor ($ind)")
        end
        if ind.op == Id continue end
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

function simplify(a::ProdOp{T, IndexOp}) where T
    args = map(simplify, a.subs)   # get the simplified factors
    c = a.coef * prod(map(prodcoef, args)) # gather all product coefs
    subs = reduce(vcat, map(prodsubs, args)) # gather all product factors
    r = map(distribute(map(sumsubs, subs)...)) do p # develop the inner sums and loop over all resulting products
        cp = c * prod(map(prodcoef, p))             # gather all product coefs for each
        sp = reduce(vcat, map(prodsubs, p))         # gather all product factors for each
        if T == Pure
            cp *= signature(sortperm(filter(fermionic, sp); by=t->t.index))  # get the fermionic order sign
        else
            cp *= signature(sortperm(filter(fermionic_left, sp); by=t->t.index))
            cp *= signature(sortperm(filter(fermionic_right, sp); by=t->t.index))
        end
        ProdOp{T, IndexOp}(cp, collect_product(sort(sp; by=t->t.index))) # return the final product ordered with operators on the same site collapsed
    end
    return collect_sum(r)
end

simplify(a::Left) = simplify_left(simplify(a.arg))
simplify(a::Right) = simplify_right(simplify(a.arg))

simplify_left(a) = apply_expr(simplify_left, a)
simplify_left(a::Indexed) = Left(a.op)(a.index...)

simplify_right(a) = apply_expr(simplify_right, a)
simplify_right(a::Indexed) = Right(a.op)(a.index...)
simplify_right(a::ProdOp) = ProdOp{Mixed, IndexOp}(conj(a.coef), map(simplify_right, a.subs))

simplify(a::Indexed) = a

simplify(a::DagOp) = simplify_dag(simplify(a.arg))

simplify_dag(a::DagOp) = a.arg
simplify_dag(a::ProdOp{T, IndexOp}) where T = ProdOp{T, IndexOp}(conj(a.coef), map(simplify_dag, a.subs))
simplify_dag(a) = apply_expr(simplify_dag, a)
simplify_dag(a::Indexed) = dag(a.op)(a.index...)


Ffactor(::Indexed{Pure, 1}) = F
Ffactor(a::Indexed{Mixed, 1}) = Ffactor(a.op)
Ffactor(::Left) = Left(F)
Ffactor(::Right) = Right(F)
Ffactor(a) = Gate(F)

insertF(op::Gate) = Gate(op.arg * F)
insertF(op::Left) = Left(op.arg * F)
insertF(op::Right) = Right(op.arg * F)
insertF(op) = op * F

"""
    insertFfactors

insert the Jordan-Wigner F operators where needed
"""
insertFfactors(a::SumOp) = apply_expr(insertFfactors, a)
insertFfactors(a::ExprOp) = insertFfactors(prodcoef(a), prodsubs(a))

function insertFfactors(c::Number, v::Vector{<:ExprIndexed{T}}) where T
    local f::ExprOp{T, 1}
    r = ExprIndexed{T}[]
    fermion_idx = 0
    for l in reverse(v)
        idx = first(l.index)
        if fermion_idx ≠ 0
            append!(r, (f(i) for i in reverse(idx + 1:fermion_idx - 1)))
            push!(r, insertF(l.op)(idx))
            if fermionic(l)
                fermion_idx = 0
            else
                f = Ffactor(l)
                fermion_idx = idx
            end
        else
            push!(r, l)
            if fermionic(l)
                f = Ffactor(l)
                fermion_idx = idx
            end
        end
    end
    if fermion_idx ≠ 0
        append!(r, f(i) for i in reverse(1:fermion_idx - 1))
    end
    return ProdOp{T, IndexOp}(c, reverse(r))
end

"""
    process(op)

simplify indexed operator and insert Jordan-Wigner F factors where needed
"""
process(a::ExprIndexed) = insertFfactors(simplify(expand_index(a)))
process(a) = map(process, a)

process(a::ExprIndexed{Pure}, ::Pure, ::Any) = process(a)
process(a::ExprIndexed{Mixed}, ::Mixed, ::Any) = process(a)
process(a::ExprIndexed{Pure}, ::Mixed, f) = process(f(a))
process(a, type::PM, f) = map(t->process(a, type, f), a)