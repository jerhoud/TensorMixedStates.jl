# GenericOp Simplification

simplify(a::Union{Identity, JW_F, Proj, JW}) = a

simplify(a::Operator{R, 1}) where R = a
simplify(a::Operator) =
    if a.expr isa GenericOp
        simplify(a.expr)
    else
        error("cannot simplify implicit multisite operator")
    end

function simplify(a::AnyProd)
    s = simplify.(a.subs)
    simplify_prod(a.coef * prod(prodcoef.(s)), reduce(vcat, prodsubs.(s)))
end

simplify(a::AnySum) =
    simplify_sum(reduce(vcat, sumsubs.(simplify.(a.subs))))

simplify(a::TensorOp{N}) where N = TensorOp{N}(simplify.(a.subs))

simplify(a::PowOp) = simplify_pow(simplfiy(a.arg), a.expo)
simplify(a::ExpOp) = simplify_exp(simplify(a.arg))
simplify(a::DagOp) = simplify_dag(simplify(a.arg))

function simplify(a::Gate)
    sarg = simplify(a.arg)
    Left(sarg) * Right(sarg)
    simplify_prod(1, [simplify_l(sarg), simplify_r(sarg)])
end

function simplify(a::Dissipator)
    sarg = simplify(a.arg)
    darg = simplify_dag(sarg)
    daga = simplify_prod(1, [darg, sarg])
    simplify_sum(simplify_prod(1, [simplify_l(sarg), simplify_r(sarg)]) -
        0.5 * simplify_l(daga) -
        0.5 * simplify_r(daga))
end

simplify(a::Left) = simplify_l(simplify(a.arg))
simplify(a::Right) = simplify_r(simplify(a.arg))

# IndexepOp simplification

function simplify(a::Evolver)
    sarg = simplify(a.arg)
    simplify_sum([simplify_l(sarg), simplify_r(sarg)])
end

simplify(a::Indexed) =
    simplify_ind(simplify(a.op), a.index...)

# GenericOp with indices Simplification

simplify_ind(::Identity, ::Int) = IndexedId
simplify_ind(a::Union{JW_F, Proj, JW, Left{1}, Right{1}}, index) = a(index)
simplify_ind(a::Union{Operator, ExpOp, PowOp, DagOp}, index) = a(index)
simplify_ind(a::Left{N}, index...) = simplify_l(simplify_ind(a.arg, index...))
simplify_ind(a::Right{N}, index...) = simplify_r(simplify_ind(a.arg, index...))


simplify_ind(a::SumGenOp, index...) =
    simplify_sum(reduce(vcat, (map(x->sumsubs(simplify_ind(x, index...)), a.subs))))

function simplify_ind(a::ProdGenOp, index...)
    s = map(x->simplify_ind(x, index...), a.subs)
    simplify_prod(a.coef * prod(prodcoef.(s)), reduce(vcat, prodsubs.(s)))
end

simplify_ind(a::TensorOp, index...) = simplify_prod(1, tensor_apply(simplify_ind, a, index...))

# helpers

is_involution(::Identity) = true
is_involution(::JW_F) = true
is_involution(a::Operator) = a.type == involution_op
is_involution(::Any) = false

simplify_pow(arg::GenericOp{R, N}, expo) where {R, N} =
    if expo == 0
        return Identity(arg)
    elseif expo == 1
        return arg
    elseif is_involution(arg)
        simplify_pow(arg, expo - 2 * floor(expo / 2))
    elseif N > 1
        error("cannot simplify power of multisite operator")
    else
        return PowOp(arg, expo)
    end

function simplify_exp(a::GenericOp{R, N}) where {R, N}
    c = prodcoef(a)
    s = prodsubs(a)
    if length(s) == 1 && is_involution(s[1])
        cosh(c) * Identity(s[1]) + sinh(c) * s[1]
    elseif N > 1
        error("cannot simplify exponential of multisite operator")
    else
        exp(a)
    end
end 


simplify_dag(a::DagOp) = a.arg
simplify_dag(a::ProdGenOp) = ProdGenOp(conj(a.coef), reverse(simplify_dag.(a.subs)))
simplify_dag(a::SumGenOp) = SumGenOp(simplify_dag.(a.subs))
simplify_dag(a::TensorOp{N}) where N = TensorOp{N}(simplify_dag.(a.subs))
simplify_dag(a::Operator) =
    if a.type == involution_op || a.type == selfadjoint_op
        a 
    else
        dag(a)
    end
simplify_dag(a::JW) = JW(simplify_dag(a.arg))
simplify_dag(a::Left) = Left(simplify_dag(a.arg))
simplify_dag(a::Right) = Right(simplify_dag(a.arg))
simplify_dag(a::PowOp) = PowOp(simplify_dag(a.arg), conj(a.expo))
simplify_dag(a::ExpOp) = ExpOp(simplify_dag(a.arg))
simplify_dag(a::Union{Identity, JW_F}) = a
simplify_dag(a::Proj) = dag(a)


simplify_l(a::GenericOp) = Left(a)

simplify_l(a::ProdIndOp{Pure}) = ProdIndOp{Mixed}(a.coef, simplify_l.(a.subs)) 
simplify_l(a::SumIndOp{Pure}) = SumIndOp{Mixed}(simplify_l.(a.subs))
simplify_l(a::Indexed{Pure, 1}) = simplify_l(a.op)(a.index...)
simplify_l(a::Multi_F{Pure}) = Multi_F{Mixed}(a.start, a.stop, true, false)

simplify_r(a::Identity) = Left(a)
simplify_r(a::GenericOp) = Right(a)

simplify_r(a::ProdIndOp{Pure}) = ProdIndOp{Mixed}(conj(a.coef), simplify_r.(a.subs)) 
simplify_r(a::SumIndOp{Pure}) = SumIndOp{Mixed}(simplify_r.(a.subs))
simplify_r(a::Indexed{Pure, 1}) = simplify_r(a.op)(a.index...)
simplify_r(a::Multi_F{Pure}) = Multi_F{Mixed}(a.start, a.stop, false, true)


function simplify_sum(v::Vector{<:GenericOp{R, N}}) where {R, N}
    subs = sort(v; by=prodsubs)
    r = GenericOp{R, N}[]
    c = 0
    p = []
    for s in subs
        ps = prodsubs(s)
        cs = prodcoef(s)
        if ps == p
            c += cs
        else
            if c ≠ 0
                push!(r, ProdGenOp{R, N}(c, p))
            end
            c = cs
            p = ps
        end
    end
    if c ≠ 0
        push!(r, ProdGenOp{R, N}(c, p))
    end
    if isempty(r)
        return 0 * Identity(subs[1])
    else
        return SumGenOp{R, N}(r)
    end
end

function simplify_sum(v::Vector{<:IndexedOp{R}}) where {R}
    subs = sort(v; by=prodsubs)
    r = IndexedOp{R}[]
    c = 0
    p = []
    for s in subs
        ps = prodsubs(s)
        cs = prodcoef(s)
        if ps == p
            c += cs
        elseif c == 0
            c = cs
            p = ps
        elseif length(p) == 1 && length(ps) == 1 && p[1].index == ps[1].index
            p = [ Indexed{R, 1}(c * p[1].op + cs * ps[1].op, p[1].index) ]
            c = 1
        else
            push!(r, ProdIndOp{R}(c, p))
            c = cs
            p = ps
        end
    end
    if c ≠ 0
        push!(r, ProdIndOp{R}(c, p))
    end
    if isempty(r)
        return Identity(subs[1])
    else
        return SumIndOp{R}(r)
    end
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

pow_base(a) = a
pow_base(a::PowOp) = a.arg
pow_expo(a) = 1
pow_expo(a::PowOp) = a.expo

function simplify_prod(c::Number, v::Vector{<:GenericOp{R, N}}) where {R, N}
    id = Identity(v[1])
    if c == 0
        return 0 * id
    end
    r = GenericOp{R, N}[]
    b = id
    e = 1
    for x in v
        bx = pow_base(x)
        ex = pow_expo(x)
        if bx == b
            e += ex
        else
            sp = simplify_pow(b, e)
            if sp ≠ id
                push!(r, sp)
            end
            b = bx
            e = ex
        end
    end
    sp = simplify_pow(b, e)
    if sp ≠ id
        push!(r, sp)
    end
    if isempty(r)
        return c * id
    else
        return ProdGenOp{R, N}(c, r)
    end
end

function simplify_prod(c::Number, v::Vector{<:IndexedOp{R}}) where R
    id = Identity(v[1])
    if c == 0
        return 0 * id
    end
    r = IndexedOp{R}[]
    b = id
    e = 1
    for x in v
        bx = pow_base(x)
        ex = pow_expo(x)
        if bx == b
            e += ex
        else
            sp = simplify_pow(b, e)
            if sp ≠ id
                push!(r, sp)
            end
            b = bx
            e = ex
        end
    end
    sp = simplify_pow(b, e)
    if sp ≠ id
        push!(r, sp)
    end
    if isempty(r)
        return c * id
    else
        return ProdGenOp{R, N}(c, r)
    end
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

function simplify(a::ProdOp{T, IndexOp}) where T
    args = map(simplify, a.subs)   # get the simplified factors
    c = a.coef * prod(map(prodcoef, args)) # gather all product coefs
    subs = reduce(vcat, map(prodsubs, args)) # gather all product factors
    r = map(distribute(map(sumsubs, subs)...)) do p # develop the inner sums and loop over all resulting products
        cp = c * prod(map(prodcoef, p))             # gather all product coefs for each
        sp = reduce(vcat, map(prodsubs, p))         # gather all product factors for each
        cp *= signature(sortperm(filter(fermionic, sp); by=t->t.index))  # get the fermionic order sign
        ProdOp{T, IndexOp}(cp, collect_product(sort(sp; by=t->t.index)))    # return the final product ordered with operators on the same site collapsed
    end
    return SumOp{T, IndexOp}(collect_sum(r))
end
