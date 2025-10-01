export simplify

# Simplification

simplify(a::Union{Identity, JW_F, Proj, JW, Operator}) = a

simplify(a::ScalarOp; kwargs...) = a.coef * simplify(a.arg; kwargs...)
simplify(a::ProdOp; kwargs...) = simplify_prod(map(x->simplify(x; kwargs...), a.subs))
simplify(a::SumOp; kwargs...) = simplify_sum(map(x->simplify(x; kwargs...), a.subs))
simplify(a::TensorOp{N}) where N = TensorOp{N}(simplify.(a.subs))

simplify(a::PowOp) = simplify_pow(simplify(a.arg), a.expo)
simplify(a::ExpOp) = simplify_exp(simplify(a.arg))
simplify(a::DagOp) = simplify_dag(simplify(a.arg))

function simplify(a::Gate)
    sarg = simplify(a.arg)
    simplify_prod([simplify_l(sarg), simplify_r(sarg)])
end

function simplify(a::Dissipator)
    sarg = simplify(a.arg)
    darg = simplify_dag(sarg)
    daga = simplify_prod([darg, sarg])
    simplify_sum([simplify_prod([simplify_l(sarg), simplify_r(sarg)]), -
        0.5 * simplify_l(daga), -
        0.5 * simplify_r(daga)])
end

simplify(a::Left) = simplify_l(simplify(a.arg))
simplify(a::Right) = simplify_r(simplify(a.arg))

function simplify(a::Evolver)
    sarg = simplify(a.arg)
    simplify_sum([simplify_l(sarg), simplify_r(sarg)])
end

simplify(a::AtIndex; kwargs...) =
    simplify_ind(simplify(a.op), a.index...; kwargs...)

# Simplification with index

simplify_ind(::Identity, ::Int; kwargs...) = Id(1)
simplify_ind(a::ScalarOp, index...; kwargs...) = a.coef * simplify_ind(a.arg, index...; kwargs...)
simplify_ind(a::Union{JW_F, Proj, JW, Left{1}, Right{1}}, index; kwargs...) = a(index)
simplify_ind(a::Union{ExpOp, PowOp, DagOp}, index; kwargs...) = a(index)
simplify_ind(a::Left, index...; kwargs...) = simplify_l(simplify_ind(a.arg, index...; kwargs...))
simplify_ind(a::Right, index...; kwargs...) = simplify_r(simplify_ind(a.arg, index...; kwargs...))

simplify_ind(a::Operator{1}, index; kwargs...) =
    if a.type == fermionic_op
        if index > 1
            Multi_F{Pure}(1, index-1, false, false) * JW(a)(index)
        else
            JW(a)(index)
        end
    else
        a(index)
    end

simplify_ind(a::Operator, index...; expand = false) =
    if expand && a.expr isa Op
        simplify(a.expr; expand)(index...)
    else
        a(index...)
    end

simplify_ind(a::SumOp, index...; kwargs...) = simplify_sum(map(x->simplify_ind(x, index...; kwargs...), a.subs))
simplify_ind(a::ProdOp, index...; kwargs...) = simplify_prod(map(x->simplify_ind(x, index...; kwargs...), a.subs))
simplify_ind(a::TensorOp, index...; kwargs...) = simplify_prod(tensor_apply(simplify_ind, a, index...; kwargs...))

# helpers

is_involution(::Identity) = true
is_involution(::JW_F) = true
is_involution(a::Operator) = a.type == involution_op
is_involution(::Op) = false

simplify_pow(a::GenericOp{Pure, N}, expo) where N =
    if expo == 0
        return MakeIdentity(a)
    elseif expo == 1
        return a
    elseif is_involution(a) && expo ≥ 2
        simplify_pow(a, expo - 2 * floor(expo / 2))
    elseif N > 1
        error("cannot simplify power of multisite operator")
    else
        return PowOp(a, expo)
    end

function simplify_exp(a::GenericOp{Pure, N}) where N
    c = prodcoef(a)
    s = prodsubs(a)
    if length(s) == 1 && is_involution(s[1])
        cosh(c) * MakeIdentity(s[1]) + sinh(c) * s[1]
    elseif N > 1
        error("cannot simplify exponential of multisite operator")
    else
        exp(a)
    end
end 


simplify_dag(a::DagOp) = a.arg
simplify_dag(a::ScalarOp) = conj(a.coef) * simplify_dag(a.arg)
simplify_dag(a::ProdOp) = ProdOp(reverse(simplify_dag.(a.subs)))
simplify_dag(a::SumOp) = SumOp(simplify_dag.(a.subs))
simplify_dag(a::TensorOp{N}) where N = TensorOp{N}(simplify_dag.(a.subs))
simplify_dag(a::Operator) =
    if a.type == involution_op || a.type == selfadjoint_op
        a 
    else
        dag(a)
    end
simplify_dag(a::JW) = JW(simplify_dag(a.arg))
simplify_dag(a::PowOp) = PowOp(simplify_dag(a.arg), conj(a.expo))
simplify_dag(a::ExpOp) = ExpOp(simplify_dag(a.arg))
simplify_dag(a::Union{Identity, JW_F}) = a
simplify_dag(a::Proj) = dag(a)


simplify_l(a::GenericOp{Pure}) = Left(a)
simplify_l(a::ScalarOp{Pure}) = a.coef * simplify_l(a.arg) 
simplify_l(a::ProdOp{Pure, Indexed}) = ProdOp(simplify_l.(a.subs)) 
simplify_l(a::SumOp{Pure, Indexed}) = SumOp(simplify_l.(a.subs))
simplify_l(a::AtIndex{Pure}) = simplify_l(a.op)(a.index...)

simplify_r(a::GenericOp{Pure}) = Right(a)
simplify_r(a::Identity) = Left(a)
simplify_r(a::ScalarOp{Pure}) = conj(a.coef) * simplify_r(a.arg) 
simplify_r(a::ProdOp{Pure, Indexed}) = ProdOp(simplify_r.(a.subs)) 
simplify_r(a::SumOp{Pure, Indexed}) = SumOp(simplify_r.(a.subs))
simplify_r(a::AtIndex{Pure}) = simplify_r(a.op)(a.index...)

simplify_sum(v::Vector) = simplify_core_sum(reduce(vcat, sumsubs.(v)))

function simplify_core_sum(v::Vector{<:Op{R, T, N}}) where {R, T, N}
    subs = sort(v; by=scalararg)
    r = Op{R, T, N}[]
    c = 0
    o = MakeIdentity{R, T, N}()
    for s in subs
        nc = scalarcoef(s)
        no = scalararg(s)
        if no == o
            c += nc
        elseif c == 0
            c = nc 
            o = no
        elseif T == Indexed && o isa AtIndex && no isa AtIndex && o.index == no.index
            o = AtIndex(simplify_sum([c * o.op, nc * no.op]), o.index)
            c = 1
            if o isa ScalarOp
                c = o.coef
                o = o.arg
            end
        else
            push!(r, c * o)
            c = nc
            o = no
        end
    end
    if c ≠ 0
        push!(r, c * o)
    end
    return SumOp(r)
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

pow_base(a::Op) = a
pow_base(a::PowOp) = a.arg
pow_expo(a::Op) = 1
pow_expo(a::PowOp) = a.expo

commut(::Op) = true
commut(::JW) = false
commut(a::DagOp) = commut(a.arg)

simplify_prod(v::Vector) =
     simplify_core_prod(prod(scalarcoef.(v)), reduce(vcat, prodsubs.(v)))

function simplify_core_prod(c::Number, v::Vector{<:GenericOp{Pure, N}}) where N
    id = MakeIdentity(v[1])
    if c == 0
        return 0 * id
    end
    r = GenericOp{Pure, N}[]
    b = id
    e = 1
    f = false
    for x in v
        bx = pow_base(x)
        ex = pow_expo(x)
        if bx isa JW_F
            f = !f
            continue
        end
        if f && !commut(bx) && isodd(ex)
            c = -c
        end
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
    if f
        push!(r, F)
    end
    return c * ProdOp(r)
end

function simplify_core_prod(c::Number, v::Vector{<:GenericOp{Mixed, N}}) where N
    id = MakeIdentity(v[1])
    if c == 0
        return 0 * id
    end
    r = GenericOp{Mixed, N}[]
    larg = map(x -> x.arg, filter(x -> x isa Left, v))
    if !isempty(larg)
        left = simplify_l(simplify_core_prod(1, larg))
        if left.arg ≠ id.arg
            push!(r, left)
        end
    end
    rarg = map(x -> x.arg, filter(x -> x isa Right, v))
    if !isempty(rarg)
        right = simplify_r(simplify_core_prod(1, rarg))
        if right.arg ≠ id.arg
            push!(r, right)
        end
    end
    return c * ProdOp(r)
end

orderprod(a::AtIndex, b::AtIndex) where R =
    if a.index == b.index
        [ simplify_core_prod(1, [a, b])(a.index...) ]
    elseif min(a.index...) > max(b.index)
        [b, a]
    else
        []
    end

function orderprod(a::AtIndex{R, N}, b::Multi_F{R}) where {R, N}
    i = min(a.index...)
    j = max(a.index...)
    if i < b.start || (N > 1 && i == b.start)
        []
    elseif i > b.stop
        [b, a]
    elseif N == 1
        [Multi_F{R}(b.start, i-1, b.left, b.right), a, Multi_F{R}(i, i, b.left, b.right), Multi_F{R}(i + 1, b.stop, b.left, b.right)]
    else
        [Multi_F{R}(b.start, i-1, b.left, b.right), a, Multi_F{R}(i, b.stop, b.left, b.right)]
    end
end

function orderprod(b::Multi_F{R}, a::AtIndex{R, 1}) where R
    i = min(a.index...)
    j = max(a.index...)
    if i < b.start || (N > 1 && i == b.start)
        [a, b]
    elseif i > b.stop
        []
    elseif N == 1
        [Multi_F{R}(b.start, i-1, b.left, b.right), a, Multi_F{R}(i, i, b.left, b.right), Multi_F{R}(i + 1, b.stop, b.left, b.right)]
    else
        [Multi_F{R}(b.start, i-1, b.left, b.right), a, Multi_F{R}(i, b.stop, b.left, b.right)]
    end
end

orderprod(a::Multi_F{R}, b::Multi_F{R}) where R = 
    if a.left == b.left && a.right == b.right && (a.stop == b.start-1 || b.stop == a.start-1)
        [ Multi_F{R}(min(a.start, b.start), max(a.stop, b.stop), a.left, a.right)]
    elseif a.stop < b.start
        []
    elseif b.stop < a.start
        [b, a]
    else
        i = min(a.stop, b.stop)
        j = max(a.start, b.start)
        [
            Multi_F{R}(a.start, min(a.stop, i-1), a.left, a.right), Multi_F{R}(b.start, bmin(b.stop, i-1), b.left, b.right),
            Multi_F{R}(i, j, a.left ⊻ b.left, a.right ⊻ b.right),
            Multi_F{R}(max(b.start, j+1), b.stop, b.left, b.right), Multi_F{R}(max(a.start, j+1), a.stop, a.left, a.right)
        ]
    end

function simplify_core_prod(c::Number, v::Vector{<:IndexedOp{R}}) where R
    id = MakeIdentity(v[1])
    if c == 0
        return 0 * id
    end
    s = map(distribute(sumsubs.(v)...)) do p # expand all inner sums
        cp = c * prod(scalarcoef.(p))
        if cp == 0
            return 0 * id
        end
        r = reduce(vcat, prodsubs.(p))
        change = true
        while change
            change = false
            start = true
            nr = IndexedOp{R}[]
            left = id
            for right in r
                if left == id
                    left = right
                    if !start
                        change = true
                    end
                else
                    f, ll, rr  = orderprod(left, right)
                    if f 
                        change = true
                        append!(nr, ll)
                        left = rr
                    else
                        push!(nr, left)
                        left = right
                    end
                end
                start = false
            end
            r = nr
        end
        return cp * ProdOp(r)
    end
    return simplify_sum(s)
end
