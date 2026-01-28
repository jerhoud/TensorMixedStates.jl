export simplify, removeMulti

"""
    simplify(op::Op; expand=false)

simplifies an operator, this is used internally by functions creating MPOs.
expand tells wheteher to expand explicit operators with their definition
"""
function simplify end

# Simplifications for collections of Operators

simplify(a; kwargs...) = map(x->simplify(x; kwargs...), a)


# Simplification principles
# One pass
# expand definitions of Gate and Evolver
# expand as go we in, simplify as we go out
# the final result should be a sum of products
# of one site (possibly complex) indexed operators with different indices in ascending order




# Simplifications for both Generic and Indexed operators

simplify(a::ScalarOp; kwargs...) = a.coef * simplify(a.arg; kwargs...)
simplify(a::ProdOp; kwargs...) = simplify_prod(map(x->simplify(x; kwargs...), a.subs))
simplify(a::SumOp; kwargs...) = simplify_sum(map(x->simplify(x; kwargs...), a.subs))
simplify(a::TensorOp{N}) where N = TensorOp{N}(simplify.(a.subs))


# Simplification of Generic Operators

simplify(a::Union{Identity, JW_F, Proj, JW, Operator, Multi_F, SetState}) = a

simplify(a::PowOp) = simplify_pow(simplify(a.arg), a.expo)
simplify(a::ExpOp) = simplify_exp(simplify(a.arg))
simplify(a::DagOp) = simplify_dag(simplify(a.arg))

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


# Simplifications of Indexed Operators

simplify(a::Multi_F) = a

function simplify(a::Gate)
    sarg = simplify(a.arg)
    simplify_prod([simplify_l(sarg), simplify_r(sarg)])
end

function simplify(a::Evolver)
    sarg = simplify(a.arg)
    simplify_sum([simplify_l(sarg), simplify_r(sarg)])
end

simplify(a::AtIndex; kwargs...) =
    simplify_ind(simplify(a.op), a.index...; kwargs...)

reindex(::Identity, ::Int) = Id(1)
reindex(op::GenericOp, i::Int...) = op(i...)

# Simplification with index
# transmit indexation as deep as possible
# to develop tensors, transform fermionic operators with JW, expand operators with expand = true


simplify_ind(a::ScalarOp, index...; kwargs...) = a.coef * simplify_ind(a.arg, index...; kwargs...)
simplify_ind(a::Union{Identity, JW_F, Proj, JW, SetState}, index; kwargs...) = a(index)
simplify_ind(a::ExpOp, index...; kwargs...) = a(index...)
simplify_ind(a::PowOp, index...; kwargs...) = simplify_pow(simplify_ind(a.arg, index...; kwargs), a.expo)
simplify_ind(a::DagOp, index...; kwargs...) = simplify_dag(simplify_ind(a.arg, index...; kwargs...))
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
        simplify_ind(a.expr, index...; expand)
    else
        a(index...)
    end

simplify_ind(a::SumOp, index...; kwargs...) = simplify_sum(map(x->simplify_ind(x, index...; kwargs...), a.subs))
simplify_ind(a::ProdOp, index...; kwargs...) = simplify_prod(map(x->simplify_ind(x, index...; kwargs...), a.subs))
simplify_ind(a::TensorOp, index...; kwargs...) = simplify_prod(tensor_apply(simplify_ind, a, index...; kwargs...))


# helpers for Generic Operators

is_involution(::Identity) = true
is_involution(::JW_F) = true
is_involution(a::Operator) = a.type == involution_op
is_involution(::Op) = false


# power simplification
# power of Generic:  X^0 => Id, X^1=> X, X^2 => Id
simplify_pow(a::GenericOp{Pure}, expo) =
    if expo == 0
        return MakeIdentity(a)
    elseif expo == 1
        return a
    elseif is_involution(a) && expo ≥ 2
        simplify_pow(a, expo - 2 * floor(expo / 2))
    else
        return PowOp(a, expo)
    end

simplify_pow(a::ScalarOp{Pure}, expo) =
    if isinteger(expo) || a.coef > 0
        a.coef^expo * simplify_pow(a.arg, expo)
    else
        (-a.coef)^expo * PowOp(-arg, expo)
    end


# power of Indexed: replace with a product if possible X(1)^2 => X(1) * X(1)
simplify_pow(a::IndexedOp, expo) =
    if expo == 0
        return 0 * MakeIdentity(a)
    elseif expo == 1
        return a
    elseif isinteger(expo)
        simplify_prod(fill(a, Integer(expo)))
    else
        a^expo
    end

simplify_pow(a::AtIndex, expo) =
    simplify_pow(a.op, expo)(a.index...)


# simplify exp : exp(0) => Id and exp(3X) => cosh(3)Id + sinh(3)X
function simplify_exp(a::GenericOp{Pure, N}) where N
    c = scalarcoef(a)
    s = prodsubs(a)
    if length(s) == 1 && is_involution(s[1])
        simplify_sum([cosh(c) * MakeIdentity(s[1]), sinh(c) * s[1]])
    else
        exp(a)
    end
end 


# dag simplification : transmit the dag as deep as possible to allow
# simplify operators according to type : dag(Id) = Id, dag(F) = F, dag(X) = X

simplify_dag(a::DagOp) = a.arg
simplify_dag(a::ScalarOp) = conj(a.coef) * simplify_dag(a.arg)
simplify_dag(a::Operator) =
    if a.type == involution_op || a.type == selfadjoint_op
        a 
    else
        dag(a)
    end
simplify_dag(a::PowOp) = PowOp(simplify_dag(a.arg), conj(a.expo))
simplify_dag(a::ExpOp) = ExpOp(simplify_dag(a.arg))
simplify_dag(a::Union{Identity, JW_F}) = a
simplify_dag(a::Union{Proj, JW}) = dag(a)


simplify_dag(a::ProdOp) = simplify_prod(reverse(simplify_dag.(a.subs)))
simplify_dag(a::SumOp) = simplify_sum(simplify_dag.(a.subs))
simplify_dag(a::TensorOp{N}) where N = TensorOp{N}(simplify_dag.(a.subs))


simplify_dag(a::AtIndex) = reindex(simplify_dag(a.op), a.index...)
simplify_dag(a::Multi_F) = a


# Left simplification
# Generic => just get the scalar factor out
# Indexed => go as deep as possible

simplify_l(a::GenericOp{Pure}) = Left(a)
simplify_l(a::ScalarOp{Pure}) = a.coef * simplify_l(a.arg)

simplify_l(a::ProdOp{Pure, Indexed}) = ProdOp(simplify_l.(a.subs)) 
simplify_l(a::SumOp{Pure, Indexed}) = SumOp(simplify_l.(a.subs))
simplify_l(a::AtIndex{Pure}) = reindex(simplify_l(a.op), a.index...)
simplify_l(a::Multi_F{Pure}) = Multi_F{Mixed}(a.start, a.stop, true, false)


# Right simplification
# Generic => get the scalar factor out and simplifies Right(Id) => Left(Id)
# Indexed => go as deep as possible

simplify_r(a::GenericOp{Pure}) = Right(a)
simplify_r(a::Identity) = Left(a)
simplify_r(a::ScalarOp{Pure}) = conj(a.coef) * simplify_r(a.arg) 

simplify_r(a::ProdOp{Pure, Indexed}) = ProdOp(simplify_r.(a.subs)) 
simplify_r(a::SumOp{Pure, Indexed}) = SumOp(simplify_r.(a.subs))
simplify_r(a::AtIndex{Pure}) = reindex(simplify_r(a.op), a.index...)
simplify_r(a::Multi_F{Pure}) = Multi_F{Mixed}(a.start, a.stop, false, true)


# sum simplification
# flatten out inner sums, order terms, collect identical terms and remove nuls
# X + (Y + Z) => X + Y + Z, X + Y + X => 2X + Y, X - X => 0
# in Indexed sums gather terms with same indices : X(1) + Y(1) => (X+Y)(1)

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
            o = reindex(simplify_sum([c * o.op, nc * no.op]), o.index...)
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


# product Simplifications
# first flatten out inner products X * (Y * Z) => X * Y * Z

simplify_prod(v::Vector) =
     simplify_core_prod(prod(scalarcoef.(v)), reduce(vcat, prodsubs.(v)))

# Generic Pure product
# gather identical factors X * X => X^2
# simplify powers using operator types Id^2 => Id, F^2 => Id, X^2 => Id
# simplify F and put them at the end with the correct sign (anticommut with JW, commut with the others)
# F*X = X*F, F*JW(C) => -JW(C)*F, F*F = Id

pow_base(a::Op) = a
pow_base(a::PowOp) = a.arg
pow_expo(a::Op) = 1
pow_expo(a::PowOp) = a.expo

commut(::Op) = true
commut(::JW) = false
commut(a::DagOp) = commut(a.arg)

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

# Generic Mixed products
# (that is only Left and Right factors)
# Gather Left and Right together
# Left(X)*Right(Y)*Left(Z) => Left(X*Z)*Right(Y)

function simplify_core_prod(c::Number, v::Vector{<:GenericOp{Mixed, N}}) where N
    id = MakeIdentity(v[1])
    if c == 0
        return 0 * id
    end
    larg = map(x -> x.arg, filter(x -> x isa Left, v))
    rarg = map(x -> x.arg, filter(x -> x isa Right, v))
    # test whether factors other than Left and Right are present (like sums)
    if length(larg) + length(rarg) ≠ length(v)
        return c * ProdOp(v)
    end
    r = GenericOp{Mixed, N}[]
    if !isempty(larg)
        left = simplify_l(simplify_prod(larg))
        if left ≠ id
            push!(r, left)
        end
    end
    if !isempty(rarg)
        right = simplify_r(simplify_prod(rarg))
        if right ≠ id
            push!(r, right)
        end
    end
    return c * ProdOp(r)
end

# Indexed product simplification
# first expand inner sums with distribute
# use orderedprod to reorder and simplify / expand factors
# in a kind of bublesort in simplify_core_prod

# order, gather and simplify factors X(1)Z(2)Y(1)Id(3) => (X*Y)(1)*Z(2)
# do the right things with Multi_F (glue, split, reduce)
# so that C(3)C(5) => Multi_F(1,2)JW(C)(3)Multi_F(1,4)JW(C)(5) => (JW(C)*F)(3)*F(4)*JW(C)(5)

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

orderprod(a::AtIndex, b::AtIndex) =
    if a.index == b.index
        [ reindex(simplify_prod([a.op, b.op]), a.index...) ]
    elseif min(a.index...) > max(b.index...)
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

function orderprod(b::Multi_F{R}, a::AtIndex{R, N}) where {R, N}
    i = min(a.index...)
    j = max(a.index...)
    if i < b.start || (N > 1 && i == b.start)
        [a, b]
    elseif i > b.stop
        []
    elseif N == 1
        [Multi_F{R}(b.start, i-1, b.left, b.right), Multi_F{R}(i, i, b.left, b.right), a, Multi_F{R}(i + 1, b.stop, b.left, b.right)]
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
        i = max(a.start, b.start)
        j = min(a.stop, b.stop)
        if a.left == b.left && a.right == b.right
            m = MakeIdentity(a)
        else
            m = Multi_F{R}(i, j, a.left ⊻ b.left, a.right ⊻ b.right)
        end
        [
            Multi_F{R}(a.start, min(a.stop, i-1), a.left, a.right), Multi_F{R}(b.start, min(b.stop, i-1), b.left, b.right),
            m,
            Multi_F{R}(max(a.start, j+1), a.stop, a.left, a.right), Multi_F{R}(max(b.start, j+1), b.stop, b.left, b.right)
        ]
    end

function simplify_core_prod(c::Number, v::Vector{<:IndexedOp{R}}) where R
    id = MakeIdentity(v[1])
    if c == 0
        return 0 * id
    end
    s = map(distribute(sumsubs.(v)...)) do p
        cp = c * prod(scalarcoef.(p))
        if cp == 0
            return 0 * id
        end
        r = reduce(vcat, prodsubs.(p))
        change = true
        while change
            change = false
            nr = IndexedOp{R}[]
            for right in r
                if right == id
                    continue
                elseif isempty(nr)
                    push!(nr, right)
                    continue
                else
                    left = nr[end]
                    t = orderprod(left, right)
                    if isempty(t)
                        push!(nr, right)
                    else
                        change = true
                        pop!(nr)
                        filter!(x->x ≠ id, t)
                        cp *= prod(scalarcoef.(t))
                        append!(nr, scalararg.(t))
                    end
                end
            end
            r = nr
        end
        return cp * ProdOp(r)
    end
    return simplify_sum(s)
end


################### removeMulti ###################

"""
    removeMulti(::Op)

transform Multi_F operators into their F equivalent
Multi_F(3, 5) => F(3)F(4)F(5)
"""
removeMulti(a::SumOp) = SumOp(removeMulti.(a.subs))
removeMulti(a::ProdOp) = ProdOp(removeMulti.(a.subs))
removeMulti(a::ScalarOp) = a.coef * removeMulti(a.arg)
removeMulti(a::AtIndex) = a
removeMulti(a::Multi_F{R}) where R = ProdOp([Multi_F{R}(i, i, a.left, a.right) for i in a.start:a.stop])