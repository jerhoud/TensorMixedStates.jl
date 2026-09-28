export simplify

"""
    simplify(op::Op)

simplifies an operator, this is used internally by functions creating MPOs.
Multi site operators defined by an expression are replaced by that expression, so that
the result is a sum of products of one site operators, which is what `PreMPO` expects.
"""
function simplify end

# Simplifications for collections of Operators

simplify(a) = map(simplify, a)


# Simplification principles
# One pass
# expand definitions of Gate and Evolver
# expand as go we in, simplify as we go out
# the final result should be a sum of products
# of one site (possibly complex) indexed operators with different indices in ascending order




# Simplifications for both Generic and Indexed operators

simplify(a::ScalarOp) = a.coef * simplify(a.arg)
simplify(a::ProdOp) = simplify_prod(map(simplify, a.subs))
simplify(a::SumOp) = simplify_sum(map(simplify, a.subs))
simplify(a::TensorOp{N}) where N = TensorOp{N}(simplify.(a.subs))


# Simplification of Generic Operators

simplify(a::Union{IdentityOp, JW_F, Proj, JW, Operator, Multi_F, SetState}) = a

simplify(a::IntPowOp) = power(simplify(a.arg), a.expo)
simplify(a::GenPowOp) = power(simplify(a.arg), a.expo)
simplify(a::ExpOp) = simplify_exp(simplify(a.arg))
simplify(a::DagOp) = simplify_dag(simplify(a.arg))
simplify(a::ModOp) = ModOp(simplify(a.arg), a.modulus)

function simplify(a::Dissipator)
    sarg = simplify(a.arg)
    darg = simplify_dag(sarg)
    daga = simplify_prod([darg, sarg])
    simplify_sum([simplify_prod([simplify_l(sarg), simplify_r(sarg)]),
                  -0.5 * simplify_l(daga), -0.5 * simplify_r(daga)])
end


simplify(a::Left) = simplify_l(simplify(a.arg))
simplify(a::Right) = simplify_r(simplify(a.arg))


# Simplifications of Indexed Operators

function simplify(a::Gate)
    sarg = simplify(a.arg)
    simplify_prod([simplify_l(sarg), simplify_r(sarg)])
end

function simplify(a::Evolver)
    sarg = simplify(a.arg)
    simplify_sum([simplify_l(sarg), simplify_r(sarg)])
end

simplify(a::AtIndex) =
    simplify_ind(simplify(a.op), a.index...)

reindex(op::GenericOp, i::Int...) = op(i...)

# Simplification with index
# transmit indexation as deep as possible
# to develop tensors, transform fermionic operators with JW, replace multi site operators by their definition


simplify_ind(a::ScalarOp, index...) = a.coef * simplify_ind(a.arg, index...)
# placed, the identity has no site, AtIndex giving the one of the whole system
simplify_ind(a::IdentityOp, index...) = a(index...)
simplify_ind(a::Union{JW_F, Proj, JW, SetState}, index) = a(index)
simplify_ind(a::ExpOp, index...) = place_function(a, index...)
simplify_ind(a::ModOp, index...) = place_function(a, index...)

# a function of an operator of one site that is not even is kept whole, and placed after other
# sites it is its part commuting with F, placed bare, plus its part anticommuting with F, which
# takes the string as C does: placed whole, it had no string at all. Each part is an operator
# of its own, the odd one fermionic, so that simplifying the result again leaves it as it is
# rather than cutting the function inside the even part once more
place_function(a, index...) =
    if length(index) == 1 && only(index) > 1 && jw_parity(a.arg) ≠ 0
        i = only(index)
        even = Operator{1}("even($a)", 0.5 * (a + F * a * F), plain_op)
        odd = Operator{1}("odd($a)", 0.5 * (a - F * a * F), fermionic_op)
        simplify_sum([simplify_ind(even, i), simplify_ind(odd, i)])
    else
        a(index...)
    end
# an integer power is its product, placed factor by factor with their strings, and any other is
# a function of the operator, placed as exp is: on one site through its matrix, split in the
# parts that commute and anticommute with F, and whole on several
simplify_ind(a::IntPowOp, index...) = simplify_prod(fill(simplify_ind(a.arg, index...), a.expo))

function simplify_ind(a::GenPowOp{Pure}, index...)
    g = simplify(a)
    p = scalararg(g)
    if !(p isa GenPowOp)
        return simplify_ind(g, index...)
    end
    return scalarcoef(g) * place_function(p, index...)
end

simplify_ind(a::GenPowOp{Mixed}, index...) = simplify(a)(index...)
simplify_ind(a::DagOp, index...) = simplify_dag(simplify_ind(a.arg, index...))
simplify_ind(a::Left, index...) = simplify_l(simplify_ind(a.arg, index...))
simplify_ind(a::Right, index...) = simplify_r(simplify_ind(a.arg, index...))

simplify_ind(a::Operator{1}, index) =
    if a.type == fermionic_op
        if index > 1
            Multi_F{Pure}(1, index-1, false, false) * JW(a)(index)
        else
            JW(a)(index)
        end
    else
        a(index)
    end

# a multi site operator has to be replaced by its definition, PreMPO only knows how to
# place one site factors. One site operators are handled by the method above and keep
# their name, their definition is read from the site when the tensor is needed.
simplify_ind(a::Operator, index...) =
    if a.expr isa Op
        simplify_ind(a.expr, index...)
    else
        a(index...)
    end

simplify_ind(a::SumOp, index...) = simplify_sum(map(x->simplify_ind(x, index...), a.subs))
simplify_ind(a::ProdOp, index...) = simplify_prod(map(x->simplify_ind(x, index...), a.subs))
simplify_ind(a::TensorOp, index...) = simplify_prod(tensor_apply(simplify_ind, a, index...))


# simplify exp : exp(0) => Id and exp(3X) => cosh(3)Id + sinh(3)X
function simplify_exp(a::GenericOp{Pure, N}) where N
    c = scalarcoef(a)
    s = prodsubs(a)
    if length(s) == 1 && is_involution(s[1])
        simplify_sum([cosh(c) * IdentityOp(s[1]), sinh(c) * s[1]])
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
# the adjoint of a product repeated is the adjoint repeated, but a non integer power goes
# through a logarithm, whose branch cut the adjoint does not respect: dag(sqrt(X)) came out
# as sqrt(X)
simplify_dag(a::IntPowOp) = power(simplify_dag(a.arg), a.expo)
simplify_dag(a::GenPowOp) = DagOp(a)
simplify_dag(a::ExpOp) = ExpOp(simplify_dag(a.arg))
simplify_dag(a::ModOp) = ModOp(-simplify_dag(a.arg), a.modulus)
# a projector is on a pure state, `matrix` refusing a mixed one, and so self adjoint
simplify_dag(a::Union{IdentityOp, JW_F, Proj}) = a
simplify_dag(a::JW) = dag(a)


simplify_dag(a::ProdOp) = simplify_prod(reverse(simplify_dag.(a.subs)))
simplify_dag(a::SumOp) = simplify_sum(simplify_dag.(a.subs))
# the adjoint of a tensor product reverses its factors once placed, which shows only when two
# of them anticommute, its sites being distinct: with a factor of odd or undefined parity it
# waits for the sites, where the product of the placed factors takes the sign
simplify_dag(a::TensorOp{N}) where N =
    if all(o -> jw_parity(o) == 0, a.subs)
        TensorOp{N}(simplify_dag.(a.subs))
    else
        DagOp(a)
    end


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
simplify_l(::IdentityOp{Pure, Indexed, 1}) = IdentityOp{Mixed, Indexed, 1}()


# Right simplification
# Generic => get the scalar factor out and simplifies Right(Id) => Left(Id)
# Indexed => go as deep as possible

simplify_r(a::GenericOp{Pure}) = Right(a)
simplify_r(a::ScalarOp{Pure}) = conj(a.coef) * simplify_r(a.arg) 

simplify_r(a::ProdOp{Pure, Indexed}) = ProdOp(simplify_r.(a.subs)) 
simplify_r(a::SumOp{Pure, Indexed}) = SumOp(simplify_r.(a.subs))
simplify_r(a::AtIndex{Pure}) = reindex(simplify_r(a.op), a.index...)
simplify_r(a::Multi_F{Pure}) = Multi_F{Mixed}(a.start, a.stop, false, true)
simplify_r(::IdentityOp{Pure, Indexed, 1}) = IdentityOp{Mixed, Indexed, 1}()


# sum simplification
# flatten out inner sums, order terms, collect identical terms and remove nuls
# X + (Y + Z) => X + Y + Z, X + Y + X => 2X + Y, X - X => 0
# in Indexed sums gather terms with same indices : X(1) + Y(1) => (X+Y)(1)
# and products that differ only by their last factor, on one site: P*X(2) + P*Y(2) => P*(X+Y)(2).
# The Jordan-Wigner string of an odd term being a factor of its own, the terms of an odd sum of
# one site share it: kept apart, they were a sum a gate refuses and an MPO carries one by one

simplify_sum(v::Vector) = simplify_core_sum(reduce(vcat, sumsubs.(v)))

same_but_last(a, b) =
    a isa ProdOp && b isa ProdOp && length(a.subs) == length(b.subs) &&
    a.subs[end] isa AtIndex && b.subs[end] isa AtIndex &&
    a.subs[end].index == b.subs[end].index && a.subs[1:end-1] == b.subs[1:end-1]

function simplify_core_sum(v::Vector{<:Op{R, T, N}}) where {R, T, N}
    subs = sort(v; by=scalararg)
    r = Op{R, T, N}[]
    c = 0
    o = IdentityOp{R, T, N}()
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
        elseif T == Indexed && same_but_last(o, no)
            l, nl = o.subs[end], no.subs[end]
            o = simplify_prod([o.subs[1:end-1]..., reindex(simplify_sum([c * l.op, nc * nl.op]), l.index...)])
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
# carry F to the right end with the sign the parity of each factor crossed gives, see jw_parity
# F*X = X*F, F*JW(C) => -JW(C)*F, F*F = Id

pow_base(a::Op) = a
pow_base(a::Union{IntPowOp, GenPowOp}) = a.arg
pow_expo(a::Op) = 1
pow_expo(a::Union{IntPowOp, GenPowOp}) = a.expo

function simplify_core_prod(c::Number, v::Vector{<:GenericOp{Pure, N}}) where N
    id = IdentityOp(v[1])
    if c == 0
        return 0 * id
    end
    w = GenericOp{Pure, N}[]
    f = false
    for x in v
        if x isa JW_F
            f = !f
            continue
        end
        p = jw_parity(x)
        if f && isnothing(p)
            # a factor of no definite parity stops the F, which are laid down just before it
            push!(w, F)
            f = false
        elseif f && p == 1
            c = -c
        end
        push!(w, x)
    end
    if f
        push!(w, F)
    end
    # neighbours of one base merge, A^p * A^q being A^(p + q), and what a merge gives has its
    # coefficient taken out of the product: (-X)^0.5 * (-X)^0.5 is -X, whose -1 left as a
    # factor kept the product from its normal form. A merge that gives F, as F^0.5 * F^0.5
    # does, has to go through the crossing of the F again
    r = GenericOp{Pure, N}[]
    again = false
    for x in w
        if !isempty(r) && pow_base(r[end]) == pow_base(x)
            y = pop!(r)
            m = power(pow_base(x), pow_expo(y) + pow_expo(x))
            c *= scalarcoef(m)
            m = scalararg(m)
            again = again || m isa JW_F
        else
            m = x
        end
        if m ≠ id
            push!(r, m)
        end
    end
    if again
        return simplify_core_prod(c, r)
    end
    return c * ProdOp(r)
end

# Generic Mixed products
# (that is only Left and Right factors)
# Gather Left and Right together
# Left(X)*Right(Y)*Left(Z) => Left(X*Z)*Right(Y)

function simplify_core_prod(c::Number, v::Vector{<:GenericOp{Mixed, N}}) where N
    id = IdentityOp(v[1])
    if c == 0
        return 0 * id
    end
    # the identity is neither a Left nor a Right, and would keep the others from gathering
    v = filter(x -> !(x isa IdentityOp), v)
    if isempty(v)
        return c * id
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
function distribute(a::Vector{<:Vector}, b::Vector, c::Vector...)
    r = Vector{Vector}(undef, length(a) * length(b))
    n = 1
    for i in a
        for j in b
            r[n] = vcat(i, [j])
            n += 1
        end
    end
    return distribute(r, c...)
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
            m = IdentityOp(a)
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
    id = IdentityOp(v[1])
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
                if right isa IdentityOp
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
                        filter!(x -> !(x isa IdentityOp), t)
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
removeMulti(a::Union{AtIndex, IdentityOp}) = a
removeMulti(a::Multi_F{R}) where R = ProdOp([Multi_F{R}(i, i, a.left, a.right) for i in a.start:a.stop])
removeMulti(a) = map(removeMulti, a)