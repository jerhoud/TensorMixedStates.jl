# simplify, which brings an operator to the normal form the MPO construction expects, a sum of
# products of one site operators ordered by site with their Jordan-Wigner strings, and
# removeMulti, which spells those strings out site by site.

export simplify

"""
    simplify(op)

the operator `op` in a normal form, as the functions building an MPO compute it.

Operators of several sites defined by an expression are replaced by it, products of sums are
expanded, coefficients and like factors are gathered, and fermionic operators of one site get
their Jordan-Wigner strings. A placed operator thus becomes a sum of products of one site
operators ordered by site, which is what `PreMPO` expects. An operator of several sites that
cannot be developed, a function of one, such as its exponential, or one defined by a matrix
and created without its sites, is kept whole, and keeps its place with respect to the
Jordan-Wigner strings. A collection of operators is simplified element by element.

# Examples

    simplify(X(1) * Z(2) * X(1))
    simplify(C(3) * dag(C)(5))
"""
function simplify end


################### Sums and products in normal form ###################

"""
    pow_base(a)
    pow_expo(a)

the base and the exponent of `a` as a power: `a` and `1` when it is not one.
"""
pow_base(a::Op) = a
pow_base(a::Union{IntPowOp, GenPowOp}) = a.arg
pow_expo(a::Op) = 1
pow_expo(a::Union{IntPowOp, GenPowOp}) = a.expo

"""
    distribute(terms...)

every product of one term taken from each of the vectors `terms`, as a vector of factors: the
expansion of a product of sums.
"""
distribute(terms::Vector...) =
    vec([ collect(reverse(p)) for p in Iterators.product(reverse(terms)...) ])

"""
    split_string(b, a, string_first)

the string `b` split around the operator `a`, placed on a site `i` inside it: the string up to
`i - 1`, the `F` of site `i` and `a`, and the rest of the string. The `F` is put before `a`
when `string_first` is true, the string having come first in the product, and after it
otherwise.
"""
function split_string(b::Multi_F{R}, a::AtIndex{R, 1}, string_first::Bool) where R
    i = only(a.index)
    piece(start, stop) = Multi_F{R}(start, stop, b.left, b.right)
    return string_first ? [piece(b.start, i - 1), piece(i, i), a, piece(i + 1, b.stop)] :
                          [piece(b.start, i - 1), a, piece(i, i), piece(i + 1, b.stop)]
end

"""
    simplify_sum(v)

the sum of the simplified operators `v`, flattened and sorted, equal terms collected and zero
ones removed: `X + (Y + Z)` gives `X + Y + Z` and `X + Y + X` gives `2X + Y`. In an indexed
sum, terms on the same sites are gathered, `X(1) + Y(1)` giving `(X + Y)(1)`, and so are
products differing only in their last factor, of one site, `P * X(2) + P * Y(2)` giving
`P * (X + Y)(2)`. The odd terms of a sum on one site, whose Jordan-Wigner string is a factor
of each, thus share it, rather than making a sum a gate refuses.
"""
simplify_sum(v::Vector) = simplify_core_sum(reduce(vcat, sumsubs.(v)))

"""
    gathered(c, o, nc, no)

the sum of the indexed terms `c * o` and `nc * no` as a single term, when they are on the same
sites or are products of the same factors but their last ones, placed on the same sites, and
`nothing` otherwise
"""
function gathered(c, o, nc, no)
    if o isa AtIndex && no isa AtIndex && o.index == no.index
        return simplify_sum([c * o.op, nc * no.op])(o.index...)
    elseif o isa ProdOp && no isa ProdOp && length(o.subs) == length(no.subs) &&
           o.subs[end] isa AtIndex && no.subs[end] isa AtIndex &&
           o.subs[end].index == no.subs[end].index && o.subs[1:end-1] == no.subs[1:end-1]
        l, nl = o.subs[end], no.subs[end]
        return simplify_prod([o.subs[1:end-1]..., simplify_sum([c * l.op, nc * nl.op])(l.index...)])
    end
    return nothing
end

"""
    simplify_core_sum(v)

the sum of the flattened terms `v`, as `simplify_sum` describes.
"""
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
            c, o = nc, no
        else
            m = T == Indexed ? gathered(c, o, nc, no) : nothing
            if isnothing(m)
                push!(r, c * o)
                c, o = nc, no
            else
                c, o = scalarcoef(m), scalararg(m)
            end
        end
    end
    if c ≠ 0
        push!(r, c * o)
    end
    return SumOp(r)
end

"""
    simplify_prod(v)

the product of the simplified operators `v`, flattened, `X * (Y * Z)` giving `X * Y * Z`, with
its coefficients gathered, and put in normal form by `simplify_core_prod`.
"""
function simplify_prod(v::Vector)
    c = prod(scalarcoef.(v))
    subs = reduce(vcat, prodsubs.(v))
    if c == 0
        return 0 * IdentityOp(subs[1])
    end
    return simplify_core_prod(c, subs)
end

"""
    simplify_core_prod(c, v)

`c` times the product of the flattened factors `v`, in normal form.

- Generic pure: the `F` are carried to the right, with the sign the parity of each factor
  crossed gives, `F * JW(C)` giving `-JW(C) * F`, see `jw_parity`, and are laid down before a
  factor of no definite parity. Neighbouring powers of one base merge, `X * X` giving `X^2`,
  an involution squared giving the identity.
- Generic mixed: the `Left` factors are gathered into one and the `Right` ones into another,
  `Left(X) * Right(Y) * Left(Z)` giving `Left(X * Z) * Right(Y)`. A product with other
  factors, such as sums, is kept as it is.
- Indexed: the sums are expanded with `distribute` and each product is sorted by site with
  `orderprod`, factors on the same sites merging and the strings `Multi_F` being glued, split
  or cancelled.
"""
function simplify_core_prod(c::Number, v::Vector{<:GenericOp{Pure, N}}) where N
    id = IdentityOp(v[1])
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

function simplify_core_prod(c::Number, v::Vector{<:GenericOp{Mixed, N}}) where N
    id = IdentityOp(v[1])
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
    for (S, args) in ((Left, larg), (Right, rarg))
        if !isempty(args)
            m = sided(S, simplify_prod(args))
            if m ≠ id
                push!(r, m)
            end
        end
    end
    return c * ProdOp(r)
end

function simplify_core_prod(c::Number, v::Vector{<:IndexedOp{R}}) where R
    s = map(distribute(sumsubs.(v)...)) do p
        cp = c * prod(scalarcoef.(p))
        r = filter(x -> !(x isa IdentityOp), reduce(vcat, prodsubs.(p)))
        change = true
        while change
            change = false
            nr = IndexedOp{R}[]
            for right in r
                t = isempty(nr) ? [] : orderprod(nr[end], right)
                if isempty(t)
                    push!(nr, right)
                else
                    change = true
                    pop!(nr)
                    # the coefficients first: a merge giving -Id left its identity in the
                    # product, the filter seeing a scalar times it
                    cp *= prod(scalarcoef.(t))
                    append!(nr, filter(x -> !(x isa IdentityOp), scalararg.(t)))
                end
            end
            r = nr
        end
        return cp * ProdOp(r)
    end
    return simplify_sum(s)
end

"""
    cannot_multiply(a)

refuse a product of the com `a` with another operator: the channels of a com are what is left
of its terms, which a product would have to multiply one by one.
"""
cannot_multiply(a::ComOp) =
    error("$a is compacted and cannot be multiplied: take the product first and compact it")

"""
    orderprod(a, b)

what replaces the product `a * b` of two placed operators or strings `Multi_F`, one step of the
sort by site of an indexed product: an empty list when the pair stays as it is, the pair
swapped when it is out of order, a single factor for two operators on the same sites or two
adjacent strings, and a string split around an operator of one site inside it. The sort
takes, for instance, `X(1) * Z(2) * Y(1)` to `(X * Y)(1) * Z(2)`, and `C(3) * C(5)`, that is
`Multi_F(1, 2) * JW(C)(3) * Multi_F(1, 4) * JW(C)(5)`, to `(JW(C) * F)(3) * F(4) * JW(C)(5)`.
"""
orderprod(a::AtIndex, b::AtIndex) =
    if a.index == b.index
        [ simplify_prod([a.op, b.op])(a.index...) ]
    elseif min(a.index...) > max(b.index...)
        [b, a]
    else
        []
    end

function orderprod(a::AtIndex{R, 1}, b::Multi_F{R}) where R
    i = only(a.index)
    return i < b.start ? [] : i > b.stop ? [b, a] : split_string(b, a, false)
end

function orderprod(b::Multi_F{R}, a::AtIndex{R, 1}) where R
    i = only(a.index)
    return i < b.start ? [a, b] : i > b.stop ? [] : split_string(b, a, true)
end

# a factor of several sites is never moved past a string, which commutes with it only if it
# holds all of its sites or none of them, and the factor is even in the first case. What is
# left of one after simplify is refused wherever it would be placed, see check_one_site, and a
# gate keeps its factors of several sites out of simplify, see prepare_gate
orderprod(::AtIndex{R}, ::Multi_F{R}) where R = []
orderprod(::Multi_F{R}, ::AtIndex{R}) where R = []

function orderprod(a::Multi_F{R}, b::Multi_F{R}) where R
    piece(s, start, stop) = Multi_F{R}(start, stop, s.left, s.right)
    if a.left == b.left && a.right == b.right && (a.stop == b.start-1 || b.stop == a.start-1)
        return [ piece(a, min(a.start, b.start), max(a.stop, b.stop)) ]
    elseif a.stop < b.start
        return []
    elseif b.stop < a.start
        return [b, a]
    end
    i = max(a.start, b.start)
    j = min(a.stop, b.stop)
    m = a.left == b.left && a.right == b.right ? IdentityOp(a) :
                                                 Multi_F{R}(i, j, a.left ⊻ b.left, a.right ⊻ b.right)
    return [ piece(a, a.start, min(a.stop, i-1)), piece(b, b.start, min(b.stop, i-1)), m,
             piece(a, max(a.start, j+1), a.stop), piece(b, max(b.start, j+1), b.stop) ]
end

orderprod(a::ComOp, ::IndexedOp) = cannot_multiply(a)
orderprod(::Union{AtIndex, Multi_F}, b::ComOp) = cannot_multiply(b)


################### Adjoint and exponential ###################

"""
    simplify_dag(a)

the adjoint of the simplified operator `a`, carried as deep into it as it goes, down to the
operators that are their own adjoint by their type: the identity, `F`, projectors,
involutions and self adjoint operators. A non integer power, and a tensor product with a
factor of odd or undefined parity, are kept under `dag`.
"""
simplify_dag(a::DagOp) = a.arg
simplify_dag(a::ScalarOp) = conj(a.coef) * simplify_dag(a.arg)
simplify_dag(a::Operator) =
    if a.type == involution_op || a.type == selfadjoint_op
        a
    else
        dag(a)
    end
# the adjoint of a product repeated is the adjoint repeated, but a non integer power goes
# through a logarithm, whose branch cut the adjoint does not respect
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


simplify_dag(a::AtIndex) = simplify_dag(a.op)(a.index...)
simplify_dag(a::Multi_F) = a
# piece by piece: the factors of a path sit on distinct sites, their strings in place
simplify_dag(a::ComOp{Pure}) = map_pieces(simplify_dag, Pure, a)

"""
    simplify_exp(a)

the exponential of the simplified pure operator `a`: `cosh(c) * Id + sinh(c) * X` when `a` is
`c * X` with `X` an involution, `exp(a)` otherwise.
"""
function simplify_exp(a::GenericOp{Pure, N}) where N
    c = scalarcoef(a)
    s = prodsubs(a)
    if length(s) == 1 && is_involution(s[1])
        simplify_sum([cosh(c) * IdentityOp(s[1]), sinh(c) * s[1]])
    else
        exp(a)
    end
end


################### Placing an operator on its sites ###################

"""
    simplify_ind(op, index...)

the generic operator `op` placed on the sites `index` and simplified, the placement carried as
deep into the expression as it goes: tensor products are split over their sites, fermionic
operators of one site get their Jordan-Wigner string, and operators of several sites defined
by an expression are replaced by it.
"""
simplify_ind(a::ScalarOp, index...) = a.coef * simplify_ind(a.arg, index...)
simplify_ind(a::IdentityOp, index...) = a(index...)
simplify_ind(a::Union{JW_F, Proj, JW, SetState}, index) = a(index)
simplify_ind(a::ExpOp, index...) = place_function(a, index...)
simplify_ind(a::ModOp, index...) = place_function(a, index...)

"""
    place_function(a, index...)

the function `a` of an operator, an exponential, a `mod` or a non integer power, placed on the
sites `index`. On a single site other than the first, when the argument is not even, it is
the sum of its part commuting with `F`, placed bare, and its part anticommuting with `F`,
which takes the Jordan-Wigner string as `C` does. Each part is an operator of its own, the odd
one fermionic, so that simplifying the result again leaves it unchanged. Otherwise the
function is placed whole.
"""
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
simplify_ind(a::Left, index...) = sided(Left, simplify_ind(a.arg, index...))
simplify_ind(a::Right, index...) = sided(Right, simplify_ind(a.arg, index...))

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
        # simplified first, as an operator placed as it is written would be: a renamed
        # exp(c * Swap) has to be developed as exp(c * Swap) itself is
        simplify_ind(simplify(a.expr), index...)
    else
        a(index...)
    end

simplify_ind(a::SumOp, index...) = simplify_sum(map(x->simplify_ind(x, index...), a.subs))
simplify_ind(a::ProdOp, index...) = simplify_prod(map(x->simplify_ind(x, index...), a.subs))
simplify_ind(a::TensorOp, index...) =
    simplify_prod([ simplify_ind(o, index[p]...) for (o, p) in zip(a.subs, factor_sites(a)) ])


################### simplify ###################

simplify(a) = map(simplify, a)

simplify(a::ScalarOp) = a.coef * simplify(a.arg)
simplify(a::ProdOp) = simplify_prod(map(simplify, a.subs))
simplify(a::SumOp) = simplify_sum(map(simplify, a.subs))
simplify(a::TensorOp{N}) where N = TensorOp{N}(simplify.(a.subs))

# a com is built by `compact` from simplified terms
simplify(a::Union{IdentityOp, JW_F, Proj, JW, Operator, Multi_F, SetState, ComOp}) = a

simplify(a::Union{IntPowOp, GenPowOp}) = power(simplify(a.arg), a.expo)
simplify(a::ExpOp) = simplify_exp(simplify(a.arg))
simplify(a::DagOp) = simplify_dag(simplify(a.arg))
simplify(a::ModOp) = ModOp(simplify(a.arg), a.modulus)

function simplify(a::Dissipator)
    sarg = simplify(a.arg)
    darg = simplify_dag(sarg)
    daga = simplify_prod([darg, sarg])
    simplify_sum([simplify_prod([sided(Left, sarg), sided(Right, sarg)]),
                  -0.5 * sided(Left, daga), -0.5 * sided(Right, daga)])
end


simplify(a::Left) = sided(Left, simplify(a.arg))
simplify(a::Right) = sided(Right, simplify(a.arg))

function simplify(a::Gate)
    sarg = simplify(a.arg)
    simplify_prod([sided(Left, sarg), sided(Right, sarg)])
end

function simplify(a::Evolver)
    sarg = simplify(a.arg)
    simplify_sum([sided(Left, sarg), sided(Right, sarg)])
end

simplify(a::AtIndex) =
    simplify_ind(simplify(a.op), a.index...)


################### removeMulti ###################

"""
    removeMulti(op)

the operator with each string `Multi_F` spelled out as one factor per site, `Multi_F(3, 5)`
giving `F(3) * F(4) * F(5)`. A collection is processed element by element.
"""
removeMulti(a::SumOp) = SumOp(removeMulti.(a.subs))
removeMulti(a::ProdOp) = ProdOp(removeMulti.(a.subs))
removeMulti(a::ScalarOp) = a.coef * removeMulti(a.arg)
removeMulti(a::Union{AtIndex, IdentityOp, ComOp}) = a
removeMulti(a::Multi_F{R}) where R = ProdOp([Multi_F{R}(i, i, a.left, a.right) for i in a.start:a.stop])
removeMulti(a) = map(removeMulti, a)
