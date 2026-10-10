# The descent of an operator placed on sites onto a system: a first pass that develops nothing,
# gathering the factors of each product, then the list of its terms, each a coefficient times
# one matrix per site, the Jordan-Wigner strings laid, which the MPO, expect and the
# measurements start from, and apply for a factor it cannot apply as a dense block.

"""
    check_index(system, i, a)

refuse the site `i` of the operator `a` when `system` has no such site
"""
function check_index(system::System, i::Int, a)
    n = length(system)
    if i < 1 || i > n
        error("$a acts on site $i, which the system does not have: it has $n sites")
    end
    return nothing
end

"""
    check_indices(system, op)

check that every site the placed operator `op` acts on is a site of `system`, naming the
factor at fault
"""
check_indices(system::System, a::AtIndex) = foreach(i -> check_index(system, i, a), a.index)
check_indices(system::System, a::Union{SumOp, ProdOp}) =
    foreach(x -> check_indices(system, x), a.subs)
check_indices(system::System, a::ScalarOp) = check_indices(system, a.arg)
check_indices(::System, ::IdentityOp) = nothing

"""
    Term

a term of an operator on a system: `coef` times the tensor product of `blocks`, the matrix of
the term on each site it acts on, sorted by site, the identity on the others
"""
struct Term
    coef::Number
    blocks::Vector{Pair{Int, Matrix}}
end

"""
    a * b

the product of two terms, site by site, `nothing` when it vanishes on a site
"""
function (a::Term * b::Term)
    da, db = Dict(a.blocks), Dict(b.blocks)
    blocks = Pair{Int, Matrix}[]
    for i in sort!(collect(union(keys(da), keys(db))))
        x, y = get(da, i, nothing), get(db, i, nothing)
        if isnothing(x) || isnothing(y)
            push!(blocks, i => something(x, y))
            continue
        end
        m = x * y
        # a product of factors that vanishes, as Sp * Sp on a spin 1/2
        if norm(m) ≤ rounding_tol * norm(x) * norm(y)
            return nothing
        end
        push!(blocks, i => m)
    end
    return Term(a.coef * b.coef, blocks)
end

"""
    products(as, bs)

the terms of the product of two sums of terms, those that vanish left out
"""
products(as::Vector{Term}, bs::Vector{Term}) =
    Term[ t for t in (a * b for a in as for b in bs) if !isnothing(t) ]

"""
    scaled(c, terms)

the terms with their coefficients multiplied by `c`
"""
scaled(c::Number, ts::Vector{Term}) = [ Term(c * t.coef, t.blocks) for t in ts ]

"""
    dag(::Term)

the adjoint of a term, site by site
"""
dag(t::Term) = Term(conj(t.coef), [ i => collect(m') for (i, m) in t.blocks ])

"""
    one_site_terms(system, i, m)

the terms of the operator of matrix `m` on site `i` of `system`: on a fermionic site, its part
of each parity, the odd one with `F` on the fermionic sites before `i`
"""
function one_site_terms(system::System, i::Int, m::AbstractMatrix)
    f = matrix(F, system[i])
    p = parity_of(m, f)
    if p == 0
        return [ Term(1., [ i => m ]) ]
    elseif p == 1
        strings = [ j => f for j in 1:i-1 for f in (matrix(F, system[j]),) if f ≠ I ]
        return [ Term(1., [ strings ; i => m ]) ]
    end
    return [ one_site_terms(system, i, (m + f * m * f) / 2) ;
             one_site_terms(system, i, (m - f * m * f) / 2) ]
end

"""
    matrix_terms(system, op, index)

the terms of the operator `op` taken by its matrix on the sites `index` of `system`: on one
site, see `one_site_terms`; on several, the terms of `split_matrix`, each a product of
operators of one site in the order of `index`
"""
function matrix_terms(system::System, a::GenericOp{IsOperator}, index)
    sites = [ system[i] for i in index ]
    m = matrix(a, sites...)
    if length(index) == 1
        return one_site_terms(system, only(index), m)
    end
    terms = Term[]
    for (c, fs) in split_matrix(string(a(index...)), m, sites)
        ts = [ Term(c, Pair{Int, Matrix}[]) ]
        for (i, f) in zip(index, fs)
            if !isnothing(f)
                ts = products(ts, one_site_terms(system, i, f))
            end
        end
        append!(terms, ts)
    end
    return terms
end

"""
    op_terms(system, op, index)

the terms of the operator `op` on `system`, as it is written: sums and products distributed, a
generic operator placed on the sites `index`. A tensor product is the product of its factors
placed, `(A ⊗ B)(i, j) = A(i) * B(j)`, a named expression is its expression, and any other
operator is taken by its matrix, see `matrix_terms`.
"""
op_terms(system::System, a::SumOp{IsOperator}, index) =
    reduce(vcat, (op_terms(system, x, index) for x in a.subs))
op_terms(system::System, a::ProdOp{IsOperator}, index) =
    foldl(products, (op_terms(system, x, index) for x in a.subs))
op_terms(system::System, a::ScalarOp{IsOperator}, index) =
    scaled(a.coef, op_terms(system, a.arg, index))
op_terms(::System, ::IdentityOp{IsOperator}, _) = [ Term(1., Pair{Int, Matrix}[]) ]
op_terms(system::System, a::AtIndex{IsOperator}, _) = op_terms(system, a.op, a.index)
op_terms(system::System, a::TensorOp{IsOperator}, index) =
    foldl(products, (op_terms(system, o, index[p]) for (o, p) in zip(a.subs, factor_sites(a))))
op_terms(system::System, a::DagOp, index) = map(dag, op_terms(system, a.arg, index))
op_terms(system::System, a::Operator, index) =
    a.expr isa GenericOp ? op_terms(system, a.expr, index) : matrix_terms(system, a, index)
op_terms(system::System, a::Union{Proj, ExpOp{IsOperator}, IntPowOp{IsOperator},
                                  GenPowOp{IsOperator}, ModOp}, index) =
    matrix_terms(system, a, index)

"""
    reduced(terms)

the terms with each block equal to a multiple of the identity taken out, every other block
divided by its first element above rounding, that element going into the coefficient, a term
that vanishes on a site dropped, and the terms of the same blocks gathered, a real coefficient
kept real
"""
function reduced(terms::Vector{Term})
    out = Term[]
    seen = Dict{Vector{Pair{Int, Matrix}}, Int}()
    for t in terms
        coef = t.coef
        blocks = Pair{Int, Matrix}[]
        for (i, m) in t.blocks
            λ = tr(m) / size(m, 1)
            if iszero(m)
                coef = 0
                break
            elseif nearly(m, Matrix(λ * I, size(m)))
                coef *= λ
                continue
            end
            k = findfirst(x -> abs(x) > rounding_tol * norm(m), m)
            coef *= m[k]
            b = m / m[k]
            b[k] = 1
            if eltype(b) <: Complex && all(x -> iszero(imag(x)), b)
                b = real(b)
            end
            push!(blocks, i => no_signed_zero(b))
        end
        if coef == 0
            continue
        end
        j = get(seen, blocks, 0)
        if j == 0
            push!(out, Term(coef, blocks))
            seen[blocks] = length(out)
        else
            out[j] = Term(out[j].coef + coef, blocks)
        end
    end
    return [ Term(isreal(t.coef) ? real(t.coef) : t.coef, t.blocks) for t in out if t.coef ≠ 0 ]
end

"""
    check_flux(system, a, terms)

refuse the operator `a` on a charged `system` when a block of its terms has no definite flux,
or when its terms do not all change the charges by the same amount
"""
function check_flux(system::System, a, terms::Vector{Term})
    if !is_charged(system)
        return nothing
    end
    total = nothing
    for t in terms
        q = sum((charge_flux(m, a, system[i]) for (i, m) in t.blocks); init = QN())
        if isnothing(total)
            total = q
        elseif q ≠ total
            error("the terms of $a change the charges by $total and by $q, so it has no " *
                  "definite flux and cannot act on a system that conserves them")
        end
    end
    return nothing
end

"""
    system_matrix(system, terms)

the matrix of a sum of terms on the whole of `system`, the last site varying fastest
"""
function system_matrix(system::System, ts::Vector{Term})
    ids = [ Matrix{Float64}(I, dim(system[k]), dim(system[k])) for k in 1:length(system) ]
    d = prod(size(m, 1) for m in ids)
    return sum((t.coef * foldl(kron, [ get(Dict(t.blocks), k, ids[k]) for k in eachindex(ids) ])
                for t in ts); init = zeros(d, d))
end

"""
    sites_of(op)

the sites the placed operator `op` acts on
"""
sites_of(a::AtIndex) = collect(a.index)
sites_of(a::Union{SumOp, ProdOp}) = union((sites_of(x) for x in a.subs)...)
sites_of(a::ScalarOp) = sites_of(a.arg)
sites_of(::IdentityOp) = Int[]

"""
    parity_on(system, op)

the fermionic parity of the placed operator `op` on `system`, 0 or 1, read off the matrices of
its factors, `nothing` when it has none
"""
function parity_on(system::System, a::AtIndex)
    sites = [ system[i] for i in a.index ]
    fs = [ matrix(F, s) for s in sites ]
    if all(f -> f == I, fs)
        return 0
    end
    return parity_of(matrix(a.op, sites...), foldl(kron, fs))
end

function parity_on(system::System, a::ProdOp)
    ps = [ parity_on(system, x) for x in a.subs ]
    return any(isnothing, ps) ? nothing : mod(sum(ps), 2)
end

function parity_on(system::System, a::SumOp)
    ps = unique(parity_on(system, x) for x in a.subs)
    return length(ps) == 1 ? only(ps) : nothing
end

parity_on(system::System, a::ScalarOp) = parity_on(system, a.arg)
parity_on(::System, ::IdentityOp) = 0

"""
    matrix_on(system, op, sites)

the matrix of the operator `op`, placed on some of `sites`, on all of them in that order, the
Jordan-Wigner strings between them included
"""
function matrix_on(system::System, a::AtIndex, sites)
    if collect(a.index) == collect(sites)
        return matrix(a.op, [ system[i] for i in sites ]...)
    end
    sub = System([ system[i] for i in sites ])
    position = Dict(i => k for (k, i) in enumerate(sites))
    return system_matrix(sub, op_terms(sub, map_sites(i -> position[i], a), ()))
end

"""
    merged(system, a, b)

the product `a * b` of two placed operators, the sites of one among those of the other, as one
operator given by its matrix on the sites of the larger, or as a number when that matrix is a
multiple of the identity
"""
function merged(system::System, a::AtIndex, b::AtIndex)
    sites = length(a.index) ≥ length(b.index) ? a.index : b.index
    x, y = matrix_on(system, a, sites), matrix_on(system, b, sites)
    m = x * y
    λ = tr(m) / size(m, 1)
    if norm(m) ≤ rounding_tol * norm(x) * norm(y)
        return scale(0, IdentityOp(a))
    elseif nearly(m, Matrix(λ * I, size(m)))
        return scale(λ, IdentityOp(a))
    end
    return Operator{length(sites)}("($(a * b))", m)(sites...)
end

"""
    dependent((sites_a, parity_a), (sites_b, parity_b))

whether two factors of a product may not change places: they share a site, or one has no
definite parity and the other is odd or has none. Two factors on other sites otherwise commute,
up to a sign when both are odd.
"""
dependent((sa, pa), (sb, pb)) = !isdisjoint(sa, sb) || !(pa == 0 || pb == 0 || pa == pb == 1)

"""
    term_count(system, op)

the number of terms of the placed operator `op` on `system`
"""
term_count(system::System, a::IndexedOp{IsOperator}) = length(reduced(op_terms(system, a, ())))

"""
    sort_factors!(fs, info)

sort the factors `fs` of a product, and their `(sites, parity)` in `info`, by their sites as
far as they (anti)commute, see `dependent`, a factor going back past a neighbour of larger
sites, and return the sign that brings, -1 for each pair of odd factors that change places
"""
function sort_factors!(fs, info)
    sign = 1
    key(k) = sort(first(info[k]))
    for k in 2:length(fs)
        for j in k:-1:2
            if !(key(j - 1) > key(j)) || dependent(info[j - 1], info[j])
                break
            end
            if info[j - 1][2] == 1 && info[j][2] == 1
                sign = -sign
            end
            fs[j - 1], fs[j] = fs[j], fs[j - 1]
            info[j - 1], info[j] = info[j], info[j - 1]
        end
    end
    return sign
end

"""
    merge_one!(system, fs, info, refused)

merge into an earlier factor of the product `fs`, see `merged`, one later factor sharing one of
its sites, the sites of one among those of the other, which goes back next to it past factors
it (anti)commutes with, provided the merge gives no more terms than the two factors, see
`term_count`; a pair whose merge would give more goes into `refused`. Return the number the
product gets multiplied by, or `nothing` when no pair can be merged.
"""
function merge_one!(system::System, fs, info, refused)
    for j in 2:length(fs)
        for i in j-1:-1:1
            # factor j goes back past those between, which it must (anti)commute with
            if i < j - 1 && dependent(info[i + 1], info[j])
                break
            end
            a, b = fs[i], fs[j]
            sa, sb = first(info[i]), first(info[j])
            if !(a isa AtIndex && b isa AtIndex) || isdisjoint(sa, sb) ||
               !(issubset(sa, sb) || issubset(sb, sa)) || (a, b) in refused
                continue
            end
            m = merged(system, a, b)
            if m isa AtIndex && term_count(system, m) > term_count(system, a) * term_count(system, b)
                push!(refused, (a, b))
                continue
            end
            sign = info[j][2] == 1 ? (-1)^count(k -> info[k][2] == 1, i+1:j-1) : 1
            deleteat!(fs, j)
            if m isa AtIndex
                fs[i] = m
                return sign
            end
            deleteat!(fs, i)
            return sign * scalarcoef(m)
        end
    end
    return nothing
end

"""
    gathered(system, op)

the placed operator `op` rewritten without developing anything, into a normal form: in each
product, the factors are sorted by their sites as far as they (anti)commute, see
`sort_factors!`, and two factors sharing a site are merged where they can be, see
`merge_one!`, until nothing moves. Sums are gathered term by term, and an operator in normal
form is left as it is.
"""
gathered(system::System, a::SumOp{IsOperator, Indexed}) =
    SumOp([ gathered(system, x) for x in a.subs ])
gathered(system::System, a::ScalarOp{IsOperator, Indexed}) =
    scale(a.coef, gathered(system, a.arg))
gathered(::System, a::Union{AtIndex{IsOperator}, IdentityOp{IsOperator, Indexed}}) = a

function gathered(system::System, a::ProdOp{IsOperator, Indexed})
    p = ProdOp(IndexedOp{IsOperator}[ gathered(system, x) for x in a.subs ])
    coef = scalarcoef(p)
    fs = IndexedOp{IsOperator}[ x for x in prodsubs(p) if !(x isa IdentityOp) ]
    refused = Set{Tuple{IndexedOp{IsOperator}, IndexedOp{IsOperator}}}()
    while coef ≠ 0
        info = [ (sites_of(x), parity_on(system, x)) for x in fs ]
        coef *= sort_factors!(fs, info)
        c = merge_one!(system, fs, info, refused)
        if isnothing(c)
            break
        end
        coef *= c
    end
    return scale(coef, ProdOp(fs))
end

"""
    lower(system, op)

the terms of the operator `op`, placed on sites, on `system`, see `Term`, once its products are
gathered, see `gathered`: the Jordan-Wigner strings are laid on the fermionic sites, and the
operator is refused on a charged system when it has no definite flux
"""
function lower(system::System, a::IndexedOp{IsOperator})
    # before any rewriting, so that a message names the operator as it was written
    check_indices(system, a)
    terms = reduced(op_terms(system, gathered(system, a), ()))
    check_flux(system, a, terms)
    return terms
end
