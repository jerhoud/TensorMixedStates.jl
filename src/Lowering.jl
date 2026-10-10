# The descent of an operator placed on sites onto a system: the list of its terms, each a
# coefficient times one matrix per site, the Jordan-Wigner strings laid, which the MPO, expect
# and the measurements start from, and apply for a factor it cannot apply as a dense block.

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
    lower(system, op)

the terms of the operator `op`, placed on sites, on `system`, see `Term`: the Jordan-Wigner
strings are laid on the fermionic sites, and the operator is refused on a charged system when
it has no definite flux
"""
function lower(system::System, a::IndexedOp{IsOperator})
    # before any rewriting, so that a message names the operator as it was written
    check_indices(system, a)
    terms = reduced(op_terms(system, a, ()))
    check_flux(system, a, terms)
    return terms
end
