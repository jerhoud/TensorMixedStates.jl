# The descent of an operator placed on sites onto a system: a first pass that develops nothing,
# gathering the factors of each product, then the list of its terms, each a coefficient times
# one matrix per site, the Jordan-Wigner strings kept apart, which the MPO, expect and the
# measurements start from, and apply for a factor it cannot apply as a dense block. A channel or
# a generator descends the same way onto the density matrices of the sites, each matrix that of
# a superoperator of one site.

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
check_indices(system::System, a::Union{Gate{Indexed}, Dissipator{Indexed}, Evolver{Indexed}}) =
    check_indices(system, a.arg)

"""
    Term

a term of an operator on a system: `coef` times the tensor product of `blocks`, the matrix of
the term on each site it acts on, sorted by site, the identity on the others, each followed by
the Jordan-Wigner strings of the term there, which the descent keeps apart, laid as blocks by
what needs them, see `laid_strings`. The term of a superoperator has the matrices of
superoperators of one site, see `matrix`, and a string on the kets, ``\\rho \\mapsto F\\rho``,
and one on the bras, ``\\rho \\mapsto \\rho F``.

- `parities`: the parity of each block, that of the kets plus twice that of the bras for a
  superoperator, which gives the sign of a string crossing it
- `kets`, `bras`: the strings, as the sorted sites where their parity changes: `F` on the sites
  from the first to the second excluded, from the third to the fourth, and so on. Each of these
  sites is site 1 or a site where the term has a block, which `compacted` relies on.
"""
struct Term
    coef::Number
    blocks::Vector{Pair{Int, Matrix}}
    parities::Vector{Int}
    kets::Vector{Int}
    bras::Vector{Int}
end

"""
    on_string(s, k)

whether site `k` lies on the string `s`, see `Term`
"""
on_string(s::Vector{Int}, k::Int) = isodd(searchsortedlast(s, k))

"""
    crossing_sign(kets, bras, blocks, parities)

the sign the strings `kets` and `bras` bring crossing `blocks`, pairs of a site and a matrix, of
parities `parities`: -1 for each block odd on the side of a string over its site, `F A = -A F`
"""
crossing_sign(kets::Vector{Int}, bras::Vector{Int}, blocks, parities) =
    (-1)^count((isodd(p) && on_string(kets, k)) ⊻ (isodd(p >> 1) && on_string(bras, k))
               for ((k, _), p) in zip(blocks, parities))

"""
    xor_strings(s, t)

the string of `F` on the sites where one of the strings `s` and `t` lies and not the other: the
sites where its parity changes are those where one of theirs does
"""
function xor_strings(s::Vector{Int}, t::Vector{Int})
    out = Int[]
    i, j = 1, 1
    while i ≤ length(s) || j ≤ length(t)
        if j > length(t) || (i ≤ length(s) && s[i] < t[j])
            push!(out, s[i])
            i += 1
        elseif i > length(s) || t[j] < s[i]
            push!(out, t[j])
            j += 1
        else
            i += 1
            j += 1
        end
    end
    return out
end

"""
    prefix(i)

the string of `F` on the sites before site `i`
"""
prefix(i::Int) = i > 1 ? [1, i] : Int[]

"""
    a * b

the product of two terms, site by site, `nothing` when it vanishes on a site. The strings of
`a` cross the blocks of `b`, each odd one under them bringing a sign, `F A = -A F`.
"""
function (a::Term * b::Term)
    da, db = Dict(a.blocks), Dict(b.blocks)
    pa, pb = Dict(zip(first.(a.blocks), a.parities)), Dict(zip(first.(b.blocks), b.parities))
    blocks = Pair{Int, Matrix}[]
    parities = Int[]
    for i in sort!(collect(union(keys(da), keys(db))))
        x, y = get(da, i, nothing), get(db, i, nothing)
        # a site without a block is even there
        push!(parities, get(pa, i, 0) ⊻ get(pb, i, 0))
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
    sign = crossing_sign(a.kets, a.bras, b.blocks, b.parities)
    return Term(sign * a.coef * b.coef, blocks, parities, xor_strings(a.kets, b.kets),
                xor_strings(a.bras, b.bras))
end

"""
    snapped(x, q)

the number `x` rounded to a multiple of `q`, a power of two, so that numbers equal up to rounding
become equal while exact ones, as 1 or 0.5, stay as they are
"""
snapped(x::Real, q::Real) = round(x / q) * q
snapped(x::Complex, q::Real) = complex(snapped(real(x), q), snapped(imag(x), q))

"""
    reduced(terms)

the terms with each block equal to a multiple of the identity taken out, every other block
divided by its first element above rounding, that element going into the coefficient, and its
elements rounded to a grid of rounding, see `snapped`, a term that vanishes on a site dropped,
and the terms of the same blocks gathered, a real coefficient kept real
"""
function reduced(terms::Vector{Term})
    out = Term[]
    seen = Dict{Tuple{Vector{Pair{Int, Matrix}}, Vector{Int}, Vector{Int}}, Int}()
    for t in terms
        coef = t.coef
        blocks = Pair{Int, Matrix}[]
        parities = Int[]
        for ((i, m), p) in zip(t.blocks, t.parities)
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
            b = snapped.(b, exp2(floor(log2(rounding_tol * norm(b)))))
            if eltype(b) <: Complex && all(x -> iszero(imag(x)), b)
                b = real(b)
            end
            push!(blocks, i => no_signed_zero(b))
            push!(parities, p)
        end
        if coef == 0
            continue
        end
        key = (blocks, t.kets, t.bras)
        j = get(seen, key, 0)
        if j == 0
            push!(out, Term(coef, blocks, parities, t.kets, t.bras))
            seen[key] = length(out)
        else
            out[j] = Term(out[j].coef + coef, blocks, parities, t.kets, t.bras)
        end
    end
    return [ Term(isreal(t.coef) ? real(t.coef) : t.coef, t.blocks, t.parities, t.kets, t.bras)
             for t in out if t.coef ≠ 0 ]
end

"""
    products(as, bs)

the terms of the product of two sums of terms, those that vanish left out
"""
products(as::Vector{Term}, bs::Vector{Term}) =
    Term[ t for t in (a * b for a in as for b in bs) if !isnothing(t) ]

"""
    tree_products(lists)

the terms of the product of the sums of terms `lists`, in their order, multiplied two by two as
a balanced tree and gathered at each level, see `reduced`: the terms of a long product do not
multiply with its factors, and a product of `n` factors of one site copies each block about
`log2(n)` times, where multiplying them in a row copies it up to `n` times
"""
function tree_products(lists::Vector)
    while length(lists) > 1
        lists = map(1:2:length(lists)) do k
            if k == length(lists)
                return lists[k]
            end
            ts = products(lists[k], lists[k + 1])
            return length(ts) > 1 ? reduced(ts) : ts
        end
    end
    return only(lists)
end

"""
    scaled(c, terms)

the terms with their coefficients multiplied by `c`
"""
scaled(c::Number, ts::Vector{Term}) =
    [ Term(c * t.coef, t.blocks, t.parities, t.kets, t.bras) for t in ts ]

"""
    dag(::Term)

the adjoint of a term of an operator, site by site: ``(AF)^\\dagger = FA^\\dagger``, a sign for
each odd block under its string
"""
dag(t::Term) =
    Term(crossing_sign(t.kets, Int[], t.blocks, t.parities) * conj(t.coef),
         [ i => collect(m') for (i, m) in t.blocks ], t.parities, t.kets, t.bras)

"""
    identity_term()

the term of the identity, with no block
"""
identity_term() = Term(1., Pair{Int, Matrix}[], Int[], Int[], Int[])

"""
    sided_term(t, u)

the term of the superoperator ``\\rho \\mapsto T\\rho U^\\dagger`` of the terms `t` and `u` of
operators, site by site, the strings of `t` on the kets and those of `u` on the bras: Left and
Right commute, so no sign comes
"""
function sided_term(t::Term, u::Term)
    dt, du = Dict(t.blocks), Dict(u.blocks)
    pt, pu = Dict(zip(first.(t.blocks), t.parities)), Dict(zip(first.(u.blocks), u.parities))
    sites = sort!(collect(union(keys(dt), keys(du))))
    blocks = map(sites) do i
        a, b = get(dt, i, nothing), get(du, i, nothing)
        if isnothing(a)
            a = one(b)
        elseif isnothing(b)
            b = one(a)
        end
        return i => kron(conj(b), a)
    end
    parities = [ get(pt, i, 0) + 2get(pu, i, 0) for i in sites ]
    return Term(t.coef * conj(u.coef), blocks, parities, t.kets, u.kets)
end

"""
    strings_on(t, k)

the strings of the term `t` on site `k`, that on the kets plus twice that on the bras
"""
strings_on(t::Term, k::Int) = on_string(t.kets, k) + 2on_string(t.bras, k)

"""
    string_matrix(R, f, g)

the matrix of the strings `g`, see `strings_on`, on a site whose `F` has the matrix `f`: `F`
for an operator, ``\\rho \\mapsto F^a \\rho F^b`` for a superoperator, `R = Mixed`
"""
string_matrix(::Type{Pure}, f::AbstractMatrix, g::Int) = isodd(g) ? f : one(f)
string_matrix(::Type{Mixed}, f::AbstractMatrix, g::Int) =
    kron(isodd(g >> 1) ? f : one(f), isodd(g) ? f : one(f))

"""
    gate_terms(terms)

the terms of the superoperator ``\\rho \\mapsto K\\rho K^\\dagger`` of the operator `K` of
terms `terms`, see `sided_term`
"""
gate_terms(ts::Vector{Term}) = Term[ sided_term(t, u) for t in ts for u in ts ]

"""
    laid_strings(R, system, t)

the term `t` with its strings laid as blocks, `F` on the fermionic sites of `system` they lie
on, after the block there, ``\\rho \\mapsto F^a \\rho F^b`` for a superoperator, `R = Mixed`
"""
function laid_strings(::Type{R}, system::System, t::Term) where R
    if isempty(t.kets) && isempty(t.bras)
        return t
    end
    bs = Dict(t.blocks)
    ps = Dict(zip(first.(t.blocks), t.parities))
    ends = [ t.kets ; t.bras ]
    for k in minimum(ends):maximum(ends) - 1
        f = matrix(F, system[k])
        g = strings_on(t, k)
        if f == I || g == 0
            continue
        end
        s = string_matrix(R, f, g)
        bs[k] = haskey(bs, k) ? bs[k] * s : s
        ps[k] = get(ps, k, 0)
    end
    sites = sort!(collect(keys(bs)))
    return Term(t.coef, [ k => bs[k] for k in sites ], [ ps[k] for k in sites ], Int[], Int[])
end

"""
    one_site_terms(R, system, i, m)

the terms of the operator of matrix `m` on site `i` of `system`: on a fermionic site, its part
of each parity, the odd one with `F` on the fermionic sites before `i`. For `Mixed`, `m` is the
matrix of a superoperator, whose part of parity `a` on the kets and `b` on the bras takes
``\\rho \\mapsto F^a \\rho F^b`` on the fermionic sites before `i`.
"""
function one_site_terms(::Type{Pure}, system::System, i::Int, m::AbstractMatrix)
    f = matrix(F, system[i])
    p = parity_of(m, f)
    if p == 0
        return [ Term(1., [ i => m ], [0], Int[], Int[]) ]
    elseif p == 1
        return [ Term(1., [ i => m ], [1], prefix(i), Int[]) ]
    end
    return [ one_site_terms(Pure, system, i, (m + f * m * f) / 2) ;
             one_site_terms(Pure, system, i, (m - f * m * f) / 2) ]
end

function one_site_terms(::Type{Mixed}, system::System, i::Int, m::AbstractMatrix)
    f = matrix(F, system[i])
    if f == I
        return [ Term(1., [ i => m ], [0], Int[], Int[]) ]
    end
    l, r = kron(one(f), f), kron(f, one(f))
    terms = Term[]
    for a in 0:1, b in 0:1
        x = (m + (-1)^a * l * m * l) / 2
        x = (x + (-1)^b * r * x * r) / 2
        if norm(x) ≤ rounding_tol * norm(m)
            continue
        end
        push!(terms, Term(1., [ i => x ], [a + 2b], a == 1 ? prefix(i) : Int[],
                          b == 1 ? prefix(i) : Int[]))
    end
    return terms
end

"""
    block_rep(op)

the representation whose matrices the terms of `op` have: `Pure` for an operator, `Mixed` for a
channel or a generator
"""
block_rep(::GenericOp{IsOperator}) = Pure
block_rep(::GenericOp) = Mixed

"""
    matrix_terms(system, op, index)

the terms of the operator `op` taken by its matrix on the sites `index` of `system`: on one
site, see `one_site_terms`; on several, the terms of `split_matrix`, each a product of
operators of one site in the order of `index`. A channel or a generator gives those of its
superoperator.
"""
function matrix_terms(system::System, a::GenericOp, index)
    R = block_rep(a)
    sites = [ system[i] for i in index ]
    m = matrix(a, sites...)
    if length(index) == 1
        return one_site_terms(R, system, only(index), m)
    end
    terms = Term[]
    for (c, fs) in split_matrix(R, string(a(index...)), m, sites)
        ts = scaled(c, [ identity_term() ])
        for (i, f) in zip(index, fs)
            if !isnothing(f)
                ts = products(ts, one_site_terms(R, system, i, f))
            end
        end
        append!(terms, ts)
    end
    return terms
end

"""
    pure_terms(system, k, index)

the terms of the operator `k`, generic and placed on the sites `index`, or placed, its products
gathered first, see `gathered`, their strings not laid
"""
pure_terms(system::System, k::GenericOp{IsOperator}, index) = reduced(op_terms(system, k, index))
pure_terms(system::System, k::IndexedOp{IsOperator}, _) =
    reduced(op_terms(system, gathered(system, k), ()))

"""
    trace_dual(system, t)

the term of the operator `M` the term `t` of a superoperator takes the trace through,
``\\mathrm{tr}\\,t(\\rho) = \\mathrm{tr}(M\\rho)``, block by block: a block on the strings
``\\rho \\mapsto F^a \\rho F^b`` gives ``F^b M F^a``, its string that of the kets and the bras
together, and a sign for an odd block under the string of the bras
"""
function trace_dual(system::System, t::Term)
    blocks = [ i => kraus_sum(b, [system[i]]) for (i, b) in t.blocks ]
    parities = [ (p & 1) ⊻ (p >> 1) for p in t.parities ]
    return Term(crossing_sign(t.bras, Int[], blocks, parities) * t.coef, blocks, parities,
                xor_strings(t.kets, t.bras), Int[])
end

"""
    op_terms(system, op, index)

the terms of the operator `op` on `system`, as it is written: sums and products distributed, a
generic operator placed on the sites `index`. A tensor product is the product of its factors
placed, `(A ⊗ B)(i, j) = A(i) * B(j)`, a named expression is its expression, and any other
operator is taken by its matrix, see `matrix_terms`. A channel or a generator gives the terms
of its superoperator: those of a gate, an evolver or a dissipator are made of the terms of
their operators on each side of the density matrix, see `sided_term`.
"""
op_terms(system::System, a::SumOp, index) =
    reduce(vcat, [ op_terms(system, x, index) for x in a.subs ])
op_terms(system::System, a::ProdOp, index) =
    tree_products([ op_terms(system, x, index) for x in a.subs ])
op_terms(system::System, a::ScalarOp, index) = scaled(a.coef, op_terms(system, a.arg, index))
op_terms(::System, ::IdentityOp, _) = [ identity_term() ]
op_terms(system::System, a::AtIndex, _) = op_terms(system, a.op, a.index)
op_terms(system::System, a::TensorOp, index) =
    tree_products([ op_terms(system, o, index[p]) for (o, p) in zip(a.subs, factor_sites(a)) ])
op_terms(system::System, a::DagOp, index) = map(dag, op_terms(system, a.arg, index))
op_terms(system::System, a::Operator, index) =
    a.expr isa GenericOp ? op_terms(system, a.expr, index) : matrix_terms(system, a, index)
op_terms(system::System, a::Union{Proj, ExpOp, IntPowOp, GenPowOp{IsOperator}, ModOp, Reset,
                                  Dephase}, index) =
    matrix_terms(system, a, index)

op_terms(system::System, a::Gate, index) = gate_terms(pure_terms(system, a.arg, index))

function op_terms(system::System, a::Evolver, index)
    ts = pure_terms(system, a.arg, index)
    one = identity_term()
    return Term[ [ sided_term(t, one) for t in ts ] ; [ sided_term(one, t) for t in ts ] ]
end

# ρ ↦ Σ K ρ K† - {M, ρ}/2, M = Σ K†K, from the terms of L, or those of a channel and the M they
# give block by block, a tensor product of superoperators giving the tensor product of their M
function op_terms(system::System, a::Dissipator, index)
    if a.arg isa Op{IsOperator}
        ts = pure_terms(system, a.arg, index)
        jumps = gate_terms(ts)
        m = reduced(products(map(dag, ts), ts))
    else
        jumps = op_terms(system, a.arg, index)
        m = reduced(map(t -> trace_dual(system, t), jumps))
    end
    one = identity_term()
    return Term[ jumps ; scaled(-0.5, [ [ sided_term(x, one) for x in m ] ;
                                        [ sided_term(one, x) for x in m ] ]) ]
end

"""
    strong_refusal(what, names)

refuse `what`, which changes the strongly conserved `names` by different amounts from one part
to another, mixing the sectors a strong symmetry keeps apart
"""
strong_refusal(what, names) =
    error("$what changes $(join(names, ", ")) by different amounts from one part to another, " *
          "which conserving it strongly forbids: drop `strong` or weaken the state")

"""
    mixed_flux(m, what, site)

the flux on the mixed index of `site` of the matrix `m` of the superoperator `what`, refused
when it has none: as changing a strongly conserved charge by different amounts, which would mix
the sectors a strong symmetry keeps apart, or as having no definite charge at all
"""
function mixed_flux(m::AbstractMatrix, what, site::AbstractSite)
    qs = conserved(site).quantities
    if isempty(qs)
        return QN()
    end
    small = rounding_tol * norm(m)
    found = unique(g.q for (g, x) in zip(element_grades(Mixed, site), vec(m)) if abs(x) > small)
    if length(found) ≤ 1
        return isempty(found) ? QN() : only(found)
    end
    st = strong_names(site)
    if !isempty(st) && allequal(weak_qn(q, st, String[]) for q in found)
        strong_refusal(what, st)
    end
    error("$what carries no definite charge of $(join(unique(q.name for q in qs), ", ")), so it " *
          "cannot act on $site")
end

"""
    block_flux(R, m, what, site)

the flux of the matrix `m` of a term of `what` on `site`, that of an operator for `Pure`, of a
superoperator for `Mixed`
"""
block_flux(::Type{Pure}, m::AbstractMatrix, what, site::AbstractSite) = charge_flux(m, what, site)
block_flux(::Type{Mixed}, m::AbstractMatrix, what, site::AbstractSite) = mixed_flux(m, what, site)

"""
    check_flux(R, system, a, terms)

refuse the operator `a` on a charged `system` when a block of its terms has no definite flux,
or when its terms do not all change the charges by the same amount, the matrices being those
of superoperators for `R = Mixed`
"""
function check_flux(::Type{R}, system::System, a, terms::Vector{Term}) where R
    if !is_charged(system)
        return nothing
    end
    st = strong_names(system)
    total = nothing
    for t in terms
        q = sum((block_flux(R, m, a, system[i]) for (i, m) in t.blocks); init = QN())
        if isnothing(total)
            total = q
        elseif q ≠ total
            if !isempty(st) && weak_qn(q, st, String[]) == weak_qn(total, st, String[])
                strong_refusal(a, st)
            end
            error("the terms of $a change the charges by $total and by $q, so it has no " *
                  "definite flux and cannot act on a system that conserves them")
        end
    end
    return nothing
end

"""
    system_matrix(R, system, terms)

the matrix of a sum of terms on the whole of `system`, the last site varying fastest, those of
a superoperator for `Mixed`, see `matrix`, their strings laid
"""
function system_matrix(::Type{R}, system::System, ts::Vector{Term}) where R
    ids = [ Matrix{Float64}(I, block_dim(R, system[k]), block_dim(R, system[k]))
            for k in 1:length(system) ]
    d = prod(size(m, 1) for m in ids)
    total = zeros(d, d)
    for x in ts
        t = laid_strings(R, system, x)
        bs = Dict(t.blocks)
        total += t.coef * foldl(kron, [ get(bs, k, ids[k]) for k in eachindex(ids) ])
    end
    return total
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
    return system_matrix(Pure, sub, op_terms(sub, map_sites(i -> position[i], a), ()))
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
    term_count(system, op)

the number of terms of the placed operator `op` on `system`
"""
term_count(system::System, a::IndexedOp{IsOperator}) = length(pure_terms(system, a, ()))

"""
    count_add!(tree, k, x)

add `x` to the count of position `k` of `tree`, a Fenwick tree of counts
"""
function count_add!(tree::Vector{Int}, k::Int, x::Int)
    while k ≤ length(tree)
        tree[k] += x
        k += k & -k
    end
    return tree
end

"""
    count_to(tree, k)

the sum of the counts of the positions 1 to `k` of the Fenwick tree `tree`
"""
function count_to(tree::Vector{Int}, k::Int)
    c = 0
    while k > 0
        c += tree[k]
        k -= k & -k
    end
    return c
end

"""
    last_counted(tree)

the last position of the Fenwick tree `tree` of counts 0 or 1 whose count is 1, 0 when there is
none
"""
function last_counted(tree::Vector{Int})
    c = count_to(tree, length(tree))
    if c == 0
        return 0
    end
    # the last position whose sum up to it is under c
    k, step = 0, prevpow(2, length(tree))
    while step > 0
        if k + step ≤ length(tree) && tree[k + step] < c
            k += step
            c -= tree[k]
        end
        step >>= 1
    end
    return k + 1
end

"""
    heap_push!(h, x, less)
    heap_pop!(h, less)

add `x` to the binary heap `h` ordered by `less`, or take its least element out
"""
function heap_push!(h::Vector, x, less)
    push!(h, x)
    i = length(h)
    while i > 1 && less(h[i], h[i >> 1])
        h[i], h[i >> 1] = h[i >> 1], h[i]
        i >>= 1
    end
    return h
end

function heap_pop!(h::Vector, less)
    top, x = h[1], pop!(h)
    if !isempty(h)
        h[1] = x
        i = 1
        while true
            m = i
            for c in (2i, 2i + 1)
                if c ≤ length(h) && less(h[c], h[m])
                    m = c
                end
            end
            if m == i
                break
            end
            h[i], h[m] = h[m], h[i]
            i = m
        end
    end
    return top
end

"""
    sort_factors!(fs, info)

sort the factors `fs` of a product, and their `(sites, parity)` in `info`, by their sites as
far as they (anti)commute, and return the sign that brings, -1 for each pair of odd factors that
change places: of the factors whose dependencies before them are placed, the one of smallest
sites goes next, so that products differing only by the order of factors that (anti)commute
take one form. Two factors depend on each other when they share a site, or when one has no
definite parity and the other is odd or has none; on other sites, they commute, up to a sign
when both are odd.
"""
function sort_factors!(fs, info)
    n = length(fs)
    order_keys = [ sort(sites) for (sites, _) in info ]
    # the factors each one waits for, enough for all it depends on: the last one on each of its
    # sites, the last of no definite parity for an odd one or one of none, and the odd ones
    # since for one of none
    waits = zeros(Int, n)
    next = [ Int[] for _ in 1:n ]
    last_on = Dict{Int, Int}()
    none, odds = 0, Int[]
    for (j, (sites, p)) in enumerate(info)
        before = Set{Int}(last_on[k] for k in sites if haskey(last_on, k))
        if p ≠ 0 && none > 0
            push!(before, none)
        end
        if isnothing(p)
            union!(before, odds)
            none, odds = j, Int[]
        elseif p == 1
            push!(odds, j)
        end
        foreach(k -> push!(next[k], j), before)
        waits[j] = length(before)
        foreach(k -> last_on[k] = j, sites)
    end
    less(a, b) = (order_keys[a], a) < (order_keys[b], b)
    ready = Int[]
    foreach(j -> heap_push!(ready, j, less), findall(iszero, waits))
    order = Int[]
    while !isempty(ready)
        j = heap_pop!(ready, less)
        push!(order, j)
        for k in next[j]
            waits[k] -= 1
            if waits[k] == 0
                heap_push!(ready, k, less)
            end
        end
    end
    # -1 for each pair of odd factors in the other order, counted in a Fenwick tree
    placed, swaps, tree = 0, 0, zeros(Int, n)
    for j in order
        if info[j][2] == 1
            swaps += placed - count_to(tree, j)
            count_add!(tree, j, 1)
            placed += 1
        end
    end
    copyto!(fs, fs[order])
    copyto!(info, info[order])
    return isodd(swaps) ? -1 : 1
end

"""
    merged_pass(system, fs, info, refused)

the factors of the product `fs`, of `(sites, parity)` `info`, each merged, from left to right,
into the nearest earlier one it does not (anti)commute with, see `sort_factors!`, if it shares a
site with it, the sites of one among those of the other, see `merged`, and the merge gives no
more terms than the two factors, see `term_count`, a pair giving more going into `refused`.
Return the factors left with their info, the number the product gets multiplied by, -1 for
each odd factor an odd one goes back past, and whether anything merged. That nearest factor is
read off the last factor on each site and the last odd one and the last one of no definite
parity, counted in Fenwick trees.
"""
function merged_pass(system::System, fs, info, refused)
    n = length(fs)
    out = Union{Nothing, IndexedOp{IsOperator}}[]
    oinfo = Tuple{Vector{Int}, Union{Nothing, Int}}[]
    on_site = Dict{Int, Vector{Int}}()
    odd, none = zeros(Int, n), zeros(Int, n)
    coef, changed = 1, false
    # the last factor still there on site s
    function last_on(s)
        ks = get!(Vector{Int}, on_site, s)
        while !isempty(ks) && isnothing(out[last(ks)])
            pop!(ks)
        end
        return isempty(ks) ? 0 : last(ks)
    end
    function tally!(k, p, x)
        if isnothing(p)
            count_add!(none, k, x)
        elseif p == 1
            count_add!(odd, k, x)
        end
    end
    for (x, (sx, px)) in zip(fs, info)
        i = maximum(last_on, sx; init = 0)
        blocking = px == 0 ? 0 : max(last_counted(none), px == 1 ? 0 : last_counted(odd))
        if i > 0 && i ≥ blocking && out[i] isa AtIndex && x isa AtIndex &&
           (issubset(oinfo[i][1], sx) || issubset(sx, oinfo[i][1])) && (out[i], x) ∉ refused
            a = out[i]
            m = merged(system, a, x)
            if m isa AtIndex && term_count(system, m) > term_count(system, a) * term_count(system, x)
                push!(refused, (a, x))
            else
                changed = true
                if px == 1 && isodd(count_to(odd, n) - count_to(odd, i))
                    coef = -coef
                end
                tally!(i, oinfo[i][2], -1)
                if m isa AtIndex
                    sm, pm = sites_of(m), parity_on(system, m)
                    for s in setdiff(sm, oinfo[i][1])
                        push!(on_site[s], i)
                    end
                    out[i], oinfo[i] = m, (sm, pm)
                    tally!(i, pm, 1)
                else
                    out[i] = nothing
                    coef *= scalarcoef(m)
                end
                continue
            end
        end
        push!(out, x)
        push!(oinfo, (sx, px))
        foreach(s -> push!(get!(Vector{Int}, on_site, s), length(out)), sx)
        tally!(length(out), px, 1)
    end
    kept = [ k for k in eachindex(out) if !isnothing(out[k]) ]
    return IndexedOp{IsOperator}[ out[k] for k in kept ], oinfo[kept], coef, changed
end

"""
    gathered(system, op)

the placed operator `op` rewritten without developing anything, into a normal form: in each
product, two factors sharing a site are merged where they can be, see `merged_pass`, until
nothing merges, then the factors are sorted by their sites as far as they (anti)commute, see
`sort_factors!`. Sums are gathered term by term, and an operator in normal form is left as it
is.
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
    info = [ (sites_of(x), parity_on(system, x)) for x in fs ]
    refused = Set{Tuple{IndexedOp{IsOperator}, IndexedOp{IsOperator}}}()
    while coef ≠ 0
        fs, info, c, changed = merged_pass(system, fs, info, refused)
        coef *= c
        if !changed
            break
        end
    end
    coef *= sort_factors!(fs, info)
    return scale(coef, ProdOp(fs))
end

"""
    lower(system, op)

the terms of the operator `op`, placed on sites, on `system`, see `Term`, once its products are
gathered, see `gathered`, their Jordan-Wigner strings kept apart, see `laid_strings`, the
operator refused on a charged system when it has no definite flux. A channel or a generator
gives the terms of its superoperator on the density matrices of the sites, refused likewise.
"""
function lower(system::System, a::IndexedOp{IsOperator})
    # before any rewriting, so that a message names the operator as it was written
    check_indices(system, a)
    terms = pure_terms(system, a, ())
    check_flux(Pure, system, a, terms)
    return terms
end

function lower(system::System, a::IndexedOp{<:Union{IsChannel, IsGenerator}})
    check_indices(system, a)
    terms = reduced(op_terms(system, a, ()))
    check_flux(Mixed, system, a, terms)
    return terms
end
