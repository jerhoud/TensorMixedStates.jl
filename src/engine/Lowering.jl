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
        error("$(obs_name(a)) acts on site $i, which the system does not have: it has $n sites")
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
    Block

the matrix `m` of a term on its site `site`, see `Term`, and its fermionic parity, that of the
kets plus twice that of the bras for a superoperator, which gives the sign of a string crossing
it
"""
struct Block
    site::Int
    m::Matrix
    parity::Int
end

(a::Block == b::Block) = a.site == b.site && a.m == b.m && a.parity == b.parity
hash(b::Block, h::UInt) = hash(b.parity, hash(b.m, hash(b.site, h)))

"""
    Term

a term of an operator on a system: `coef` times the tensor product of its `blocks`, see `Block`,
one on each site it acts on, sorted by site, the identity on the others, each followed by the
Jordan-Wigner strings of the term there, which the descent keeps apart, laid as blocks by what
needs them, see `laid_strings`. The term of a superoperator has the matrices of superoperators
of one site, see `matrix`, and a string on the kets, ``\\rho \\mapsto F\\rho``, and one on the
bras, ``\\rho \\mapsto \\rho F``.

`kets` and `bras` are the strings, as the sorted sites where their parity changes: `F` on the
sites from the first to the second excluded, from the third to the fourth, and so on. Each of
these sites is site 1 or a site where the term has a block, which `string_factors` relies on.
"""
struct Term
    coef::Number
    blocks::Vector{Block}
    kets::Vector{Int}
    bras::Vector{Int}
end

"""
    on_string(s, k)

whether site `k` lies on the string `s`, see `Term`
"""
on_string(s::Vector{Int}, k::Int) = isodd(searchsortedlast(s, k))

"""
    crossing_sign(kets, bras, blocks)

the sign the strings `kets` and `bras` bring crossing `blocks`: -1 for each block odd on the
side of a string over its site, `F A = -A F`
"""
crossing_sign(kets::Vector{Int}, bras::Vector{Int}, blocks) =
    (-1)^count((isodd(b.parity) && on_string(kets, b.site)) ⊻
               (isodd(b.parity >> 1) && on_string(bras, b.site)) for b in blocks)

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
    paired_blocks(as, bs)

the blocks `as` and `bs`, sorted by site, paired on each site where one of them has a block,
in order, `nothing` standing for the one that has none
"""
function paired_blocks(as::Vector{Block}, bs::Vector{Block})
    out = Tuple{Union{Block, Nothing}, Union{Block, Nothing}}[]
    i, j, na, nb = 1, 1, length(as), length(bs)
    while i ≤ na || j ≤ nb
        k = min(i ≤ na ? as[i].site : typemax(Int), j ≤ nb ? bs[j].site : typemax(Int))
        x = i ≤ na && as[i].site == k ? as[i] : nothing
        y = j ≤ nb && bs[j].site == k ? bs[j] : nothing
        push!(out, (x, y))
        i += !isnothing(x)
        j += !isnothing(y)
    end
    return out
end

"""
    a * b

the product of two terms, site by site, `nothing` when it vanishes on a site. The strings of
`a` cross the blocks of `b`, each odd one under them bringing a sign, `F A = -A F`.
"""
function (a::Term * b::Term)
    blocks = Block[]
    for (x, y) in paired_blocks(a.blocks, b.blocks)
        if isnothing(x) || isnothing(y)
            push!(blocks, something(x, y))
            continue
        end
        m = x.m * y.m
        # a product of factors that vanishes, as Sp * Sp on a spin 1/2
        if norm(m) ≤ rounding_tol * norm(x.m) * norm(y.m)
            return nothing
        end
        push!(blocks, Block(x.site, m, x.parity ⊻ y.parity))
    end
    sign = crossing_sign(a.kets, a.bras, b.blocks)
    return Term(sign * a.coef * b.coef, blocks, xor_strings(a.kets, b.kets),
                xor_strings(a.bras, b.bras))
end

"""
    representatives(blocks)

for each of the matrices `blocks`, a factor and a representative, the block being the factor
times its representative up to rounding, the blocks proportional up to rounding sharing theirs:
the first of them in an order of their own, divided by its first element of largest modulus and
cleared of rounding. The factor of another block is its element there, the whole block then
checked against the representative, so that no two blocks have to choose the same element, as
in the deparallelisation of Hubig, McCulloch and Schollwöck, Phys. Rev. B 95, 035129 (2017).
They are sorted by the projection of their direction on fixed random weights, the same for a
block and its multiples, which proportional blocks have within a window of the width of
rounding, each compared with the representatives of those of its window before it alone
"""
function representatives(blocks::Vector{Matrix})
    factors = Vector{Number}(undef, length(blocks))
    reps = similar(blocks)
    # a representative before it is cleared of rounding, and the position of its 1
    raw = similar(blocks)
    pivots = zeros(Int, length(blocks))
    for n in unique(length.(blocks))
        group = findall(b -> length(b) == n, blocks)
        # fixed, so that an operator always gives the same MPO, and drawn from a generator of
        # their own, which leaves the draws of the user as they are
        w = randn(Xoshiro(20261015), 2n)
        v = complex.(w[1:n], w[n+1:2n])
        p = [ abs2(sum(v[j] * b[j] for j in 1:n)) / sum(abs2, b) for b in blocks[group] ]
        δ = 4 * sum(abs2, v) * rounding_tol
        order = sortperm(p)
        for (a, g) in enumerate(order)
            i = group[g]
            m = blocks[i]
            for c in a-1:-1:1
                if p[g] - p[order[c]] > δ
                    break
                end
                j = group[order[c]]
                k = pivots[j]
                if size(raw[j]) == size(m) && norm(m - m[k] * raw[j]) ≤ rounding_tol * norm(m)
                    factors[i], reps[i], raw[i], pivots[i] = m[k], reps[j], raw[j], k
                    break
                end
            end
            if pivots[i] == 0
                k = argmax(abs.(vec(m)))
                r = m / m[k]
                r[k] = 1
                b = copy(r)
                # rounding between charges, which an eigendecomposition leaves behind, is no
                # charge
                b[abs.(b) .≤ rounding_tol * norm(b)] .= 0
                # as the phases a split by singular values leaves on its factors, see
                # `split_matrix`
                if eltype(b) <: Complex && all(x -> abs(imag(x)) ≤ rounding_tol * norm(b), b)
                    b = real(b)
                end
                factors[i], reps[i], raw[i], pivots[i] = m[k], no_signed_zero(b), r, k
            end
        end
    end
    return factors, reps
end

"""
    reduced(terms)

the terms with each block equal to a multiple of the identity taken out, every other block
written as a factor, which goes into the coefficient, times the representative of the blocks
proportional to it up to rounding, see `representatives`, a term that vanishes on a site
dropped, and the terms of the same blocks gathered, a real coefficient kept real
"""
function reduced(terms::Vector{Term})
    normed = Tuple{Number, Vector{Block}, Term}[]
    for t in terms
        coef = t.coef
        blocks = Block[]
        for x in t.blocks
            m = x.m
            λ = tr(m) / size(m, 1)
            if iszero(m)
                coef = 0
                break
            elseif nearly(m, Matrix(λ * I, size(m)))
                coef *= λ
                continue
            end
            push!(blocks, x)
        end
        if coef ≠ 0
            push!(normed, (coef, blocks, t))
        end
    end
    factors, reps = representatives(Matrix[ x.m for (_, bs, _) in normed for x in bs ])
    out = Term[]
    seen = Dict{Tuple{Vector{Block}, Vector{Int}, Vector{Int}}, Int}()
    r = 0
    for (coef, bs, t) in normed
        blocks = Block[]
        for x in bs
            r += 1
            coef *= factors[r]
            push!(blocks, Block(x.site, reps[r], x.parity))
        end
        key = (blocks, t.kets, t.bras)
        j = get(seen, key, 0)
        if j == 0
            push!(out, Term(coef, blocks, t.kets, t.bras))
            seen[key] = length(out)
        else
            out[j] = Term(out[j].coef + coef, blocks, t.kets, t.bras)
        end
    end
    return [ Term(isreal(t.coef) ? real(t.coef) : t.coef, t.blocks, t.kets, t.bras)
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
scaled(c::Number, ts::Vector{Term}) = [ Term(c * t.coef, t.blocks, t.kets, t.bras) for t in ts ]

# the adjoint of a term of an operator, site by site, (AF)† = FA†, a sign for each odd block
# under its string; a comment, not a docstring, which would put a method of no user in the
# reference of the exported dag
dag(t::Term) =
    Term(crossing_sign(t.kets, Int[], t.blocks) * conj(t.coef),
         [ Block(x.site, collect(x.m'), x.parity) for x in t.blocks ], t.kets, t.bras)

"""
    identity_term()

the term of the identity, with no block
"""
identity_term() = Term(1., Block[], Int[], Int[])

"""
    sided_term(t, u)

the term of the superoperator ``\\rho \\mapsto T\\rho U^\\dagger`` of the terms `t` and `u` of
operators, site by site, the strings of `t` on the kets and those of `u` on the bras: Left and
Right commute, so no sign comes
"""
function sided_term(t::Term, u::Term)
    blocks = map(paired_blocks(t.blocks, u.blocks)) do (x, y)
        a = isnothing(x) ? one(y.m) : x.m
        b = isnothing(y) ? one(x.m) : y.m
        return Block(something(x, y).site, kron(conj(b), a),
                     (isnothing(x) ? 0 : x.parity) + 2(isnothing(y) ? 0 : y.parity))
    end
    return Term(t.coef * conj(u.coef), blocks, t.kets, u.kets)
end

"""
    side_terms(side, terms)

the terms of the superoperator ``\\rho \\mapsto K\\rho``, for `side = :left`, or
``\\rho \\mapsto \\rho K^\\dagger``, for `:right`, of the operator `K` of terms `terms`, see
`sided_term`
"""
function side_terms(side::Symbol, ts::Vector{Term})
    one = identity_term()
    return Term[ side == :left ? sided_term(t, one) : sided_term(one, t) for t in ts ]
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
    bs = Dict(x.site => x for x in t.blocks)
    ends = [ t.kets ; t.bras ]
    for k in minimum(ends):maximum(ends) - 1
        f = matrix(F, system[k])
        g = strings_on(t, k)
        if f == I || g == 0
            continue
        end
        s = string_matrix(R, f, g)
        bs[k] = haskey(bs, k) ? Block(k, bs[k].m * s, bs[k].parity) : Block(k, s, 0)
    end
    return Term(t.coef, [ bs[k] for k in sort!(collect(keys(bs))) ], Int[], Int[])
end

"""
    grade_parts(R, m, site)

the parts of the matrix `m` on `site`, of an operator for `Pure` or of a superoperator for
`Mixed`, each made of its elements of one grade, see `element_grades`, as pairs `(p, part)`, `p`
the fermionic parity of the part: `m` alone when its elements have a single grade. What
rounding leaves of another grade is dropped, as charged indices drop it
"""
function grade_parts(::Type{R}, m::AbstractMatrix, site::AbstractSite) where R
    gs = element_grades(R, site)
    small = rounding_tol * norm(m)
    found = unique(g for (g, x) in zip(gs, vec(m)) if abs(x) > small)
    if length(found) ≤ 1
        return [ (isempty(found) ? 0 : only(found).p, m) ]
    end
    return [ (g.p, reshape([ h == g ? x : zero(x) for (h, x) in zip(gs, vec(m)) ], size(m)))
             for g in found ]
end

"""
    one_site_terms(R, system, i, m)

the terms of the operator of matrix `m` on site `i` of `system`, its part of each grade, charge
and parity, see `grade_parts`, an odd part taking `F` on the fermionic sites before `i`: for
`Mixed`, `m` is the matrix of a superoperator, whose part of parity `a` on the kets and `b` on
the bras takes ``\\rho \\mapsto F^a \\rho F^b`` there. A block connecting several charges, as
that of `Sx`, gives terms of definite charge, which a sum may cancel.
"""
one_site_terms(::Type{R}, system::System, i::Int, m::AbstractMatrix) where R =
    [ Term(1., [ Block(i, x, p) ], isodd(p) ? prefix(i) : Int[], p ≥ 2 ? prefix(i) : Int[])
      for (p, x) in grade_parts(R, m, system[i]) ]

"""
    block_rep(op)

the representation whose matrices the terms of `op` have: `Pure` for an operator, `Mixed` for a
channel or a generator
"""
block_rep(::GenericOp{IsOperator}) = Pure
block_rep(::GenericOp) = Mixed

"""
    matrix_terms(system, op, index; name)

the terms of the operator `op` taken by its matrix on the sites `index` of `system`: on one
site, see `one_site_terms`; on several, the terms of `split_matrix`, each a product of
operators of one site in the order of `index`, `name` naming it in a refusal. A channel or a
generator gives those of its superoperator.
"""
function matrix_terms(system::System, a::GenericOp, index; name::String = string(a(index...)))
    R = block_rep(a)
    sites = [ system[i] for i in index ]
    m = matrix(a, sites...)
    if length(index) == 1
        return one_site_terms(R, system, only(index), m)
    end
    terms = Term[]
    for (c, fs) in split_matrix(R, name, m, sites)
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
    blocks = [ Block(x.site, kraus_sum(x.m, [system[x.site]]), (x.parity & 1) ⊻ (x.parity >> 1))
               for x in t.blocks ]
    return Term(crossing_sign(t.bras, Int[], blocks) * t.coef, blocks, xor_strings(t.kets, t.bras),
                Int[])
end

"""
    op_terms(system, op, index)

the terms of the operator `op` on `system`, as it is written: sums and products distributed, a
generic operator placed on the sites `index`. A tensor product is the product of its factors
placed, `(A ⊗ B)(i, j) = A(i) * B(j)`, and a named expression is its expression. A generic sum
or a named expression of one site is taken by its matrix instead, a term for each parity rather
than one for each operand, and so is any other operator, see `matrix_terms`. A channel or a
generator gives the terms of its superoperator: those of a gate, an evolver or a dissipator are
made of the terms of their operators on each side of the density matrix, see `sided_term`, and
those of a dephasing are the gates of its projectors when they are few, see `projector_gates`.
"""
op_terms(system::System, a::SumOp, index) =
    reduce(vcat, [ op_terms(system, x, index) for x in a.subs ])
op_terms(system::System, a::SumOp{IsOperator, Generic, 1}, index) = matrix_terms(system, a, index)
op_terms(system::System, a::ProdOp, index) =
    tree_products([ op_terms(system, x, index) for x in a.subs ])
op_terms(system::System, a::ScalarOp, index) = scaled(a.coef, op_terms(system, a.arg, index))
op_terms(::System, ::IdentityOp, _) = [ identity_term() ]
op_terms(system::System, a::AtIndex, _) = op_terms(system, a.op, a.index)
op_terms(system::System, a::TensorOp, index) =
    tree_products([ op_terms(system, o, index[p]) for (o, p) in zip(a.subs, factor_sites(a)) ])
op_terms(system::System, a::DagOp, index) = map(dag, op_terms(system, a.arg, index))
op_terms(system::System, a::Operator{N}, index) where N =
    N > 1 && a.expr isa GenericOp ? op_terms(system, a.expr, index) : matrix_terms(system, a, index)
op_terms(system::System, a::Union{Proj, ExpOp, IntPowOp, GenPowOp{IsOperator}, ModOp, Reset},
         index) =
    matrix_terms(system, a, index)

op_terms(system::System, a::Gate, index) = gate_terms(pure_terms(system, a.arg, index))

# the projectors are split as operators, of dimension d, where the superoperator, of dimension
# d², took a minute to split on six qubits
"""
    projector_gates(system, a, index)

the terms of the dephasing `a` on the sites `index` of `system` as the gates of its projectors,
`nothing` on a single site or when they would outnumber the dimension of its superoperator, as
for an operator given by a generic matrix
"""
function projector_gates(system::System, a::Dephase, index)
    if length(index) == 1
        return nothing
    end
    sites = [ system[i] for i in index ]
    values, _ = eigenspaces(matrix(a.arg, sites...), "Dephase", a.arg,
                            " on " * join(sites, " ⊗ "))
    name = string(a(index...))
    tss = [ reduced(matrix_terms(system, Proj(a.arg => λ), index; name)) for λ in values ]
    if sum(ts -> length(ts)^2, tss) > prod(s -> dim(s)^2, sites)
        return nothing
    end
    return reduce(vcat, map(gate_terms, tss))
end

function op_terms(system::System, a::Dephase, index)
    ts = projector_gates(system, a, index)
    return isnothing(ts) ? matrix_terms(system, a, index) : ts
end

function op_terms(system::System, a::Evolver, index)
    ts = pure_terms(system, a.arg, index)
    return Term[ side_terms(:left, ts) ; side_terms(:right, ts) ]
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
    return Term[ jumps ; scaled(-0.5, [ side_terms(:left, m) ; side_terms(:right, m) ]) ]
end

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
    names = join(unique(q.name for q in qs), ", ")
    charge_refusal(what, strong_names(site), found,
                   "$(obs_name(what)) carries no definite charge of $names, so it cannot act " *
                   "on $site")
end

"""
    block_flux(R, m, what, site)

the flux of the matrix `m` of a term of `what` on `site`, that of an operator for `Pure`, of a
superoperator for `Mixed`
"""
block_flux(::Type{Pure}, m::AbstractMatrix, what, site::AbstractSite) = charge_flux(m, what, site)
block_flux(::Type{Mixed}, m::AbstractMatrix, what, site::AbstractSite) = mixed_flux(m, what, site)

"""
    normalized_norm(m)

the Hilbert-Schmidt norm of the matrix `m` over that of the identity, 1 for any unitary, so
that the norm of a tensor product, the product of those of its factors, does not grow with
their number
"""
normalized_norm(m::AbstractMatrix) = norm(m) / sqrt(size(m, 1))

"""
    normalized_size(t)

the coefficient of the term `t` times the normalized norms of its blocks, see `normalized_norm`,
next to which what rounding leaves of a sum of terms is measured
"""
normalized_size(t::Term) = abs(t.coef) * prod(normalized_norm(x.m) for x in t.blocks; init = 1.0)

"""
    term_flux(R, system, a, t)

the amount the term `t` of the operator `a` changes the charges of `system` by, the sum of the
fluxes of its blocks, refused when one has none, see `block_flux`
"""
term_flux(::Type{R}, system::System, a, t::Term) where R =
    sum((block_flux(R, x.m, a, system[x.site]) for x in t.blocks); init = QN())

"""
    flux_refusal(system, a, fluxes)

refuse the operator `a`, whose parts change the charges of `system` by the different amounts
`fluxes`, see `FluxError`, or by amounts differing only in what it conserves strongly, see
`strong_refusal`
"""
flux_refusal(system::System, a, fluxes) =
    charge_refusal(a, strong_names(system), fluxes,
                   "the terms of $(obs_name(a)) change the charges by $(fluxes[1]) and by " *
                   "$(fluxes[2]), so it has no definite flux and cannot act on a system that " *
                   "conserves them")

"""
    check_flux(R, system, a, terms)

the terms of the operator `a` on a charged `system`, the matrices being those of superoperators
for `R = Mixed`, refused when a block has no definite flux or when the terms do not all change
the charges by the same amount, see `FluxError`. A term of rounding size next to the largest
one, what is left of a sum cancelling up to rounding, is dropped rather than counted, see
`normalized_size`.
"""
function check_flux(::Type{R}, system::System, a, terms::Vector{Term}) where R
    if !is_charged(system)
        return terms
    end
    fluxes = [ term_flux(R, system, a, t) for t in terms ]
    sizes = normalized_size.(terms)
    small = rounding_tol * maximum(sizes; init = 0.)
    main = unique(q for (q, s) in zip(fluxes, sizes) if s > small)
    if length(main) > 1
        flux_refusal(system, a, main)
    end
    return [ t for (t, q) in zip(terms, fluxes) if isempty(main) || q == only(main) ]
end

"""
    system_matrix(R, system, terms, sites = 1:length(system), orders = nothing)

the matrix of a sum of terms on `sites` of `system`, which hold all their blocks, the last site
varying fastest, their strings laid, that of a superoperator for `Mixed`, see `matrix`, the
elements of each site in the order `orders` gives for it, as on its mixed index, see
`mixed_elements`
"""
function system_matrix(::Type{R}, system::System, ts::Vector{Term}, sites = 1:length(system),
                       orders = nothing) where R
    ids = [ Matrix{Float64}(I, block_dim(R, system[k]), block_dim(R, system[k])) for k in sites ]
    d = prod(size(m, 1) for m in ids)
    total = zeros(d, d)
    for x in ts
        t = laid_strings(R, system, x)
        bs = Dict(b.site => b.m for b in t.blocks)
        ms = [ get(bs, k, ids[j]) for (j, k) in enumerate(sites) ]
        if !isnothing(orders)
            ms = [ m[o, o] for (m, o) in zip(ms, orders) ]
        end
        total += t.coef * foldl(kron, ms)
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
Jordan-Wigner strings between them included: its products and sums taken as matrices, a factor
on exactly `sites` by its matrix, any other by its terms on these sites alone, see
`system_matrix`
"""
function matrix_on(system::System, a::AtIndex, sites)
    if collect(a.index) == collect(sites)
        return matrix(a.op, [ system[i] for i in sites ]...)
    end
    sub = System([ system[i] for i in sites ])
    position = Dict(i => k for (k, i) in enumerate(sites))
    return system_matrix(Pure, sub, op_terms(sub, map_sites(i -> position[i], a), ()))
end

matrix_on(system::System, a::ProdOp, sites) = prod(matrix_on(system, x, sites) for x in a.subs)
matrix_on(system::System, a::SumOp, sites) = sum(matrix_on(system, x, sites) for x in a.subs)
matrix_on(system::System, a::ScalarOp, sites) = a.coef * matrix_on(system, a.arg, sites)

function matrix_on(system::System, ::IdentityOp, sites)
    d = prod(dim(system[k]) for k in sites)
    return Matrix{Float64}(I, d, d)
end

"""
    placed_matrix(m, magnitude, name, sites)

the operator `name` of matrix `m` placed on `sites`, or a number when `m` vanishes next to
`magnitude`, the size of what it was computed from, or is a multiple of the identity
"""
function placed_matrix(m::AbstractMatrix, magnitude::Real, name::String, sites)
    λ = tr(m) / size(m, 1)
    if norm(m) ≤ rounding_tol * magnitude
        return scale(0, IdentityOp{IsOperator, Indexed, 1}())
    elseif nearly(m, Matrix(λ * I, size(m)))
        return scale(λ, IdentityOp{IsOperator, Indexed, 1}())
    end
    return Operator{length(sites)}(name, m)(sites...)
end

"""
    merged(system, a, b)

the product `a * b` of two placed operators, the sites of one among those of the other, as one
operator given by its matrix on the sites of the larger, see `placed_matrix`
"""
function merged(system::System, a::AtIndex, b::AtIndex)
    sites = length(a.index) ≥ length(b.index) ? a.index : b.index
    x, y = matrix_on(system, a, sites), matrix_on(system, b, sites)
    return placed_matrix(x * y, norm(x) * norm(y), "($(a * b))", sites)
end

"""
    term_count(system, op)

the number of terms of the placed operator `op` on `system`
"""
term_count(system::System, a::IndexedOp{IsOperator}) = length(pure_terms(system, a, ()))

"""
    or_on_plain(f, xs...)

`f(xs...)`, or, when that refuses an operator changing the charges by no definite amount, see
`FluxError`, `f` of the states or the systems `xs` conserving nothing, see `conserving_nothing`
"""
function or_on_plain(f, xs...)
    try
        return f(xs...)
    catch e
        if !(e isa FluxError)
            rethrow()
        end
    end
    return f(map(conserving_nothing, xs)...)
end

"""
    conserving_nothing(x)

the state `x` conserving nothing, see `plain_form`, or the system `x` with the same sites
conserving nothing
"""
conserving_nothing(x::State) = plain_form(x)
conserving_nothing(x::System) = weaken(x, ())

"""
    more_terms(system, m, a, x)

whether the factor `m`, the merge of the factors `a` and `x`, has more terms on `system` than
`a` and `x` together; counted on the system conserving nothing when one of them has no definite
charge, as a merge may have on its way to a product that has one
"""
function more_terms(system::System, m, a, x)
    cm, ca, cx = or_on_plain(s -> (term_count(s, m), term_count(s, a), term_count(s, x)), system)
    return cm > ca * cx
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
            if m isa AtIndex && more_terms(system, m, a, x)
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
    as_product(a)

the placed operator `a`, a placed tensor product being taken as the product of its factors
placed, `(A ⊗ B)(i, j) == A(i) * B(j)`
"""
as_product(a::AtIndex{IsOperator}) =
    if a.op isa TensorOp
        ProdOp(IndexedOp{IsOperator}[ as_product(o(a.index[r]...))
                                      for (o, r) in zip(a.op.subs, factor_sites(a.op)) ])
    else
        a
    end
as_product(a) = a

"""
    gathered(system, op)

the placed operator `op` rewritten without developing anything, into a normal form: in each
product, a tensor product is taken as the product of its factors, see `as_product`, which then
merge one by one rather than through the dense matrix of the whole, two factors sharing a site
are merged where they can be, see `merged_pass`, until nothing merges, then the factors are
sorted by their sites as far as they (anti)commute, see `sort_factors!`. Sums are gathered term
by term, a sum on a single site becoming the operator of its matrix, as a merged product does,
and an operator in normal form is left as it is.
"""
function gathered(system::System, a::SumOp{IsOperator, Indexed})
    s = SumOp([ gathered(system, x) for x in a.subs ])
    sites = sites_of(s)
    if !(s isa SumOp) || length(sites) ≠ 1
        return s
    end
    ms = [ matrix_on(system, x, sites) for x in s.subs ]
    return placed_matrix(sum(ms), sum(norm, ms), "($s)", sites)
end
gathered(system::System, a::ScalarOp{IsOperator, Indexed}) =
    scale(a.coef, gathered(system, a.arg))
gathered(::System, a::Union{AtIndex{IsOperator}, IdentityOp{IsOperator, Indexed}}) = a

function gathered(system::System, a::ProdOp{IsOperator, Indexed})
    p = ProdOp(IndexedOp{IsOperator}[ gathered(system, as_product(x)) for x in a.subs ])
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
    return lower(system, a, gathered(system, a))
end

# from `g`, the operator `a` gathered
function lower(system::System, a::IndexedOp{IsOperator}, g::IndexedOp{IsOperator})
    p = scalararg(g)
    if is_charged(system) && p isa ProdOp && allunique(reduce(vcat, map(sites_of, p.subs)))
        # a product of factors on distinct sites cancels nothing: a factor changing the charges
        # by different amounts is refused before the product multiplies its terms
        parts = [ check_flux(Pure, system, a, reduced(op_terms(system, x, ()))) for x in p.subs ]
        terms = reduced(scaled(scalarcoef(g), tree_products(parts)))
    else
        terms = reduced(op_terms(system, g, ()))
    end
    return check_flux(Pure, system, a, terms)
end

function lower(system::System, a::IndexedOp{<:Union{IsChannel, IsGenerator}})
    check_indices(system, a)
    return check_flux(Mixed, system, a, reduced(op_terms(system, a, ())))
end
