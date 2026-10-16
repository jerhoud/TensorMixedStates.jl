# The operators compacted on a system, see `CompactedOp`, in which an operator is measured and
# laid as an MPO, and their algebra: the adjoint, the product, the linear combinations and the
# powers, each reduced to as few channels as its rank, see `reduced_block`. An operator is
# compacted from its terms, but for the products of sums that would develop into many terms,
# made of the compactions of their factors multiplied, and so are the gates and the dissipators
# of operators of many terms, from the sides of these.

"""
    term_size(t)

the coefficient of the term `t` times the operator norms of its blocks, which bounds its norm,
the strings of `F` changing none, and so the expectation values of the term on a normalized state
"""
term_size(t::Term) = abs(t.coef) * prod(opnorm(b.m) for b in t.blocks; init = 1.0)

"""
    term_width(system, t)

a bound of the width of the spectrum of ``(t + t^\\dagger)/2`` for the term `t` of an operator on
a pure state of `system`: the bounds of its terms add up to one of a Hermitian operator, by the
inequalities of Weyl. It is exact for a term of at most two sites, as `apply` takes them dense,
and twice the norm of the term beyond; with a Jordan-Wigner string, it is twice the largest
eigenvalue, in absolute value, of the Hermitian part of the term with the `F` of the string on
the site of its first block, the `F` of the sites between, ±1, changing nothing; a constant
has none
"""
function term_width(system::System, t::Term)
    if isempty(t.blocks)
        return 0.
    elseif length(t.blocks) > 2
        return 2abs(t.coef) * prod(opnorm(b.m) for b in t.blocks)
    end
    m = t.coef * foldl(kron, [ on_string(t.kets, b.site) ? b.m * matrix(F, system[b.site]) : b.m
                               for b in t.blocks ])
    λ = eigvals(Hermitian((m + m') / 2))
    return isempty(t.kets) ? last(λ) - first(λ) : 2maximum(abs, λ)
end

"""
    terms_flux(R, system, a, terms)

the flux of the terms of the operator `a`, all the same once checked, see `check_flux`:
`nothing` when there is none, and no charge on a system conserving nothing
"""
terms_flux(::Type{R}, system::System, a, terms::Vector{Term}) where R =
    isempty(terms) ? nothing : is_charged(system) ? term_flux(R, system, a, first(terms)) : QN()

"""
    CompactedOp{R}

an operator compacted on `system`, see `compacted`, a superoperator for `R = Mixed`: its
constant, its terms of one site, a matrix for each site, and its block of channels, the bounds of
its norm and of the width of its spectrum, see `term_size` and `term_width`, its flux, `nothing`
when it has no term, and its terms when it was made of them alone. Each time function of a
`PreMPO` holds one.
"""
struct CompactedOp{R}
    system::System
    constant::Number
    singles::Dict{Int, Matrix}
    block::Union{Nothing, ComBlock}
    norm::Float64
    width::Float64
    flux::Union{Nothing, QN}
    terms::Union{Nothing, Vector{Term}}
end

"""
    op_eltype(c)

the type of the values of the operator compacted `c`, `Float64` at least: real for a real
operator, whose MPO is real, cheaper to contract
"""
op_eltype(c::CompactedOp) =
    promote_type(Float64, typeof(c.constant), (eltype(m) for m in values(c.singles))...,
                 isnothing(c.block) ? Float64 : eltype(c.block))

"""
    zero_op(R, system)

the operator zero compacted on `system`, of no term
"""
zero_op(::Type{R}, system::System) where R =
    CompactedOp{R}(system, 0, Dict{Int, Matrix}(), nothing, 0.0, 0.0, nothing, nothing)

"""
    from_terms(R, system, a, terms, limits = Limits())

the operator `a` compacted on `system` from its `terms`, see `compacted`, truncated by `limits`,
the bounds of its norm and of its width the sums of those of its terms, the width that of a pure
state alone, `a` naming the operator in the refusals
"""
function from_terms(::Type{R}, system::System, a, terms::Vector{Term},
                    limits::Limits = Limits()) where R
    width = R === Pure ? sum(t -> term_width(system, t), terms; init = 0.0) : 0.0
    return CompactedOp{R}(system, compacted(R, system, terms, limits)...,
                          sum(term_size, terms; init = 0.0), width,
                          terms_flux(R, system, a, terms), terms)
end

"""
    channel_count(c, k)

the number of channels of the block of `c` on the link on the right of site `k`, none outside
the block
"""
function channel_count(c::CompactedOp, k::Int)
    b = c.block
    if isnothing(b) || !(b.start ≤ k < b.start + length(b.pieces) - 1)
        return 0
    end
    return b.dims[k - b.start + 1]
end

"""
    site_moves(c, k)

the moves of the operator compacted `c` on site `k`, as their matrices, between the states of
the links on either side: a channel, 0 the term not yet begun or -1 the term finished, the
identity going from 0 to 0 and from -1 to -1, the constant a term of one site on the first site
"""
function site_moves(c::CompactedOp{R}, k::Int) where R
    d = block_dim(R, c.system[k])
    T = op_eltype(c)
    id = Matrix{T}(I, d, d)
    moves = Dict{Tuple{Int, Int}, Matrix{T}}((0, 0) => id, (-1, -1) => id)
    if k == 1 && !iszero(c.constant)
        add!(moves, (0, -1), c.constant * id)
    end
    if haskey(c.singles, k)
        add!(moves, (0, -1), Matrix{T}(c.singles[k]))
    end
    b = c.block
    if !isnothing(b) && b.start ≤ k < b.start + length(b.pieces)
        j = k - b.start + 1
        for ((l, r), comb) in b.pieces[j]
            add!(moves, (l == 0 ? 0 : l, r == 0 ? -1 : r),
                 Matrix{T}(combination_matrix(b, j, comb, d)))
        end
    end
    return moves
end

"""
    adjoint(c)

the adjoint of the operator compacted `c`, its pieces taken adjoint, of the opposite flux
"""
function Base.adjoint(c::CompactedOp{R}) where R
    b = c.block
    block = isnothing(b) ? nothing :
        ComBlock(b.start, copy(b.dims),
                 [ typeof(p)(lr => typeof(comb)(a => conj(x) for (a, x) in comb)
                             for (lr, comb) in p) for p in b.pieces ],
                 [ [ collect(m') for m in atoms ] for atoms in b.atoms ])
    return CompactedOp{R}(c.system, conj(c.constant),
                          Dict{Int, Matrix}(k => collect(m') for (k, m) in c.singles), block,
                          c.norm, c.width, isnothing(c.flux) ? nothing : -c.flux, nothing)
end

"""
    from_moves(R, system, T, dims, moves, limits)

the terms of one site and the block of channels of the operator compacted on `system` whose
pieces on each site `k` are `moves(k)`, as `(l, r, m, x)`: `x` times the matrix `m` from channel
`l` of the link on the left to channel `r` of the link on the right, `dims` their numbers, 0 for
an opening or a closing, a term of one site when both are; the matrices of each site numbered as
its atoms, the block reduced, see `reduced_block`, and truncated by `limits`
"""
function from_moves(::Type{R}, sys::System, ::Type{T}, dims::Vector{Int}, moves,
                    limits::Limits) where {R, T}
    n = length(sys)
    singles = Dict{Int, Matrix{T}}()
    pieces = [ SitePieces{T}() for _ in 1:n ]
    atoms = [ Matrix{T}[] for _ in 1:n ]
    for k in 1:n
        ids = Dict{Matrix{T}, Int}()
        for (l, r, m, x) in moves(k)
            if l == 0 && r == 0
                add!(singles, k, x * m)
                continue
            end
            a = get!(ids, m) do
                push!(atoms[k], m)
                return length(atoms[k])
            end
            add!(get!(Combination{T}, pieces[k], (l, r)), a, T(x))
        end
    end
    block = n == 1 ? nothing :
            reduced_block(R, sys, ComBlock{T}(1, dims, pieces, atoms), singles, limits)
    return singles, block
end

"""
    product(a, b, limits = Limits())

the product of the operators compacted `a` and `b`, of one system: the channels of a link the
pairs of its states in `a` and in `b`, see `site_moves`, the pieces the products of theirs, the
block reduced, see `reduced_block`, to as few channels as the rank of the product across each
link and truncated by `limits`, at a cost linear in the number of sites where the terms of the
product grow as the product of theirs. The bound of its norm is the product of theirs, and that
of the width of its spectrum twice that, the width of the Hermitian part of an operator being at
most twice its norm.
"""
function product(a::CompactedOp{R}, b::CompactedOp{R}, limits::Limits = Limits()) where R
    sys = a.system
    if isnothing(a.flux) || isnothing(b.flux)
        return zero_op(R, sys)
    end
    n = length(sys)
    T = promote_type(op_eltype(a), op_eltype(b))
    ends = ((0, 0), (-1, -1))
    channels = map(1:n-1) do k
        sa, sb = [0; 1:channel_count(a, k); -1], [0; 1:channel_count(b, k); -1]
        pairs = [ (x, y) for x in sa for y in sb if (x, y) ∉ ends ]
        return Dict(p => i for (i, p) in enumerate(pairs))
    end
    function moves(k)
        out = Tuple{Int, Int, Matrix{T}, T}[]
        for ((a1, b1), m1) in site_moves(a, k), ((a2, b2), m2) in site_moves(b, k)
            from, to = (a1, a2), (b1, b2)
            # nothing comes from the left of the first site, nothing open goes past the last
            if (from == to && from in ends) || (k == 1 && from ≠ (0, 0)) ||
               (k == n && to ≠ (-1, -1))
                continue
            end
            m = m1 * m2
            # a product vanishing but for rounding, which would be an atom of its own
            if norm(m) > rounding_tol * norm(m1) * norm(m2)
                push!(out, (from == (0, 0) ? 0 : channels[k-1][from],
                            to == (-1, -1) ? 0 : channels[k][to], m, one(T)))
            end
        end
        return out
    end
    singles, block = from_moves(R, sys, T, [ length(c) for c in channels ], moves, limits)
    bound = a.norm * b.norm
    return CompactedOp{R}(sys, zero(T), singles, block, bound, 2bound, a.flux + b.flux, nothing)
end

"""
    combined(system, parts, limits = Limits(), what = "the sum")

the sum of the operators compacted on `system` of `parts`, pairs of an operator and its
coefficient: the channels of a link those of all of them side by side, the coefficient of each
on its openings, the block reduced, see `reduced_block`, and truncated by `limits`. The bound of
its norm is the sum of theirs times the coefficients, and that of the width of its spectrum the
sum of theirs, by the inequalities of Weyl, a complex coefficient `c` giving `2|c|` times the
norm. Refused when they change the charges by different amounts, see `flux_refusal`, `what`
naming the operator.
"""
function combined(sys::System, parts::Vector{<:Tuple{CompactedOp{R}, Number}},
                  limits::Limits = Limits(), what = "the sum") where R
    kept = [ (c, x) for (c, x) in parts if !iszero(x) && !isnothing(c.flux) ]
    if isempty(kept)
        return zero_op(R, sys)
    end
    fluxes = unique(c.flux for (c, _) in kept)
    if length(fluxes) > 1
        flux_refusal(sys, what, fluxes)
    end
    n = length(sys)
    T = promote_type((op_eltype(c) for (c, _) in kept)..., (typeof(x) for (_, x) in kept)...)
    # the channels of the i-th operator on link k numbered after those of the ones before
    offsets = [ cumsum([ 0; [ channel_count(c, k) for (c, _) in kept ] ]) for k in 1:n-1 ]
    # the moves of each operator but the identity, its coefficient on its openings
    moves(k) = [ (from == 0 ? 0 : offsets[k-1][i] + from, to == -1 ? 0 : offsets[k][i] + to, m,
                  from == 0 ? T(x) : one(T))
                 for (i, (c, x)) in enumerate(kept) for ((from, to), m) in site_moves(c, k)
                 if !(from == to && from in (0, -1)) ]
    singles, block = from_moves(R, sys, T, last.(offsets), moves, limits)
    return CompactedOp{R}(sys, zero(T), singles, block, sum(abs(x) * c.norm for (c, x) in kept),
                          sum(isreal(x) ? abs(x) * c.width : 2abs(x) * c.norm for (c, x) in kept),
                          only(fluxes), nothing)
end

"""
    power(c, p)

the `p`-th power of the operator compacted `c` by squaring, the power 4 being `(c c)(c c)`, two
products, and the power 0 the identity
"""
function power(c::CompactedOp{R}, p::Integer) where R
    if p < 0
        error("a PreMPO has no power $p: it has no inverse")
    elseif p == 0
        return CompactedOp{R}(c.system, 1.0, Dict{Int, Matrix}(), nothing, 1.0, 0.0, QN(), nothing)
    elseif p == 1
        return c
    end
    h = power(c, p ÷ 2)
    s = product(h, h)
    return isodd(p) ? product(s, c) : s
end

"""
    vanishes(c)

whether the operator compacted `c` is zero but for rounding: no term, or no channel and terms
of one site of rounding size next to the bound of its norm
"""
vanishes(c::CompactedOp) = isnothing(c.flux) ||
    ((isnothing(c.block) || all(iszero, c.block.dims)) &&
     abs(c.constant) ≤ rounding_tol * c.norm &&
     all(m -> normalized_norm(m) ≤ rounding_tol * c.norm, values(c.singles)))

"""
    term_estimate(a)

the number of terms of the placed operator `a` developed as it is written, a bound, terms that
merge counted apart
"""
term_estimate(a::SumOp) = sum(term_estimate, a.subs)
term_estimate(a::ProdOp) = prod(term_estimate, a.subs)
term_estimate(a::ScalarOp) = term_estimate(a.arg)
term_estimate(::Op) = 1

"""
    by_blocks(a)

whether the placed operator `a`, gathered, see `gathered`, is a product whose terms would
outnumber those of its factors together, see `term_estimate`, so that it is better made of the
compactions of its factors multiplied, see `compacted_op`
"""
function by_blocks(a::IndexedOp)
    p = scalararg(a)
    if !(p isa ProdOp)
        return false
    end
    counts = term_estimate.(p.subs)
    return prod(counts) > sum(counts)
end

"""
    addends(a)

the operators the placed operator `a` is the sum of, its coefficient left out, `a` itself when it
is no sum
"""
function addends(a::IndexedOp)
    s = scalararg(a)
    return s isa SumOp ? s.subs : [s]
end

"""
    blocked_product(fs, factor)

the product of the operators `fs`, `factor(f)` giving the compaction of the factor `f`, equal
factors in a row raised by squaring, see `power`
"""
function blocked_product(fs, factor)
    c = nothing
    k = 1
    while k ≤ length(fs)
        run = findfirst(f -> f ≠ fs[k], fs[k:end])
        n = isnothing(run) ? length(fs) - k + 1 : run - 1
        p = power(factor(fs[k]), n)
        c = isnothing(c) ? p : product(c, p)
        k += n
    end
    return c
end

"""
    add_part!(parts, c, coef)

add to `parts`, to be summed by `combined`, the operator compacted `c` times `coef`, unless it
vanishes, see `vanishes`: a product of definite fluxes that vanishes has none
"""
function add_part!(parts::Vector, c::CompactedOp, coef::Number)
    if !vanishes(c)
        push!(parts, (c, coef))
    end
    return parts
end

"""
    compacted_op(system, op, limits = Limits(); side = nothing)

the operator `op`, placed on sites, compacted on `system` and truncated by `limits`, see
`truncated_block`: from its terms, see `lower`, but for its products that would develop into
more terms than their factors have, see `by_blocks`. Those are made of the compactions of their
factors, made the same way, multiplied, equal factors in a row by squaring, see `power`, and
added to that of the rest in one reduction, see `combined`. With `side`, `:left` or `:right`,
the superoperator ``\\rho \\mapsto op\\,\\rho`` or ``\\rho \\mapsto \\rho\\,op^\\dagger``, see
`side_terms`, a product made of the sides of its factors in the same order on both sides,
``\\rho (AB)^\\dagger`` being ``(\\rho B^\\dagger) A^\\dagger``. Refused when it changes the
charges by no definite amount, see `FluxError`.
"""
function compacted_op(system::System, a::IndexedOp{IsOperator}, limits::Limits = Limits();
                      side::Union{Nothing, Symbol} = nothing)
    check_indices(system, a)
    g = gathered(system, a)
    summands = addends(g)
    R = isnothing(side) ? Pure : Mixed
    if !any(by_blocks, summands)
        terms = lower(system, a, g)
        if !isnothing(side)
            terms = check_flux(Mixed, system, a, side_terms(side, terms))
        end
        return from_terms(R, system, a, terms, limits)
    end
    # the side of the right takes the adjoint of the operator, and so of its coefficients
    coef(x) = side == :right ? conj(x) : x
    parts = Tuple{CompactedOp{R}, Number}[]
    rest = filter(x -> !by_blocks(x), summands)
    if !isempty(rest)
        add_part!(parts, compacted_op(system, SumOp(rest); side), coef(scalarcoef(g)))
    end
    for x in filter(by_blocks, summands)
        factors = blocked_product(prodsubs(scalararg(x)), f -> compacted_op(system, f; side))
        add_part!(parts, factors, coef(scalarcoef(g) * scalarcoef(x)))
    end
    return combined(system, parts, limits, a)
end

"""
    many_terms(k)

whether the placed operator `k` has more than two terms, or a product that would develop into
more terms than its factors have, see `by_blocks`, so that a gate or a dissipator of `k` is
better made of its sides, see `compacted_super`
"""
many_terms(k::IndexedOp) = term_estimate(k) > 2 || any(by_blocks, addends(k))

"""
    super_by_blocks(a)

whether the superoperator `a`, placed on sites, is better made of the sides of its operators,
see `compacted_super`: a gate or a dissipator of an operator of many terms, see `many_terms`, an
evolver of one with a product that would develop, or a product of channels whose terms would
outnumber those of its factors together, a gate of `k` having the square of the terms of `k`
"""
function super_by_blocks(a::IndexedOp)
    s = scalararg(a)
    if s isa Union{Gate{Indexed}, Dissipator{Indexed}}
        return s.arg isa IndexedOp{IsOperator} && many_terms(s.arg)
    elseif s isa Evolver{Indexed}
        return any(by_blocks, addends(s.arg))
    elseif s isa ProdOp
        counts = [ scalararg(x) isa Gate{Indexed} ? term_estimate(scalararg(x).arg)^2 : 1
                   for x in s.subs ]
        return prod(counts) > sum(counts)
    end
    return false
end

"""
    compacted_super(system, a, limits = Limits())

the superoperator `a`, a channel or a generator placed on sites, compacted on `system` and
truncated by `limits`, see `truncated_block`: from its terms, see `lower`, but for what is better
made of the sides of its operators, see `super_by_blocks` and `compacted_op`, multiplied, the
adjoint of a side being the side of the adjoint. The evolver of `K` is ``K\\rho + \\rho
K^\\dagger``, its gate ``K\\rho K^\\dagger``, the dissipator of `L` ``L\\rho L^\\dagger -
(L^\\dagger L\\rho + \\rho L^\\dagger L)/2``, and a product of channels that of their compactions.
"""
function compacted_super(system::System, a::IndexedOp, limits::Limits = Limits())
    check_indices(system, a)
    c = scalarcoef(a)
    summands = addends(a)
    if !any(super_by_blocks, summands)
        return from_terms(Mixed, system, a, lower(system, a), limits)
    end
    parts = Tuple{CompactedOp{Mixed}, Number}[]
    rest = filter(x -> !super_by_blocks(x), summands)
    if !isempty(rest)
        add_part!(parts, compacted_super(system, SumOp(rest)), c)
    end
    for x in filter(super_by_blocks, summands)
        y, k = scalararg(x), c * scalarcoef(x)
        if y isa ProdOp
            add_part!(parts, blocked_product(y.subs, f -> compacted_super(system, f)), k)
            continue
        end
        l, r = compacted_op(system, y.arg; side = :left), compacted_op(system, y.arg; side = :right)
        if y isa Evolver
            add_part!(add_part!(parts, l, k), r, k)
            continue
        end
        add_part!(parts, product(l, r), k)
        if y isa Dissipator
            add_part!(parts, product(l', l), -k / 2)
            add_part!(parts, product(r', r), -k / 2)
        end
    end
    return combined(system, parts, limits, a)
end
