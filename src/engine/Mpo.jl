# The MPO of an operator on the system of a state: PreMPO, which lays its terms, compacted, site
# by site, and make_mpo, which builds the MPO from it. On a mixed state, the MPO of a channel or
# a generator, its superoperator laid on the mixed indices.

export PreMPO, make_mpo

"""
    mixed_elements(system, i)

the element of the density matrix of site `i` of `system` at each position of its mixed index,
numbered as in the matrices of superoperators, the ket faster
"""
function mixed_elements(system::System, i::Int)
    k = siteind(Mixed, system, i)
    e = ket_bras(system, i)
    d = dim(system[i])
    return invperm(vec([ findfirst(!iszero, Array(e(m, n), k)) for m in 1:d, n in 1:d ]))
end

"""
    block_orders(R, system)

the order on the index of each site of the elements of the matrix of an operator there: none
to take for `Pure`, see `mixed_elements` for `Mixed`, read once for each kind of site, which
alone gives it
"""
block_orders(::Type{Pure}, ::System) = Vector{Int}[]

function block_orders(::Type{Mixed}, system::System)
    orders = Dict{AbstractSite, Vector{Int}}()
    return [ get!(() -> mixed_elements(system, i), orders, system[i]) for i in 1:length(system) ]
end

"""
    PreMPO(state, op)

the operator `op`, placed on sites, laid on the system of `state` as the pieces of its MPO,
its terms compacted into as few channels as a triangular MPO can have; `make_mpo` builds the
MPO from it. A vector of operators gives each its own time function, see `make_mpo`. On a pure
state, a generator with no dissipator, ``\\rho \\mapsto K\\rho + \\rho K^\\dagger``, gives
the MPO of `K`, which the state follows as ``d\\psi/dt = K\\psi``. On a mixed state, a channel
or a generator gives the MPO of its superoperator, and an operator `K` that of `Evolver(K)`.

`pre'`, `a * b`, `a ^ p`, `a + b`, `a - b`, `c * a` and `sum(v)`, of `PreMPO`s of one system
each of a single time function, are the `PreMPO`s of the adjoint, of the product, of the power
and of the linear combinations of their operators: built with as few channels as their rank
across each link, at a cost linear in the number of sites, where the terms of a product grow as
the product of theirs. `PreMPO(state, op)` builds so the products of sums of `op` that would
develop into more terms than their factors have.
"""
struct PreMPO{R <: PM}
    system::System
    linkdims::Vector{Int}
    # the pieces on each site, see `add_term!`
    terms::Vector{Vector{Tuple{Int, Int, ITensor, Int, Matrix}}}
    # the number of time functions the terms take, one per operator of a vector: it cannot be
    # counted from `terms`, a whole operator vanishing on its sites
    nterms::Int
    orders::Vector{Vector{Int}}
    # the last MPOs built for each method of the solvers, with what they were built for, see
    # `cached_mpos`
    cache::Dict{Symbol, Tuple{Any, Any}}
    lock::ReentrantLock
    # for a pure state, a bound of the width of the spectrum of the operator, see `term_width`
    width::Base.RefValue{Float64}
    # a bound of the norm of the operator, see `term_size`, which a product multiplies
    norm::Base.RefValue{Float64}
    function PreMPO{R}(system::System, nterms::Int = 1) where R
        n = length(system)
        return new{R}(system, fill(1, n - 1),
                      [ Tuple{Int, Int, ITensor, Int, Matrix}[] for _ in 1:n ], nterms,
                      block_orders(R, system), Dict{Symbol, Tuple{Any, Any}}(), ReentrantLock(),
                      Ref(0.), Ref(0.))
    end
end

"""
    on_index(pre, k, m)

the matrix `m` of an operator on site `k`, its rows and columns in the order of the index of
the site, see `block_orders`
"""
on_index(::PreMPO{Pure}, ::Int, m::Matrix) = m
on_index(pre::PreMPO{Mixed}, k::Int, m::Matrix) = m[pre.orders[k], pre.orders[k]]

"""
    add_term!(pre, k, l, r, m, ref)

add to site `k` of `pre`, and return it, the piece of matrix `m` going from channel `l` to
channel `r` and taking the time function `ref`, with its tensor on the index of the site, both
in the order of that index
"""
function add_term!(pre::PreMPO{R}, k::Int, l::Int, r::Int, m::Matrix, ref::Int) where R
    idx = siteind(R, pre.system, k)
    m = on_index(pre, k, m)
    push!(pre.terms[k], (l, r, charged_itensor(m, [idx', dag(idx)]), ref, m))
    return pre
end

"""
    check_coefs(pre, coefs)

`coefs` as a `Vector{Float64}`, refused unless it holds one real value per time function of
`pre`
"""
function check_coefs(pre::PreMPO, coefs::AbstractVector)
    if length(coefs) ≠ pre.nterms
        error("an evolver of $(pre.nterms) terms takes as many time functions, got $(length(coefs))")
    end
    for (k, c) in enumerate(coefs)
        if !(c isa Real)
            error("time functions take real values, and term $k got $c: " *
                  "write f * A as real(f) * A + imag(f) * (im * A)")
        end
    end
    return Vector{Float64}(coefs)
end

"""
    term_size(t)

the coefficient of the term `t` times the operator norms of its blocks, which bounds its norm,
the strings of `F` changing none, and so the expectation values of the term on a normalized state
"""
term_size(t::Term) = abs(t.coef) * prod(opnorm(m) for (_, m) in t.blocks; init = 1.0)

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
        return 2abs(t.coef) * prod(opnorm(m) for (_, m) in t.blocks)
    end
    m = t.coef * foldl(kron, [ on_string(t.kets, k) ? m * matrix(F, system[k]) : m
                               for (k, m) in t.blocks ])
    λ = eigvals(Hermitian((m + m') / 2))
    return isempty(t.kets) ? last(λ) - first(λ) : 2maximum(abs, λ)
end

"""
    lay_compacted!(pre, constant, singles, b, ref)

add to `pre`, and return it, an operator compacted, see `compacted`, taking the time function
`ref`: the constant on the first site of the others, so that the MPO is the identity outside
their sites, see `mpo_window`, the terms of one site `singles` on their site, and the channels
of the block `b` numbered after those their links already have
"""
function lay_compacted!(pre::PreMPO{R}, constant::Number, singles::AbstractDict,
                        b::Union{Nothing, ComBlock}, ref::Int) where R
    sys = pre.system
    if !iszero(constant)
        ks = [ collect(keys(singles)) ; isnothing(b) ? Int[] : [b.start] ]
        k = isempty(ks) ? 1 : minimum(ks)
        d = block_dim(R, sys[k])
        add_term!(pre, k, 1, 1, Matrix(constant * I, d, d), ref)
    end
    for k in sort!(collect(keys(singles)))
        add_term!(pre, k, 1, 1, singles[k], ref)
    end
    if isnothing(b)
        return pre
    end
    ks = b.start .+ (0:length(b.pieces) - 1)
    ld = pre.linkdims
    offs = [ ld[k] for k in ks[1:end-1] ]
    for (j, k) in enumerate(ks[1:end-1])
        ld[k] += b.dims[j]
    end
    for (j, k) in enumerate(ks), ((l, r), comb) in sort!(collect(b.pieces[j]); by = first)
        m = combination_matrix(b, j, comb, block_dim(R, sys[k]))
        if !iszero(m)
            add_term!(pre, k, l == 0 ? 1 : offs[j-1] + l, r == 0 ? 1 : offs[j] + r, m, ref)
        end
    end
    return pre
end

"""
    lay_terms!(pre, terms, ref)

add to `pre`, and return it, the terms of an operator, compacted, see `compacted` and
`lay_compacted!`, taking the time function `ref`
"""
function lay_terms!(pre::PreMPO{R}, terms::Vector{Term}, ref::Int) where R
    if R === Pure
        pre.width[] += sum(t -> term_width(pre.system, t), terms; init = 0.)
    end
    pre.norm[] += sum(term_size, terms; init = 0.)
    return lay_compacted!(pre, compacted(R, pre.system, terms)..., ref)
end

"""
    natural(pre, k, m)

the matrix `m` of site `k` of `pre`, in the order of the index of the site, back in the order of
the site, see `on_index`
"""
natural(::PreMPO{Pure}, ::Int, m::Matrix) = m

function natural(pre::PreMPO{Mixed}, k::Int, m::Matrix)
    n = similar(m)
    n[pre.orders[k], pre.orders[k]] = m
    return n
end

"""
    site_moves(pre, k)

the pieces of site `k` of `pre`, of a single time function, as the matrices, in the order of the
site, of the moves between the states of the links on either side: a channel, 0 the term not
yet begun or -1 the term finished, the identity going from 0 to 0 and from -1 to -1, see
`window_mpo`
"""
function site_moves(pre::PreMPO{R}, k::Int) where R
    d = block_dim(R, pre.system[k])
    id = Matrix{ComplexF64}(I, d, d)
    moves = Dict{Tuple{Int, Int}, Matrix{ComplexF64}}((0, 0) => id, (-1, -1) => id)
    for (l, r, _, _, m) in pre.terms[k]
        add!(moves, (l == 1 ? 0 : l, r == 1 ? -1 : r), natural(pre, k, m))
    end
    return moves
end

"""
    adjoint(pre)

the `PreMPO` of the adjoint of the operator of `pre`, its pieces taken adjoint
"""
function Base.adjoint(pre::PreMPO{R}) where R
    out = PreMPO{R}(pre.system, pre.nterms)
    out.linkdims .= pre.linkdims
    for k in eachindex(pre.terms), (l, r, _, ref, m) in pre.terms[k]
        add_term!(out, k, l, r, collect(natural(pre, k, m)'), ref)
    end
    out.width[] = pre.width[]
    out.norm[] = pre.norm[]
    return out
end

"""
    a * b

the `PreMPO` of the product of the operators of `a` and `b`, of one system and representation
and of a single time function each: the channels of a link the pairs of its states in `a` and
in `b`, see `site_moves`, the pieces the products of theirs, the block reduced as that of the
terms would be, see `reduced_block`, to as few channels as the rank of the product across each
link, at a cost linear in the number of sites where the terms of the product grow as the
product of theirs. The bound of its norm is the product of theirs, and that of the width of its
spectrum, which a search of excited states takes for its weight, twice that, the width of the
Hermitian part of an operator being at most twice its norm.
"""
function (a::PreMPO{R} * b::PreMPO{R}) where {R <: PM}
    if a.system ≠ b.system
        error("the PreMPOs of a product are on two systems")
    elseif a.nterms ≠ 1 || b.nterms ≠ 1
        error("a PreMPO of several time functions is not multiplied")
    end
    sys = a.system
    n = length(sys)
    T = ComplexF64
    ends = ((0, 0), (-1, -1))
    channels = map(1:n-1) do k
        sa, sb = [0; 2:a.linkdims[k]; -1], [0; 2:b.linkdims[k]; -1]
        pairs = [ (x, y) for x in sa for y in sb if (x, y) ∉ ends ]
        return Dict(p => i for (i, p) in enumerate(pairs))
    end
    singles = Dict{Int, Matrix{T}}()
    pieces = [ SitePieces{T}() for _ in 1:n ]
    atoms = [ Matrix{T}[] for _ in 1:n ]
    for k in 1:n
        ids = Dict{Matrix{T}, Int}()
        for ((a1, b1), m1) in site_moves(a, k), ((a2, b2), m2) in site_moves(b, k)
            from, to = (a1, a2), (b1, b2)
            # nothing comes from the left of the first site, nothing open goes past the last
            if (from == to && from in ends) || (k == 1 && from ≠ (0, 0)) || (k == n && to ≠ (-1, -1))
                continue
            end
            m = m1 * m2
            # a product vanishing but for rounding, which would be an atom of its own
            if norm(m) ≤ rounding_tol * norm(m1) * norm(m2)
                continue
            elseif from == (0, 0) && to == (-1, -1)
                add!(singles, k, m)
                continue
            end
            x = get!(ids, m) do
                push!(atoms[k], m)
                return length(atoms[k])
            end
            l = from == (0, 0) ? 0 : channels[k-1][from]
            r = to == (-1, -1) ? 0 : channels[k][to]
            add!(get!(Combination{T}, pieces[k], (l, r)), x, one(T))
        end
    end
    block = n == 1 ? nothing :
            reduced_block(R, sys, ComBlock{T}(1, [ length(c) for c in channels ], pieces, atoms), singles)
    out = lay_compacted!(PreMPO{R}(sys), zero(T), singles, block, 1)
    out.norm[] = a.norm[] * b.norm[]
    out.width[] = 2out.norm[]
    return out
end

(::PreMPO * ::PreMPO) = error("the PreMPOs of a product are of two representations")

"""
    combined(pres, coefs)

the terms of one site, a matrix for each site, and the block of channels, reduced, see
`reduced_block`, of the sum of the operators of the `PreMPO`s `pres`, of a single time function
each, times `coefs`: the channels of a link those of all of them side by side, a coefficient on
the openings of its `PreMPO`
"""
function combined(pres::Vector{<:PreMPO{R}}, coefs::Vector{<:Number}) where R
    sys = first(pres).system
    n = length(sys)
    T = ComplexF64
    # the channels of the i-th PreMPO on link k numbered after those of the ones before
    offsets = [ cumsum([ 0; [ p.linkdims[k] - 1 for p in pres ] ]) for k in 1:n-1 ]
    singles = Dict{Int, Matrix{T}}()
    pieces = [ SitePieces{T}() for _ in 1:n ]
    atoms = [ Matrix{T}[] for _ in 1:n ]
    for k in 1:n
        ids = Dict{Matrix{T}, Int}()
        for (i, (p, c)) in enumerate(zip(pres, coefs)), ((from, to), m) in site_moves(p, k)
            if iszero(c) || (from == to && from in (0, -1))
                continue
            elseif from == 0 && to == -1
                add!(singles, k, c * m)
                continue
            end
            x = get!(ids, m) do
                push!(atoms[k], m)
                return length(atoms[k])
            end
            l = from == 0 ? 0 : offsets[k-1][i] + from - 1
            r = to == -1 ? 0 : offsets[k][i] + to - 1
            add!(get!(Combination{T}, pieces[k], (l, r)), x, from == 0 ? T(c) : one(T))
        end
    end
    block = n == 1 ? nothing :
            reduced_block(R, sys, ComBlock{T}(1, last.(offsets), pieces, atoms), singles)
    return singles, block
end

"""
    linear_combination(pres, coefs)

the `PreMPO` of the sum of the operators of `pres` times `coefs`, see `combined`: the bound of
its norm is the sum of theirs times the coefficients, and that of the width of its spectrum the
sum of theirs, by the inequalities of Weyl, a complex coefficient `c` giving `2|c|` times the
norm, the width of the Hermitian part of an operator being at most twice its norm
"""
function linear_combination(pres::Vector{<:PreMPO{R}}, coefs::Vector{<:Number}) where R
    sys = first(pres).system
    if any(p -> p.system ≠ sys, pres)
        error("the PreMPOs of a sum are on two systems")
    elseif any(p -> p.nterms ≠ 1, pres)
        error("a PreMPO of several time functions is not added or scaled")
    end
    out = lay_compacted!(PreMPO{R}(sys), zero(ComplexF64), combined(pres, coefs)..., 1)
    out.norm[] = sum(abs(c) * p.norm[] for (p, c) in zip(pres, coefs))
    out.width[] = sum(isreal(c) ? abs(c) * p.width[] : 2abs(c) * p.norm[] for (p, c) in zip(pres, coefs))
    return out
end

(a::PreMPO{R} + b::PreMPO{R}) where {R <: PM} = linear_combination([a, b], [1, 1])
(a::PreMPO{R} - b::PreMPO{R}) where {R <: PM} = linear_combination([a, b], [1, -1])
(::PreMPO + ::PreMPO) = error("the PreMPOs of a sum are of two representations")
(::PreMPO - ::PreMPO) = error("the PreMPOs of a sum are of two representations")
(c::Number * a::PreMPO) = linear_combination([a], [c])
(a::PreMPO * c::Number) = c * a
(a::PreMPO / c::Number) = inv(c) * a
-(a::PreMPO) = -1 * a
# reduced once, where `+` would reduce each partial sum
Base.sum(pres::AbstractVector{<:PreMPO{R}}) where {R <: PM} =
    linear_combination(collect(pres), ones(length(pres)))

"""
    a ^ p

the `PreMPO` of the `p`-th power of the operator of `a` by squaring, `a^4` being `(a * a)^2`, two
products, and `a^0` the identity
"""
function (a::PreMPO{R} ^ p::Integer) where {R <: PM}
    if p < 0
        error("a PreMPO has no power $p: it has no inverse")
    elseif p == 0
        id = lay_compacted!(PreMPO{R}(a.system), one(ComplexF64), Dict{Int, Matrix}(), nothing, 1)
        id.norm[] = 1.0
        return id
    elseif p == 1
        return a
    end
    h = a ^ (p ÷ 2)
    return isodd(p) ? h * h * a : h * h
end

# a literal negative exponent would go through `inv`, which a PreMPO has no method of
Base.literal_pow(::typeof(^), a::PreMPO, ::Val{p}) where p = a ^ p

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
`PreMPO`s of its factors multiplied, see `compacted_op`
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
    vanishes(pre)

whether the operator of `pre` is zero but for rounding: no channel and pieces of rounding size
next to the bound of its norm
"""
vanishes(pre::PreMPO) = all(==(1), pre.linkdims) &&
    all(normalized_norm(m) ≤ rounding_tol * pre.norm[] for ts in pre.terms for (_, _, _, _, m) in ts)

"""
    gathered_summands(g)

the operators the placed operator `g`, gathered, see `gathered`, is the sum of, its coefficient
left out
"""
function gathered_summands(g::IndexedOp)
    s = scalararg(g)
    return s isa SumOp ? s.subs : [s]
end

"""
    CompactedOp

an operator placed on sites compacted on a system, see `compacted_op`: its constant, its terms of
one site and its block of several, see `compacted`, the bounds of its norm and of the width of its
spectrum, see `term_size` and `term_width`, its flux, `nothing` when it has no term, and its terms
when it was made of them alone
"""
struct CompactedOp
    constant::Number
    singles::Dict
    block::Union{Nothing, ComBlock}
    norm::Float64
    width::Float64
    flux::Union{Nothing, QN}
    terms::Union{Nothing, Vector{Term}}
end

"""
    laid(system, c)

the `PreMPO`, of a single time function, of the operator compacted `c` on `system`
"""
function laid(system::System, c::CompactedOp)
    pre = lay_compacted!(PreMPO{Pure}(system), c.constant, c.singles, c.block, 1)
    pre.norm[] = c.norm
    pre.width[] = c.width
    return pre
end

"""
    compacted_op(system, op)

the operator `op`, placed on sites, compacted on `system`, see `CompactedOp`: from its terms,
see `lower`, but for its products that would develop into more terms than their factors have,
see `by_blocks`. Those are made of the `PreMPO`s of their factors, compacted the same way and
multiplied, equal factors in a row by squaring, see `^`, and added to that of the rest in one
reduction, see `combined`. Refused when it changes the charges by no definite amount, see
`FluxError`.
"""
function compacted_op(system::System, a::IndexedOp{IsOperator})
    check_indices(system, a)
    g = gathered(system, a)
    summands = gathered_summands(g)
    if !any(by_blocks, summands)
        terms = lower(system, a, g)
        flux = isempty(terms) ? nothing :
               is_charged(system) ? term_flux(Pure, system, a, first(terms)) : QN()
        return CompactedOp(compacted(Pure, system, terms)..., sum(term_size, terms; init = 0.0),
                           sum(t -> term_width(system, t), terms; init = 0.0), flux, terms)
    end
    parts, coefs, fluxes = PreMPO{Pure}[], Number[], QN[]
    rest = filter(x -> !by_blocks(x), summands)
    if !isempty(rest)
        c = compacted_op(system, SumOp(rest))
        if !isnothing(c.flux)
            push!(parts, laid(system, c))
            push!(coefs, scalarcoef(g))
            push!(fluxes, c.flux)
        end
    end
    for x in filter(by_blocks, summands)
        fs = prodsubs(scalararg(x))
        pre, flux = nothing, QN()
        k = 1
        while k ≤ length(fs)
            run = findfirst(f -> f ≠ fs[k], fs[k:end])
            n = isnothing(run) ? length(fs) - k + 1 : run - 1
            c = compacted_op(system, fs[k])
            if isnothing(c.flux)
                pre = nothing
                break
            end
            p = laid(system, c) ^ n
            pre = isnothing(pre) ? p : pre * p
            flux = sum(fill(c.flux, n); init = flux)
            k += n
        end
        # a product of definite fluxes that vanishes has none
        if !isnothing(pre) && !vanishes(pre)
            push!(parts, pre)
            push!(coefs, scalarcoef(g) * scalarcoef(x))
            push!(fluxes, flux)
        end
    end
    if length(unique(fluxes)) > 1
        flux_refusal(system, a, unique(fluxes))
    end
    if isempty(parts)
        return CompactedOp(zero(ComplexF64), Dict{Int, Matrix{ComplexF64}}(), nothing, 0.0, 0.0,
                           nothing, nothing)
    end
    singles, block = combined(parts, coefs)
    return CompactedOp(zero(ComplexF64), singles, block,
                       sum(abs(c) * p.norm[] for (p, c) in zip(parts, coefs)),
                       sum(isreal(c) ? abs(c) * p.width[] : 2abs(c) * p.norm[] for (p, c) in zip(parts, coefs)),
                       first(fluxes), nothing)
end

"""
    followed_operator(a)

the operator `K` a pure state follows as ``d\\psi/dt = K\\psi`` under the generator `a`,
refused when `a` holds a dissipator; an operator is itself, and a channel is refused
"""
function followed_operator(a::IndexedOp{IsGenerator})
    if !is_pure(a)
        error("$a holds a dissipator, which a pure state cannot follow: mix the state")
    end
    return pure_part(a)
end

followed_operator(a::IndexedOp{IsOperator}) = a
followed_operator(a::IndexedOp{IsChannel}) =
    error("$a is a channel, which acts on a density matrix: mix the state")

"""
    superoperator(a)

the channel or the generator `a`, an operator given for a density matrix being its evolver
"""
superoperator(a::IndexedOp{IsOperator}) = Evolver(a)
superoperator(a::IndexedOp) = a

PreMPO(state::State{Pure}, a::IndexedOp{IsOperator}) =
    laid(state.system, compacted_op(state.system, a))

function PreMPO(state::State{Pure}, as::Vector{<:IndexedOp{IsOperator}})
    pre = PreMPO{Pure}(state.system, length(as))
    for (ref, a) in enumerate(as)
        c = compacted_op(state.system, a)
        lay_compacted!(pre, c.constant, c.singles, c.block, ref)
        pre.norm[] += c.norm
        pre.width[] += c.width
    end
    return pre
end

PreMPO(state::State{Pure}, a::IndexedOp) = PreMPO(state, followed_operator(a))
PreMPO(state::State{Pure}, as::Vector{<:IndexedOp}) =
    PreMPO(state, IndexedOp{IsOperator}[ followed_operator(a) for a in as ])

PreMPO(state::State{Mixed}, a::IndexedOp) =
    lay_terms!(PreMPO{Mixed}(state.system), lower(state.system, superoperator(a)), 1)

function PreMPO(state::State{Mixed}, as::Vector{<:IndexedOp})
    pre = PreMPO{Mixed}(state.system, length(as))
    for (ref, a) in enumerate(as)
        lay_terms!(pre, lower(state.system, superoperator(a)), ref)
    end
    return pre
end

"""
    mpo_eltype(::PreMPO)

the element type of the tensors of the MPO built from `pre`, `Float64` at least, and real for
real matrices, the coefficients being real
"""
mpo_eltype(pre::PreMPO) =
    promote_type(Float64, mapreduce(t -> eltype(t[5]), promote_type, Iterators.flatten(pre.terms);
                                    init = Bool))

"""
    charge_mismatch(system, a, b, shown)

raise the error for an MPO on `system` whose terms do not all carry the same charge, those of
an operator or of the operators of a vector, `a` and `b` being the charges of two of them,
printed when `shown`, or of two partial terms entering the same channel; when they differ only
by quantities conserved strongly, that of `strong_refusal`
"""
function charge_mismatch(sys::System, a::QN, b::QN, shown::Bool)
    st = strong_names(sys)
    if only_strong(st, (a, b))
        strong_refusal("the MPO", st)
    end
    error("the terms of the MPO do not all carry the same charge, " *
          (shown ? "$(-a) and $(-b), " : "") *
          "so it has no definite flux and cannot be put on a system that conserves it")
end

"""
    channel_counts(pre)

the number of channels the terms of `pre` take on each link, from the one on the left of the
first site to the one on the right of the last, the term not yet begun included
"""
channel_counts(pre::PreMPO) = [ 1; pre.linkdims; 1 ]

"""
    mpo_charges(pre, coefs)

the charge of every channel of every link of the MPO, `q[i + 1][k]` for channel `k` of the
link on the right of site `i`: minus the sum of the fluxes of the factors on its left. The
first channel, the term not yet begun, has no charge; the last, the term finished, has minus
the flux of the operator, which every term must share, as must the pieces entering a channel.
Terms whose time function is zero are left out.
"""
function mpo_charges(pre::PreMPO, coefs::Vector{Float64})
    n = length(pre.system)
    tm = pre.terms
    q = [ Union{Nothing, QN}[ QN(); fill(nothing, d) ] for d in channel_counts(pre) ]
    total = nothing
    for i in 1:n
        for (l, r, u, ref, _) in tm[i]
            if coefs[ref] == 0
                continue
            end
            left = q[i][l]
            if isnothing(left)
                error("bug: channel $l of link $(i-1) has no charge")
            end
            # the link loses what the operator carries, which is the sign that makes the
            # tensor of the site a flux of zero once its two ends are oriented
            c = left - flux(u)
            if r == 1
                if !isnothing(total) && total ≠ c
                    charge_mismatch(pre.system, total, c, true)
                end
                total = c
            else
                if !isnothing(q[i+1][r]) && q[i+1][r] ≠ c
                    charge_mismatch(pre.system, q[i+1][r], c, false)
                end
                q[i+1][r] = c
            end
        end
    end
    total = something(total, QN())
    for qi in q
        qi[end] = total
        # a channel no term goes through, kept so that the numbering does not move
        replace!(qi, nothing => total)
    end
    return q
end

"""
    mpo_links(pre, coefs, charges, extra)

the links of the MPO of `pre`, each with `extra` channels besides those of `channel_counts`,
one for the term finished in `make_mpo`. On a charged system the channels carry the charges
`charges(pre, coefs)` gives.
"""
function mpo_links(pre::PreMPO, coefs::Vector{Float64}, charges, extra::Int)
    ds = channel_counts(pre) .+ extra
    if !is_charged(pre.system)
        return [ Index(d, "Link,l=$(i-1)") for (i, d) in enumerate(ds) ]
    end
    q = charges(pre, coefs)
    return [ Index([ q[i][k] => 1 for k in 1:d ]...; tags = "Link,l=$(i-1)")
             for (i, d) in enumerate(ds) ]
end

"""
    close_ends(ts, links, k)

the MPO of the tensors `ts` of the sites, its outer links `links[1]` and `links[end]` fixed on
channel 1, the term not yet begun, and on channel `k`, the term finished.
"""
function close_ends(ts::Vector{ITensor}, links::Vector{<:Index}, k::Int)
    ts[1] *= onehot(links[1] => 1)
    ts[end] *= onehot(dag(links[end]) => k)
    return MPO(ts)
end

"""
    site_array(elt, idx, llink, rlink)

the array, zero, of the tensor of the site of index `idx` in an MPO of links `llink` and
`rlink`, its axes those of `idx'`, `dag(idx)`, `llink` and `rlink`, to fill with `add_block!`
"""
site_array(elt::Type, idx::Index, llink::Index, rlink::Index) =
    zeros(elt, dim(idx), dim(idx), dim(llink), dim(rlink))

"""
    add_block!(a, l, r, m[, c])

add `c` times the matrix `m` of a one site operator, rows on `idx'`, to the block of channels
`l` and `r` of the array `a` of `site_array`, and return `a`
"""
function add_block!(a::Array{<:Number, 4}, l::Int, r::Int, m::AbstractMatrix, c::Number = 1)
    @views a[:, :, l, r] .+= c .* m
    return a
end

"""
    site_tensor(a, idx, llink, rlink)

the tensor of the array `a` of `site_array`; on a charged system only the blocks that are not
zero are kept, and they must all have the flux the charges of the links give
"""
site_tensor(a::Array{<:Number, 4}, idx::Index, llink::Index, rlink::Index) =
    ITensor(a, idx', dag(idx), dag(llink), rlink)

"""
    mpo_window(pre)

the sites from the first to the last on which `pre` has pieces, outside which its MPO is the
identity, `nothing` when it has none
"""
function mpo_window(pre::PreMPO)
    ks = findall(!isempty, pre.terms)
    return isempty(ks) ? nothing : first(ks):last(ks)
end

"""
    window_mpo(pre, coefs, sites)

the MPO of `pre`, `coefs` the values of its time functions, on the consecutive `sites` alone,
its outer links fixed on the term not yet begun and on the term finished, see `close_ends`: that
of the whole system restricted to `sites` when it is the identity outside, see `mpo_window`
"""
function window_mpo(pre::PreMPO{R}, coefs::Vector{Float64}, sites::UnitRange{Int}) where R
    sys = pre.system
    elt = mpo_eltype(pre)
    dims = channel_counts(pre)
    links = mpo_links(pre, coefs, mpo_charges, 1)
    ts = map(sites) do i
        idx = siteind(R, sys, i)
        llink, rlink = links[i], links[i+1]
        a = site_array(elt, idx, llink, rlink)
        id = Matrix{elt}(I, dim(idx), dim(idx))
        add_block!(a, 1, 1, id)
        add_block!(a, 1 + dims[i], 1 + dims[i+1], id)
        for (l, r, _, ref, m) in pre.terms[i]
            c = coefs[ref]
            if c ≠ 0
                # the coefficient of a term goes on its closing piece alone, the only one
                # with r == 1
                if r == 1
                    add_block!(a, l, 1 + dims[i+1], m, c)
                else
                    add_block!(a, l, r, m)
                end
            end
        end
        return site_tensor(a, idx, llink, rlink)
    end
    return close_ends(ts, links[first(sites):last(sites)+1], 2)
end

"""
    make_mpo(::PreMPO[, coefs])
    make_mpo(::State, op)

the MPO of the operator `op` on the system of the state, given with the state or as
`PreMPO(state, op)`. For a vector of operators, built from its `PreMPO`, `coefs` holds the
value of each time function, one real number per operator (default `[1.]`, a single operator).

# Examples

    mpo = make_mpo(state, sum(Z(i) * Z(i + 1) for i in 1:9))
    mpo = make_mpo(ρ, -im * sum(Z(i) * Z(i + 1) for i in 1:9) + sum(Dissipator(Sm)(i) for i in 1:10))
"""
make_mpo(pre::PreMPO, coefs::AbstractVector = [1.]) =
    window_mpo(pre, check_coefs(pre, coefs), 1:length(pre.system))

make_mpo(state::State, a) = make_mpo(PreMPO(state, a))

"""
    cached_mpos(build, pre, method, key)

the MPOs `build()` gives for the solver `method` and the parameters `key`, kept in `pre` for the
next call: only the last ones of each method, other parameters replacing them
"""
function cached_mpos(build, pre::PreMPO, method::Symbol, key)
    lock(pre.lock) do
        k, v = get(pre.cache, method, (nothing, nothing))
        if isnothing(v) || !isequal(k, key)
            v = build()
            pre.cache[method] = (key, v)
        end
        return v
    end
end
