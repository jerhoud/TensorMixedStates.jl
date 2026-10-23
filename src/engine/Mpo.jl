# The MPO of an operator on the system of a state: PreMPO, which lays its terms, compacted, site
# by site, and make_mpo, which builds the MPO from it. On a mixed state, the MPO of a channel or
# a generator, its superoperator laid on the mixed indices.

export PreMPO, make_mpo

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
    Piece

a piece of an MPO laid on a site of a `PreMPO`, see `add_piece!`: its matrix `m`, in the order of
the index of the site, from the channel `from` of the link on its left to the channel `to` of
that on its right, `1` being the term not yet begun on the left and the term finished on the
right, see `window_mpo`, taking the time function `ref`, and the charge it carries, see
`block_flux`, which `mpo_charges` reads
"""
struct Piece
    from::Int
    to::Int
    ref::Int
    m::Matrix
    flux::QN
end

"""
    PreMPO(state, op; limits = Limits())
    PreMPO(state, [A, B => g]; limits = Limits())

the operator `op`, placed on sites, laid on the system of `state` as the pieces of its MPO,
its terms compacted into as few channels as a triangular MPO can have; `make_mpo` builds the
MPO from it. A vector of terms, each an operator or an operator with its real function of time,
is the evolver depending on time ``A + g(t) B`` of `evolve`, each term laid apart, see
`make_mpo`. On a pure state, a generator with no dissipator,
``\\rho \\mapsto K\\rho + \\rho K^\\dagger``, gives the MPO of `K`, which the state follows as
``d\\psi/dt = K\\psi``. On a mixed state, a channel or a generator gives the MPO of its
superoperator, and an operator `K` that of `Evolver(K)`.

`limits` truncates the MPO in its triangular form, see `Limits`: across each link, the singular
values of the part of the operator acting on both sides, for the Hilbert-Schmidt norm
normalized on each site, are cut as those of a bond of a state, the bond dimension of the MPO
being two more than its channels; by default nothing is cut but rounding.

`pre'` is the `PreMPO` of the adjoint, its time functions kept, and `a * b`, `a ^ p`, `a + b`,
`a - b`, `c * a` and `sum(v)`, of `PreMPO`s of one system each of a single term not depending on
time, those of the product, of the power and of the linear combinations of their operators:
built with as few channels as their rank across each link, at a cost linear in the number of
sites, where the terms of a product grow as the product of theirs. `PreMPO(state, op)` builds so
the products of sums of `op` that would develop into more terms than their factors have, and, on
a mixed state, its gates and dissipators of operators of many terms from the sides ``K\\rho``
and ``\\rho K^\\dagger`` of these.
"""
struct PreMPO{R <: PM}
    system::System
    # the operator of each term, compacted, which its pieces lay, see `lay_compacted!`
    parts::Vector{CompactedOp{R}}
    # the function of time of each term, `nothing` for an operator that does not depend on time
    functions::Union{Nothing, Vector{Any}}
    linkdims::Vector{Int}
    # the pieces on each site, see `add_piece!`
    pieces::Vector{Vector{Piece}}
    orders::Vector{Vector{Int}}
    # the last MPOs built for each method of the solvers, with what they were built for, see
    # `cached_mpos`
    cache::Dict{Symbol, Tuple{Any, Any}}
    lock::ReentrantLock
    # for a pure state, a bound of the width of the spectrum of the operator, see `term_width`
    width::Float64
    # a bound of the norm of the operator, see `term_size`
    norm::Float64
end

# its system, and its number of terms when they depend on time, its pieces being internal
show(io::IO, pre::PreMPO{R}) where R =
    print(io, "PreMPO{", R, "}(", pre.system,
          isnothing(pre.functions) ? "" : ", $(length(pre.parts)) terms", ")")

"""
    on_index(pre, k, m)

the matrix `m` of an operator on site `k`, its rows and columns in the order of the index of
the site, see `block_orders`
"""
on_index(::PreMPO{Pure}, ::Int, m::Matrix) = m
on_index(pre::PreMPO{Mixed}, k::Int, m::Matrix) = m[pre.orders[k], pre.orders[k]]

"""
    add_piece!(pre, k, l, r, m, ref)

add to site `k` of `pre`, and return it, the piece of matrix `m`, in the order of the site, going
from channel `l` to channel `r` and taking the time function `ref`, see `Piece`
"""
function add_piece!(pre::PreMPO{R}, k::Int, l::Int, r::Int, m::Matrix, ref::Int) where R
    sys = pre.system
    flux = is_charged(sys) ? block_flux(R, m, "the MPO", sys[k]) : QN()
    push!(pre.pieces[k], Piece(l, r, ref, on_index(pre, k, m), flux))
    return pre
end

"""
    check_coefs(pre, coefs)

`coefs` as a `Vector{Float64}`, refused unless it holds one real value per term of `pre`, that
of its time function
"""
function check_coefs(pre::PreMPO, coefs::AbstractVector)
    if length(coefs) ≠ length(pre.parts)
        error("a PreMPO of $(length(pre.parts)) terms takes as many coefficients, got " *
              "$(length(coefs))")
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
    constant_coefs(pre)

the coefficients of the terms of `pre`, all one, refused when `pre` depends on time
"""
function constant_coefs(pre::PreMPO)
    if !isnothing(pre.functions)
        error("the PreMPO depends on time: give the values of its time functions, " *
              "make_mpo(pre, coefs)")
    end
    return ones(length(pre.parts))
end

"""
    lay_compacted!(pre, c, ref)

add to `pre`, and return it, the operator compacted `c`, see `CompactedOp`, taking the time
function `ref`: the constant on the first site of the others, so that the MPO is the identity
outside their sites, see `mpo_window`, the terms of one site on their site, and the channels of
the block numbered after those their links already have
"""
function lay_compacted!(pre::PreMPO{R}, c::CompactedOp{R}, ref::Int) where R
    sys = pre.system
    b = c.block
    if !iszero(c.constant)
        ks = [ collect(keys(c.singles)) ; isnothing(b) ? Int[] : [b.start] ]
        k = isempty(ks) ? 1 : minimum(ks)
        d = block_dim(R, sys[k])
        add_piece!(pre, k, 1, 1, Matrix(c.constant * I, d, d), ref)
    end
    for k in sort!(collect(keys(c.singles)))
        add_piece!(pre, k, 1, 1, c.singles[k], ref)
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
            add_piece!(pre, k, l == 0 ? 1 : offs[j-1] + l, r == 0 ? 1 : offs[j] + r, m, ref)
        end
    end
    return pre
end

# the PreMPO on `system` of the operators compacted `parts`, one per term, laid in their order,
# see `lay_compacted!`, and of their time functions `functions`; a comment, not a docstring,
# which would put a method of no user in the reference of the exported PreMPO
function PreMPO(system::System, parts::Vector{CompactedOp{R}},
                functions::Union{Nothing, Vector{Any}} = nothing) where R
    n = length(system)
    pre = PreMPO{R}(system, parts, functions, fill(1, n - 1), [ Piece[] for _ in 1:n ],
                    block_orders(R, system), Dict{Symbol, Tuple{Any, Any}}(), ReentrantLock(),
                    sum(c -> c.width, parts; init = 0.0), sum(c -> c.norm, parts; init = 0.0))
    for (ref, c) in enumerate(parts)
        lay_compacted!(pre, c, ref)
    end
    return pre
end

"""
    single_part(pres, what)

the operator compacted of each of the `PreMPO`s `pres`, refused when they are on two systems or
one depends on time or has several terms, `what` naming the operation in the refusals
"""
function single_part(pres, what::String)
    sys = first(pres).system
    if any(p -> p.system ≠ sys, pres)
        error("the PreMPOs of a $what are on two systems")
    end
    done = what == "product" ? "multiplied" : "added or scaled"
    for p in pres
        if !isnothing(p.functions)
            error("a PreMPO depending on time is not $done")
        elseif length(p.parts) ≠ 1
            error("a PreMPO of several terms is not $done")
        end
    end
    return [ only(p.parts) for p in pres ]
end

"""
    adjoint(pre)

the `PreMPO` of the adjoint of the operator of `pre`, each term apart, its time functions, of
real values, kept
"""
Base.adjoint(pre::PreMPO) = PreMPO(pre.system, [ c' for c in pre.parts ], pre.functions)

"""
    a * b

the `PreMPO` of the product of the operators of `a` and `b`, of one system and representation
and of a single time function each, see `product`
"""
(a::PreMPO{R} * b::PreMPO{R}) where {R <: PM} = product(a, b, Limits())

(::PreMPO * ::PreMPO) = error("the PreMPOs of a product are of two representations")

"""
    product(a, b, limits)

`a * b`, truncated by `limits`, see `truncated_block`
"""
function product(a::PreMPO{R}, b::PreMPO{R}, limits::Limits) where R
    x, y = single_part([a, b], "product")
    return PreMPO(a.system, [product(x, y, limits)])
end

"""
    linear_combination(pres, coefs)

the `PreMPO` of the sum of the operators of `pres` times `coefs`, see `combined`
"""
function linear_combination(pres::Vector{<:PreMPO{R}}, coefs::Vector{<:Number}) where R
    sys = first(pres).system
    return PreMPO(sys, [combined(sys, collect(zip(single_part(pres, "sum"), coefs)))])
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
products, and `a^0` the identity, see `power`
"""
(a::PreMPO{R} ^ p::Integer) where {R <: PM} =
    PreMPO(a.system, [power(only(single_part([a], "product")), p)])

# a literal negative exponent would go through `inv`, which a PreMPO has no method of
Base.literal_pow(::typeof(^), a::PreMPO, ::Val{p}) where p = a ^ p

"""
    followed_operator(a)

the operator `K` a pure state follows as ``d\\psi/dt = K\\psi`` under the generator `a`,
refused when `a` holds a dissipator; an operator is itself, and a channel is refused
"""
function followed_operator(a::IndexedOp{IsGenerator})
    if !is_pure(a)
        error("$(obs_name(a)) holds a dissipator, which a pure state cannot follow: mix the state")
    end
    return pure_part(a)
end

followed_operator(a::IndexedOp{IsOperator}) = a
followed_operator(a::IndexedOp{IsChannel}) =
    error("$(obs_name(a)) is a channel, which acts on a density matrix: mix the state")

"""
    superoperator(a)

the channel or the generator `a`, an operator given for a density matrix being its evolver
"""
superoperator(a::IndexedOp{IsOperator}) = Evolver(a)
superoperator(a::IndexedOp) = a

"""
    split_terms(terms)

the operators of the vector of terms `terms`, see `PreMPO`, and their time functions, `nothing`
when no term has one, a constant term among others taking `Returns(1.)`; the vectors paired of
1.x, `[A, B] => [f, g]`, are refused
"""
function split_terms(terms::AbstractVector)
    ops = IndexedOp[]
    fs = Any[]
    for t in terms
        if t isa IndexedOp
            push!(ops, t)
            push!(fs, nothing)
        elseif t isa Pair && first(t) isa IndexedOp && !(last(t) isa AbstractVector)
            push!(ops, first(t))
            push!(fs, last(t))
        elseif t isa Pair && (first(t) isa Union{AbstractVector, PreMPO} ||
                              last(t) isa AbstractVector)
            error("an evolver depending on time gives each term its function, [A, B => g] for " *
                  "A + g(t) B, where 1.x paired vectors, [A, B] => [f, g]")
        else
            error("$(repr(t)) is not a term of an evolver: an operator A, or A => f, f its " *
                  "function of time")
        end
    end
    return ops, all(isnothing, fs) ? nothing : Any[ isnothing(f) ? Returns(1.) : f for f in fs ]
end

# a single operator, or a single term, being a vector of one
PreMPO(state::State, a::Union{IndexedOp, Pair}; limits::Limits = Limits()) =
    PreMPO(state, [a]; limits)

function PreMPO(state::State{Pure}, terms::AbstractVector; limits::Limits = Limits())
    sys = state.system
    ops, fs = split_terms(terms)
    parts = CompactedOp{Pure}[ compacted_op(sys, followed_operator(a), limits) for a in ops ]
    return PreMPO(sys, parts, fs)
end

function PreMPO(state::State{Mixed}, terms::AbstractVector; limits::Limits = Limits())
    sys = state.system
    ops, fs = split_terms(terms)
    parts = CompactedOp{Mixed}[ compacted_super(sys, superoperator(a), limits) for a in ops ]
    return PreMPO(sys, parts, fs)
end

"""
    truncated(pre, limits)

`pre` truncated by `limits`, see `truncated_block`, each of its terms apart, its time functions
and the bounds of its norm and of its width kept, and `pre` itself for `Limits()`
"""
function truncated(pre::PreMPO, limits::Limits)
    if limits == Limits()
        return pre
    end
    return PreMPO(pre.system, [ combined(pre.system, [(c, 1)], limits) for c in pre.parts ],
                  pre.functions)
end

"""
    mpo_eltype(::PreMPO)

the element type of the tensors of the MPO built from `pre`, `Float64` at least, and real for
real matrices, the coefficients being real
"""
mpo_eltype(pre::PreMPO) =
    promote_type(Float64, mapreduce(p -> eltype(p.m), promote_type, Iterators.flatten(pre.pieces);
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
    q = [ Union{Nothing, QN}[ QN(); fill(nothing, d) ] for d in channel_counts(pre) ]
    total = nothing
    for i in 1:n
        for p in pre.pieces[i]
            l, r = p.from, p.to
            if coefs[p.ref] == 0
                continue
            end
            left = q[i][l]
            if isnothing(left)
                error("bug: channel $l of link $(i-1) has no charge")
            end
            # the link loses what the operator carries, which is the sign that makes the
            # tensor of the site a flux of zero once its two ends are oriented
            c = left - p.flux
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
    ks = findall(!isempty, pre.pieces)
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
        for p in pre.pieces[i]
            c = coefs[p.ref]
            if c ≠ 0
                # the coefficient of a term goes on its closing piece alone, the only one
                # with p.to == 1
                if p.to == 1
                    add_block!(a, p.from, 1 + dims[i+1], p.m, c)
                else
                    add_block!(a, p.from, p.to, p.m)
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
`PreMPO(state, op)`. For a vector of terms, built from its `PreMPO`, `coefs` holds the value of
the time function of each, one real number per term, one for all by default, which a `PreMPO`
depending on time refuses.

# Examples

    mpo = make_mpo(state, sum(Z(i) * Z(i + 1) for i in 1:9))
    mpo = make_mpo(ρ, -im * sum(Z(i) * Z(i + 1) for i in 1:9) +
                      sum(Dissipator(Sm)(i) for i in 1:10))
"""
make_mpo(pre::PreMPO, coefs::AbstractVector = constant_coefs(pre)) =
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
