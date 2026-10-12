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
    function PreMPO{R}(system::System, nterms::Int = 1) where R
        n = length(system)
        return new{R}(system, fill(1, n - 1),
                      [ Tuple{Int, Int, ITensor, Int, Matrix}[] for _ in 1:n ], nterms,
                      block_orders(R, system), Dict{Symbol, Tuple{Any, Any}}(), ReentrantLock())
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
    lay_terms!(pre, terms, ref)

add to `pre`, and return it, the terms of an operator, compacted, see `compacted`, taking the
time function `ref`: the constant on the first site of the others, so that the MPO is the
identity outside their sites, see `mpo_window`, the terms of one site on their site, and the
channels of the terms of several sites numbered after those their links already have
"""
function lay_terms!(pre::PreMPO{R}, terms::Vector{Term}, ref::Int) where R
    sys = pre.system
    constant, singles, b = compacted(R, sys, terms)
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
    followed_operator(a)

the operator `K` a pure state follows as ``d\\psi/dt = K\\psi`` under the generator `a`,
refused when `a` holds a dissipator
"""
function followed_operator(a::IndexedOp{IsGenerator})
    if !is_pure(a)
        error("$a holds a dissipator, which a pure state cannot follow: mix the state")
    end
    return pure_part(a)
end

"""
    superoperator(a)

the channel or the generator `a`, an operator given for a density matrix being its evolver
"""
superoperator(a::IndexedOp{IsOperator}) = Evolver(a)
superoperator(a::IndexedOp) = a

PreMPO(state::State{Pure}, a::IndexedOp{IsOperator}) =
    lay_terms!(PreMPO{Pure}(state.system), lower(state.system, a), 1)

function PreMPO(state::State{Pure}, as::Vector{<:IndexedOp{IsOperator}})
    pre = PreMPO{Pure}(state.system, length(as))
    for (ref, a) in enumerate(as)
        lay_terms!(pre, lower(state.system, a), ref)
    end
    return pre
end

PreMPO(state::State{Pure}, a::IndexedOp{IsGenerator}) = PreMPO(state, followed_operator(a))
PreMPO(state::State{Pure}, as::Vector{<:IndexedOp{IsGenerator}}) =
    PreMPO(state, map(followed_operator, as))

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
