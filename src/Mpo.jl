# The MPO of an operator on the system of a state: PreMPO, which lays its terms, compacted, site
# by site, and make_mpo, which builds the MPO from it.

export PreMPO, make_mpo

"""
    PreMPO(state, op)

the operator `op`, placed on sites, laid on the system of `state` as the pieces of its MPO,
its terms compacted into as few channels as a triangular MPO can have; `make_mpo` builds the
MPO from it. A vector of operators gives each its own time function, see `make_mpo`.
"""
struct PreMPO{R <: PM}
    system::System
    linkdims::Vector{Int}
    # the pieces on each site, see `add_term!`
    terms::Vector{Vector{Tuple{Int, Int, ITensor, Int, Matrix}}}
    # the number of time functions the terms take, one per operator of a vector: it cannot be
    # counted from `terms`, a whole operator vanishing on its sites
    nterms::Int
    function PreMPO{R}(system::System, nterms::Int = 1) where R
        n = length(system)
        return new{R}(system, fill(1, n - 1),
                      [ Tuple{Int, Int, ITensor, Int, Matrix}[] for _ in 1:n ], nterms)
    end
end

"""
    add_term!(pre, k, l, r, m, ref)

add to site `k` of `pre`, and return it, the piece of matrix `m` going from channel `l` to
channel `r` and taking the time function `ref`, with its tensor on the index of the site
"""
function add_term!(pre::PreMPO{R}, k::Int, l::Int, r::Int, m::Matrix, ref::Int) where R
    idx = siteind(R, pre.system, k)
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
time function `ref`: the constant on the first site, the terms of one site on their site, and
the channels of the terms of several sites numbered after those their links already have
"""
function lay_terms!(pre::PreMPO, terms::Vector{Term}, ref::Int)
    sys = pre.system
    constant, singles, b = compacted(sys, terms)
    if !iszero(constant)
        add_term!(pre, 1, 1, 1, Matrix(constant * I, dim(sys[1]), dim(sys[1])), ref)
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
        m = combination_matrix(b, j, comb, dim(sys[k]))
        if !iszero(m)
            add_term!(pre, k, l == 0 ? 1 : offs[j-1] + l, r == 0 ? 1 : offs[j] + r, m, ref)
        end
    end
    return pre
end

PreMPO(state::State{Pure}, a::IndexedOp{IsOperator}) =
    lay_terms!(PreMPO{Pure}(state.system), lower(state.system, a), 1)

function PreMPO(state::State{Pure}, as::Vector{<:IndexedOp{IsOperator}})
    pre = PreMPO{Pure}(state.system, length(as))
    for (ref, a) in enumerate(as)
        lay_terms!(pre, lower(state.system, a), ref)
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

raise the error for an operator on `system` whose terms do not all carry the same charge, `a`
and `b` being the charges of two of them, printed when `shown`, or of two partial terms
entering the same channel
"""
function charge_mismatch(sys::System, a::QN, b::QN, shown::Bool)
    st = strong_names(sys)
    only_strong = !isempty(st) && weak_qn(a, st, String[]) == weak_qn(b, st, String[])
    error("the terms of this operator do not all carry the same charge, " *
          (shown ? "$(-a) and $(-b), " : "") *
          (only_strong ?
              "which only a strong symmetry tells apart: drop `strong` or weaken the state" :
              "so it has no definite flux and cannot be put on a system that conserves it"))
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
    make_mpo(::PreMPO[, coefs])
    make_mpo(::State, op)

the MPO of the operator `op` on the system of the state, given with the state or as
`PreMPO(state, op)`. For a vector of operators, built from its `PreMPO`, `coefs` holds the
value of each time function, one real number per operator (default `[1.]`, a single operator).

# Examples

    mpo = make_mpo(state, sum(Z(i) * Z(i + 1) for i in 1:9))
"""
function make_mpo(pre::PreMPO{R}, coefs::AbstractVector = [1.]) where R
    coefs = check_coefs(pre, coefs)
    sys = pre.system
    elt = mpo_eltype(pre)
    dims = channel_counts(pre)
    links = mpo_links(pre, coefs, mpo_charges, 1)
    ts = map(1:length(sys)) do i
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
    return close_ends(ts, links, 2)
end

make_mpo(state::State, a) = make_mpo(PreMPO(state, a))
