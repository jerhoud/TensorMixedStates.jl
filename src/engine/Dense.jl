# The terms of an operator, a channel or a generator acting on a dense state: the vector of a pure
# state, or the density matrix of a mixed one, on the basis of the product states ordered as
# `kron` orders them, the first site varying the slowest.

"""
the dimensions of the sites of `system`, in their order
"""
site_dims(system::System) = [ dim(s) for s in system.sites ]

"""
the numbers of elements after and before the axis of site `k` in an array of the sites of
dimensions `dims` in the order of `kron`, the first site the slowest
"""
around(dims::Vector{Int}, k::Int) = (prod(dims[k+1:end]; init = 1), prod(dims[1:k-1]; init = 1))

"""
the array `x` in the order of `kron` on one axis per site of dimensions `dims`, in their order
"""
on_sites(x::AbstractArray, dims::Vector{Int}) =
    permutedims(reshape(x, reverse(dims)...), length(dims):-1:1)

"""
the array `a` of one axis per site, in their order, as the vector in the order of `kron`
"""
in_kron(a::AbstractArray) = vec(permutedims(a, ndims(a):-1:1))

"""
the signs ``(-1)^{\\sum p_k p_l}`` of the product states of `sites`, in the order of `kron`, the
sum over the pairs `(k, l)` of `pairs`, `p` the fermionic parity of the value of a site;
`nothing` when no pair joins two fermionic sites
"""
function pair_signs(sites::Vector{<:AbstractSite}, pairs)
    n = length(sites)
    parities = [ collect(fermion_parities(s)) for s in sites ]
    fermionic = [ (k, l) for (k, l) in pairs if any(isodd, parities[k]) && any(isodd, parities[l]) ]
    if isempty(fermionic)
        return nothing
    end
    on_axis(v, k) = reshape(v, ntuple(j -> j == k ? length(v) : 1, n))
    exponent = zeros(Int, [ dim(s) for s in sites ]...)
    for (k, l) in fermionic
        exponent = exponent .+ on_axis(parities[k], k) .* on_axis(parities[l], l)
    end
    return in_kron(@. ifelse(iseven(exponent), 1.0, -1.0))
end

"""
the array `y` plus `coef` times the array `x`, both seen as `(a, d, c)`, with the matrix `m`
applied to the middle axis of `x`
"""
function axis_added!(y::AbstractArray, coef::Number, m::AbstractMatrix, x::AbstractArray, a::Int,
                     d::Int, c::Int)
    t = reshape(x, a, d, c)
    u = reshape(y, a, d, c)
    # element by element of `m`, which holds many zeros, without permuting the array
    for i in 1:d, j in 1:d
        if !iszero(m[j, i])
            f = coef * m[j, i]
            @views u[:, j, :] .+= f .* t[:, i, :]
        end
    end
    return y
end

"""
the array `x`, seen as `(a, d, c)`, with the matrix `m` applied to its middle axis, as a vector
"""
axis_acted(m::AbstractMatrix, x::AbstractArray, a::Int, d::Int, c::Int) =
    axis_added!(zeros(promote_type(eltype(m), eltype(x)), a * d * c), true, m, x, a, d, c)

"""
the terms `ts`, their strings laid, prepared once as a function of a vector of the sites of
dimensions `dims` followed by `extra` columns, the slowest, applying them to it; with `diagonal`,
their diagonal terms summed into one diagonal, kept as large as the vector
"""
function prepared_on_axes(::Type{R}, system::System, ts::Vector{Term}, dims::Vector{Int};
                          diagonal::Bool = true) where R
    laid = [ laid_strings(R, system, t) for t in ts ]
    blocks = [ [ b for b in t.blocks if !isone(b.m) ] for t in laid ]
    summed = [ diagonal && all(b -> isdiag(b.m), bs) for bs in blocks ]
    T = reduce(promote_type, [ [ typeof(t.coef) for t in laid ];
                               [ eltype(b.m) for bs in blocks for b in bs ] ]; init = Bool)
    D = prod(dims)
    sizes = [ around(dims, k) for k in eachindex(dims) ]
    # the diagonal terms summed into one diagonal, applied in one pass
    d = nothing
    if any(summed)
        d = zeros(T, D)
        u = Vector{T}(undef, D)
        for (t, bs, s) in zip(laid, blocks, summed)
            if s
                fill!(u, t.coef)
                for b in bs
                    after, before = sizes[b.site]
                    reshape(u, after, dims[b.site], before) .*= reshape(diag(b.m), 1, :, 1)
                end
                d .+= u
            end
        end
    end
    rest = [ (t.coef, bs) for (t, bs, s) in zip(laid, blocks, summed) if !s ]
    deep = any(r -> length(r[2]) > 1, rest)
    return function (x::AbstractVector, extra::Int)
        S = promote_type(T, eltype(x))
        y = Vector{S}(undef, length(x))
        if isnothing(d)
            fill!(y, zero(S))
        else
            reshape(y, D, extra) .= d .* reshape(x, D, extra)
        end
        # the blocks of each other term but its last through two buffers, the last added to `y`
        buffers = deep ? (similar(y), similar(y)) : nothing
        for (coef, bs) in rest
            if isempty(bs)
                y .+= coef .* x
                continue
            end
            z = x
            for (p, b) in enumerate(bs)
                k = b.site
                after, before = sizes[k]
                if p == length(bs)
                    axis_added!(y, coef, b.m, z, after, dims[k], before * extra)
                else
                    w = buffers[isodd(p) ? 1 : 2]
                    fill!(w, zero(S))
                    z = axis_added!(w, true, b.m, z, after, dims[k], before * extra)
                end
            end
        end
        return y
    end
end

"""
    prepared_terms(R, system, terms; diagonal = true)

the terms `terms` of an operator, see `lower`, prepared once as a function applying them to a
vector or to the kets of a matrix, for `Pure`, or those of a superoperator as a function applying
them to a density matrix, for `Mixed`; `diagonal = false` keeps no array as large as the state,
for a function prepared among many
"""
function prepared_terms(::Type{<:Ket}, system::System, ts::Vector{Term}; diagonal::Bool = true)
    dims = site_dims(system)
    D = prod(dims)
    act = prepared_on_axes(Pure, system, ts, dims; diagonal)
    return function (x::AbstractVecOrMat)
        if size(x, 1) ≠ D
            error("a dense state of the system has $D rows, not $(size(x, 1))")
        end
        return reshape(act(vec(x), size(x, 2)), size(x))
    end
end

function prepared_terms(::Type{Mixed}, system::System, ts::Vector{Term}; diagonal::Bool = true)
    dims = site_dims(system)
    D = prod(dims)
    # the elements of each site together, as the blocks of a superoperator take them
    kets, bras = pair_positions(dims)
    act = prepared_on_axes(Mixed, system, ts, dims .^ 2; diagonal)
    return function (x::AbstractMatrix)
        if size(x) ≠ (D, D)
            error("a density matrix of the system is $D × $D, not $(size(x, 1)) × $(size(x, 2))")
        end
        w = similar(vec(x))
        for b in 1:D, a in 1:D
            w[1 + kets[a] + bras[b]] = x[a, b]
        end
        y = act(w, 1)
        return [ y[1 + k + b] for k in kets, b in bras ]
    end
end

"""
    dense_acted(R, system, terms, x)

the terms `terms` applied once to `x`, see `prepared_terms`
"""
dense_acted(::Type{R}, system::System, ts::Vector{Term}, x::AbstractVecOrMat) where R =
    prepared_terms(R, system, ts)(x)

"""
the trace of the operator of terms `ts`, see `lower`, times the density matrix `ρ` of `system`:
for each term, its Jordan-Wigner strings laid, each product state `a` and each `b` the nonzero
elements of its blocks reach from `a`, the product of these elements times `ρ[b, a]`
"""
function dense_trace(system::System, ts::Vector{Term}, ρ::AbstractMatrix)
    dims = site_dims(system)
    strides = [ first(around(dims, k)) for k in eachindex(dims) ]
    total = zero(ComplexF64)
    for t in ts
        laid = laid_strings(Pure, system, t)
        blocks = [ b for b in laid.blocks if !isone(b.m) ]
        # the nonzero elements of each row of each block, as (column, value)
        rows = [ [ [ (j, b.m[i, j]) for j in axes(b.m, 2) if !iszero(b.m[i, j]) ]
                   for i in axes(b.m, 1) ] for b in blocks ]
        function walk(a, p, b, v)
            if p > length(blocks)
                return v * ρ[b + 1, a + 1]
            end
            k = blocks[p].site
            ak = (a ÷ strides[k]) % dims[k]
            return sum((walk(a, p + 1, b + (j - 1 - ak) * strides[k], v * x)
                        for (j, x) in rows[p][ak + 1]); init = zero(ComplexF64))
        end
        total += laid.coef * sum(a -> walk(a, 1, a, one(ComplexF64)), 0:size(ρ, 1) - 1)
    end
    return total
end

"""
    summed_action(R, system, parts, terms_of, prepared; diagonal = true)

the sum of `parts` as a function on the arrays of representation `R`: the terms `terms_of(x)` of
the parts that have some gathered into one list, see `prepared_terms`, the others by
`prepared(x)`
"""
function summed_action(::Type{R}, system::System, parts, terms_of, prepared;
                       diagonal::Bool = true) where R
    terms = Term[]
    rest = Function[]
    for x in parts
        u = terms_of(x)
        if isnothing(u)
            push!(rest, prepared(x))
        else
            append!(terms, u)
        end
    end
    merged = reduced(terms)
    on = isempty(merged) ? zero : prepared_terms(R, system, merged; diagonal)
    return y -> on(y) + sum(f -> f(y), rest; init = zero(y))
end

"""
    prepared_acted(system, g; diagonal = true)

the placed operator `g`, gathered, as a function applying it to a vector or to the kets of a
matrix, by its terms where it develops, otherwise factor by factor; its flux is not checked; see
`prepared_terms` for `diagonal`
"""
function prepared_acted(system::System, g::IndexedOp{IsOperator}; diagonal::Bool = true)
    ts = developed(system, g)
    if !isnothing(ts)
        return prepared_terms(Pure, system, ts; diagonal)
    end
    c, s = scalarcoef(g), scalararg(g)
    if s isa ProdOp
        # the factors prepared together, each without a diagonal as large as the state
        fs = [ prepared_acted(system, f; diagonal = false) for f in s.subs ]
        return x -> c * foldr((f, y) -> f(y), fs; init = x)
    elseif s isa SumOp
        f = summed_action(Pure, system, s.subs, t -> developed(system, t),
                          t -> prepared_acted(system, t; diagonal); diagonal)
        return x -> c * f(x)
    end
    error("$(obs_name(g)) neither develops nor is a sum or a product")
end

"""
    gathered_acted(system, g, x)

the placed operator `g`, gathered, applied to the vector `x`, or to the kets of the matrix `x`,
see `prepared_acted`
"""
gathered_acted(system::System, g::IndexedOp{IsOperator}, x::AbstractVecOrMat) =
    prepared_acted(system, g)(x)
