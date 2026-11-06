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
the array `x`, seen as `(a, d, c)`, with the matrix `m` applied to its middle axis, as a vector
"""
function axis_acted(m::AbstractMatrix, x::AbstractArray, a::Int, d::Int, c::Int)
    t = reshape(x, a, d, c)
    y = zeros(promote_type(eltype(m), eltype(x)), a, d, c)
    # element by element of `m`, which holds many zeros, without permuting the array
    for i in 1:d, j in 1:d
        if !iszero(m[j, i])
            @views y[:, j, :] .+= m[j, i] .* t[:, i, :]
        end
    end
    return vec(y)
end

"""
the sum of the terms `ts` applied to `x`, of `dims` on the sites and `extra` columns, the slowest
"""
function terms_acted(::Type{R}, system::System, ts::Vector{Term}, x::AbstractVector,
                     dims::Vector{Int}, extra::Int) where R
    return sum(ts; init = zero(x)) do t
        laid = laid_strings(R, system, t)
        foldl(laid.blocks; init = x) do y, b
            k = b.site
            after, before = around(dims, k)
            isone(b.m) ? y : axis_acted(b.m, y, after, dims[k], before * extra)
        end * laid.coef
    end
end

"""
    dense_acted(R, system, terms, x)

the terms `terms` of an operator applied to the vector `x` or to the kets of the matrix `x`, for
`Pure`, or those of a superoperator applied to the density matrix `x`, for `Mixed`, see `lower`
"""
function dense_acted(::Type{Pure}, system::System, ts::Vector{Term}, x::AbstractVecOrMat)
    dims = site_dims(system)
    if size(x, 1) ≠ prod(dims)
        error("a dense state of the system has $(prod(dims)) rows, not $(size(x, 1))")
    end
    return reshape(terms_acted(Pure, system, ts, vec(x), dims, size(x, 2)), size(x))
end

function dense_acted(::Type{Mixed}, system::System, ts::Vector{Term}, x::AbstractMatrix)
    dims = site_dims(system)
    D = prod(dims)
    if size(x) ≠ (D, D)
        error("a density matrix of the system is $D × $D, not $(size(x, 1)) × $(size(x, 2))")
    end
    # the elements of each site together, as the blocks of a superoperator take them
    kets, bras = pair_positions(dims)
    w = similar(vec(x))
    for b in 1:D, a in 1:D
        w[1 + kets[a] + bras[b]] = x[a, b]
    end
    y = terms_acted(Mixed, system, ts, w, dims .^ 2, 1)
    return [ y[1 + k + b] for k in kets, b in bras ]
end

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
    gathered_acted(system, g, x)

the placed operator `g`, gathered, see `gathered`, applied to the vector `x`, or to the kets of
the matrix `x`: through its terms when it develops, see `developed`, otherwise its products
factor by factor from the right and its sums term by term, without developing them. Its flux is
not checked: see `operator_acted` on a charged system.
"""
function gathered_acted(system::System, g::IndexedOp{IsOperator}, x::AbstractVecOrMat)
    ts = developed(system, g)
    if !isnothing(ts)
        return dense_acted(Pure, system, ts, x)
    end
    c, s = scalarcoef(g), scalararg(g)
    if s isa ProdOp
        return c * foldr((f, y) -> gathered_acted(system, f, y), s.subs; init = x)
    elseif s isa SumOp
        return c * sum(t -> gathered_acted(system, t, x), s.subs)
    end
    error("$(obs_name(g)) neither develops nor is a sum or a product")
end
