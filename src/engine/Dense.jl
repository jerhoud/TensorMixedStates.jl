# The terms of an operator, a channel or a generator acting on a dense state: the vector of a pure
# state, or the density matrix of a mixed one, on the basis of the product states ordered as
# `kron` orders them, the first site varying the slowest.

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
            isone(b.m) ? y : axis_acted(b.m, y, prod(dims[k+1:end]; init = 1), dims[k],
                                        prod(dims[1:k-1]; init = 1) * extra)
        end * laid.coef
    end
end

"""
    dense_acted(R, system, terms, x)

the terms `terms` of an operator applied to the vector `x` or to the kets of the matrix `x`, for
`Pure`, or those of a superoperator applied to the density matrix `x`, for `Mixed`, see `lower`
"""
function dense_acted(::Type{Pure}, system::System, ts::Vector{Term}, x::AbstractVecOrMat)
    dims = [ dim(s) for s in system.sites ]
    if size(x, 1) ≠ prod(dims)
        error("a dense state of the system has $(prod(dims)) rows, not $(size(x, 1))")
    end
    return reshape(terms_acted(Pure, system, ts, vec(x), dims, size(x, 2)), size(x))
end

function dense_acted(::Type{Mixed}, system::System, ts::Vector{Term}, x::AbstractMatrix)
    dims = [ dim(s) for s in system.sites ]
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
