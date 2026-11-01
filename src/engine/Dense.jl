# The terms of an operator, a channel or a generator acting on a dense state: the vector of a pure
# state, or the density matrix of a mixed one, on the basis of the product states ordered as
# `kron` orders them, the first site varying the slowest.

"""
the array `x`, seen as `(a, d, c)`, with the matrix `m` applied to its middle axis, as a vector
"""
function axis_acted(m::AbstractMatrix, x::AbstractArray, a::Int, d::Int, c::Int)
    v = reshape(permutedims(reshape(x, a, d, c), (2, 1, 3)), d, :)
    return vec(permutedims(reshape(m * v, d, a, c), (2, 1, 3)))
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
            axis_acted(b.m, y, prod(dims[k+1:end]; init = 1), dims[k],
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
