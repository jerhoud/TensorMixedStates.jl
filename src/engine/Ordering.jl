# The Fenwick tree of counts and the binary heap with which the gathering of a product orders
# its factors, see `sort_factors!`.

"""
    count_add!(tree, k, x)

add `x` to the count of position `k` of `tree`, a Fenwick tree of counts
"""
function count_add!(tree::Vector{Int}, k::Int, x::Int)
    while k ≤ length(tree)
        tree[k] += x
        k += k & -k
    end
    return tree
end

"""
    count_to(tree, k)

the sum of the counts of the positions 1 to `k` of the Fenwick tree `tree`
"""
function count_to(tree::Vector{Int}, k::Int)
    c = 0
    while k > 0
        c += tree[k]
        k -= k & -k
    end
    return c
end

"""
    heap_push!(h, x, less)
    heap_pop!(h, less)

add `x` to the binary heap `h` ordered by `less`, or take its least element out
"""
function heap_push!(h::Vector, x, less)
    push!(h, x)
    i = length(h)
    while i > 1 && less(h[i], h[i >> 1])
        h[i], h[i >> 1] = h[i >> 1], h[i]
        i >>= 1
    end
    return h
end

function heap_pop!(h::Vector, less)
    top, x = h[1], pop!(h)
    if !isempty(h)
        h[1] = x
        i = 1
        while true
            m = i
            for c in (2i, 2i + 1)
                if c ≤ length(h) && less(h[c], h[m])
                    m = c
                end
            end
            if m == i
                break
            end
            h[i], h[m] = h[m], h[i]
            i = m
        end
    end
    return top
end
