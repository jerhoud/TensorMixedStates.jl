# The Fenwick tree of counts and the binary heap with which the gathering of a product orders
# its factors and finds the one a factor merges into, see `sort_factors!` and `merged_pass`.

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
    last_counted(tree)

the last position of the Fenwick tree `tree` of counts 0 or 1 whose count is 1, 0 when there is
none
"""
function last_counted(tree::Vector{Int})
    c = count_to(tree, length(tree))
    if c == 0
        return 0
    end
    # the last position whose sum up to it is under c
    k, step = 0, prevpow(2, length(tree))
    while step > 0
        if k + step ≤ length(tree) && tree[k + step] < c
            k += step
            c -= tree[k]
        end
        step >>= 1
    end
    return k + 1
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
