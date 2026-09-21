export graph_base_size, line_graph, circle_graph, complete_graph, square_lattice

"""
    graph_base_size(::Vector{Tuple{Int, Int}})

return the largest vertex number of a graph, that is the number of sites a system needs to
host it. Note that this is neither the size of the graph, which is its number of edges, nor
its order, which is its number of vertices: a graph may leave a vertex out of every edge,
and the system still has to hold that site.
"""
graph_base_size(g::Vector{Tuple{Int, Int}}) = maximum(maximum, g)

"""
    line_graph(n)

return the graph 1-2, 2-3, ..., (n-1)-n
"""
line_graph(n::Int) =
    [(i, i+1) for i in 1:(n-1)]

"""
    circle_graph(n)

return the graph 1-2, 2-3, ..., (n-1)-n, n-1
"""
circle_graph(n::Int) =
    [line_graph(n); [(n, 1)]]

"""
    complete_graph(n)

return the complete graph with n vertices 1-2, 1-3, ..., 1-n, 2-3, ..., 2-n, ...
"""
complete_graph(n::Int) =
    [(i, j) for i in 1:(n-1) for j in (i+1):n]

"""
    square_lattice(nx, ny)
    square_lattice(n)

return the `nx` by `ny` square lattice, its sites numbered along a snake running down the
columns: the first column is numbered top to bottom, the second one bottom to top, and so
on. This is the usual ordering for a tensor network, because the chain then advances along
`x` while staying only `ny` sites wide, like a thickened line: vertical bonds join
consecutive sites and horizontal ones are never further apart than `2ny - 1`. Put the
short side of the lattice in `ny`, since that is what bounds the bond dimension.

With a single argument the lattice is `n` by `n`.

# Examples

    square_lattice(3, 2)     # 1-2, 3-4, 5-6, 1-4, 2-3, 4-5, 3-6

that is the lattice

    1 4 5
    2 3 6
"""
function square_lattice(nx::Int, ny::Int)
    # site (r, c) of the snake, columns numbered from 1
    pos(r, c) = (c - 1) * ny + (isodd(c) ? r : ny - r + 1)
    g = Tuple{Int, Int}[]
    for c in 1:nx, r in 1:ny-1
        push!(g, minmax(pos(r, c), pos(r + 1, c)))
    end
    for c in 1:nx-1, r in 1:ny
        push!(g, minmax(pos(r, c), pos(r, c + 1)))
    end
    return g
end

square_lattice(n::Int) = square_lattice(n, n)
