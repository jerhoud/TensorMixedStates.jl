# Graphs as lists of edges, line, circle, complete or square lattice, for building hamiltonians
# and graph states.

export graph_base_size, line_graph, circle_graph, complete_graph, square_lattice

"""
    graph_base_size(::Vector{Tuple{Int, Int}})

the largest vertex of a graph, which is the number of sites a system needs to host it. It is
neither the size of the graph, its number of edges, nor its order, its number of vertices: a
vertex may be in no edge, and the system still has to hold that site.

# Examples

    graph_base_size([(1, 2), (2, 5)])     # 5
"""
graph_base_size(g::Vector{Tuple{Int, Int}}) = maximum(maximum, g)

"""
    line_graph(n)

the chain of `n` vertices, as a vector of edges: `[(1, 2), (2, 3), ..., (n-1, n)]`.
"""
line_graph(n::Int) =
    [(i, i+1) for i in 1:(n-1)]

"""
    circle_graph(n)

the ring of `n` vertices, as a vector of edges: `[(1, 2), (2, 3), ..., (n-1, n), (n, 1)]`.
The ring of two vertices joins them twice, `[(1, 2), (2, 1)]`, as the two bonds of a periodic
chain of two sites; a ring of fewer vertices is refused.
"""
function circle_graph(n::Int)
    if n < 2
        error("a ring has two vertices at least, and circle_graph was given $n")
    end
    return [line_graph(n); [(n, 1)]]
end

"""
    complete_graph(n)

the complete graph of `n` vertices, as a vector of edges:
`[(1, 2), (1, 3), ..., (1, n), (2, 3), ..., (n-1, n)]`.

# Examples

    graph_state(complete_graph(10); limits = Limits(maxdim = 10))
"""
complete_graph(n::Int) =
    [(i, j) for i in 1:(n-1) for j in (i+1):n]

"""
    square_lattice(nx, ny)
    square_lattice(n)

the `nx` by `ny` square lattice, as a vector of edges, its sites numbered along a snake
running down the columns: the first column top to bottom, the second bottom to top, and so on.
The chain thus runs through the lattice column by column, in a band `ny` sites wide: vertical
bonds join consecutive sites and horizontal ones are at most `2ny - 1` sites apart. Put the short side in `ny`, since it
bounds the bond dimension. With a single argument the lattice is `n` by `n`.

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
