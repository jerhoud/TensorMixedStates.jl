# Graphs as lists of edges, line, circle, complete or square lattice, for building Hamiltonians
# and graph states.

export graph_base_size, line_graph, circle_graph, complete_graph, square_lattice

"""
    graph_base_size(::Vector{Tuple{Int, Int}})

the largest vertex of a graph, the number of sites a system needs to host it, 0 for a graph
without edges. It is not the number of vertices: a vertex may be in no edge.

# Examples

    graph_base_size([(1, 2), (2, 5)])     # 5
"""
graph_base_size(g::Vector{Tuple{Int, Int}}) = maximum(maximum, g; init = 0)

"""
    line_graph(n)

the chain of `n` vertices, as a vector of edges: `[(1, 2), (2, 3), ..., (n-1, n)]`.
"""
line_graph(n::Int) =
    [(i, i+1) for i in 1:(n-1)]

"""
    check_ring(n, what)

refuse a ring of fewer than two vertices, `what` naming the function given it
"""
function check_ring(n::Int, what)
    if n < 2
        error("a ring has two vertices at least, and $what was given $n")
    end
    return nothing
end

"""
    circle_graph(n)

the ring of `n` vertices, as a vector of edges: `[(1, 2), (2, 3), ..., (n-1, n), (n, 1)]`.
The ring of two vertices joins them twice, `[(1, 2), (2, 1)]`, as the two bonds of a periodic
chain of two sites; a ring of fewer vertices is refused.
"""
function circle_graph(n::Int)
    check_ring(n, "circle_graph")
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
    square_lattice(nx, ny; xperiodic = false, yperiodic = false)
    square_lattice(n; xperiodic = false, yperiodic = false)

the `nx` by `ny` square lattice, as a vector of edges, its sites numbered along a snake
running down the columns: the first column top to bottom, the second bottom to top, and so on.
Horizontal bonds join sites at most `2ny - 1` apart: put the short side in `ny`, since it
bounds the bond dimension. With a single argument the lattice is `n` by `n`.

`yperiodic` joins the last row to the first, making a cylinder, its bonds within the same bound;
`xperiodic` joins the last column to the first, by bonds across the whole lattice, which cost
far more: put a single periodic direction in `y`. The bonds closing a direction come after the
others, and, as in `circle_graph`, a direction of length 2 has its bonds twice, and one of
length 1 cannot be periodic.

# Examples

    square_lattice(3, 2)                    # 1-2, 3-4, 5-6, 1-4, 2-3, 4-5, 3-6
    square_lattice(4, 3; yperiodic = true)  # a cylinder of circumference 3

that is, for the first, the lattice

    1 4 5
    2 3 6
"""
function square_lattice(nx::Int, ny::Int; xperiodic::Bool = false, yperiodic::Bool = false)
    # site (r, c) of the snake, columns numbered from 1
    pos(r, c) = (c - 1) * ny + (isodd(c) ? r : ny - r + 1)
    g = Tuple{Int, Int}[]
    for c in 1:nx, r in 1:ny-1
        push!(g, minmax(pos(r, c), pos(r + 1, c)))
    end
    for c in 1:nx-1, r in 1:ny
        push!(g, minmax(pos(r, c), pos(r, c + 1)))
    end
    if yperiodic
        check_ring(ny, "square_lattice with yperiodic")
        append!(g, [ minmax(pos(ny, c), pos(1, c)) for c in 1:nx ])
    end
    if xperiodic
        check_ring(nx, "square_lattice with xperiodic")
        append!(g, [ minmax(pos(r, nx), pos(r, 1)) for r in 1:ny ])
    end
    return g
end

square_lattice(n::Int; xperiodic::Bool = false, yperiodic::Bool = false) =
    square_lattice(n, n; xperiodic, yperiodic)
