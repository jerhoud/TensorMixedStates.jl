export graph_base_size, line_graph, circle_graph, complete_graph

"""
    graph_base_size(::Vector{Tuple{Int, Int}})

return the maximum vertex number in a graph
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

return the complete graph with n vertices 1-2, 1-3, 2-3, 1-4, 2-4, 3-4 ...
"""
complete_graph(n::Int) =
    [(i, j) for i in 1:(n-1) for j in (i+1):n]
