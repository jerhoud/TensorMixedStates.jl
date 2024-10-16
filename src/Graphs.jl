export graph_base_size, line_graph, circle_graph, complete_graph

graph_base_size(g::Vector{Tuple{Int, Int}}) = maximum(maximum, g)

line_graph(n::Int) =
    [(i, i+1) for i in 1:(n-1)]

circle_graph(n::Int) =
    vcat(line_graph(n),[(n, 1)])

complete_graph(n::Int) =
    [(i, j) for i in 1:(n-1) for j in (i+1):n]
