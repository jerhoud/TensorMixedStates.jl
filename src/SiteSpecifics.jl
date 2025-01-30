import ITensors: state
export graph_state, create_graph_state

# example on how to defined a mixed state
# state(::StateName"Half", ::SiteType"MixedQubit") = [0.5 0 0 0.5]

function graph_state(tp::Union{TPure, TMixed}, g::Vector{Tuple{Int, Int}}; cutoff=1e-16, kwargs...)
    n = graph_base_size(g)
    s = System(n, "Qubit")
    state = State(tp, s, "+")
    gates = prod(CZ(i, j) for (i, j) in g)
    state = apply(gates, state; cutoff, kwargs...)
    return state
end

create_graph_state(tp::Union{TPure, TMixed}, g::Vector{Tuple{Int, Int}}; limits = no_limits) = 
    [
        CreateState(
            name = "Creating initial state |++...++> for graph state",
            type = tp,
            system = System(graph_base_size(g), "Qubit"),
            state = "+",
        ),
        Gates(
            name = "Applying gates CZ for building graph state",
            gates = prod(CZ(i, j) for (i, j) in g),
            limits = limits
        )
    ]