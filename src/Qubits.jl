export Qubits

struct Qubit <: AbstractSite end

dim(::Qubit) = 2

controlled(op::ExprOp{Pure, N}) where N =
    TensorOp{Pure, N+1}(ExprOp[[Proj("Up")] ; [Id for _ in 1:N]]) + (Proj("Dn") ⊗ op)

@def_states(Qubit(),
[
    ["Up", "Z+"] => [1., 0.],
    ["Dn", "Z-"] => [0., 1.],
    ["+", "X+"] => [1., 1.] / √2,
    ["-", "X-"] => [1., -1.] / √2,
    ["i", "Y+"] => [1., im] / √2,
    ["-i", "Y-"] => [1., -im] / √2,
])


@def_operators(Qubit(),
[
    F = Id,
    X = [0. 1. ; 1. 0.],
    Y = [0. -im; im 0.],
    Z = [1. 0. ; 0. -1.],
    Sp = [0. 1. ; 0. 0.],
    Sm = [0. 0. ; 1. 0.],
    Sx = X / 2,
    Sy = Y / 2,
    Sz = Z / 2,
    S2 = 0.75 * Id,
    H = [1. 1. ; 1. -1] / √2,
    S = [1. 0. ; 0. im],
    Swap = (Id ⊗ Id + X ⊗ X + Y ⊗ Y + Z ⊗ Z) / 2
])

function graph_state(tp::PM, g::Vector{Tuple{Int, Int}}; limits::Limits=Limits(cutoff=1.e-16), kwargs...)
    n = graph_base_size(g)
    s = System(n, Qubit())
    state = State(tp, s, "+")
    gates = prod(CZ(i, j) for (i, j) in g)
    state = apply(gates, state; cutoff, kwargs...)
    return state
end

create_graph_state(tp::PM, g::Vector{Tuple{Int, Int}}; limits = Limits()) = 
    [
        CreateState(
            name = "Creating initial state |++...++> for graph state",
            type = tp,
            system = System(graph_base_size(g), Qubit()),
            state = "+",
        ),
        Gates(
            name = "Applying gates CZ for building graph state",
            gates = prod(CZ(i, j) for (i, j) in g),
            limits = limits
        )
    ]

module Qubits

import ..Qubit, ..controlled, ..graph_state, ..create_graph_state, ..X, ..Y, ..Z, ..Sx, ..Sy, ..Sz, ..S2, ..Sp, ..Sm, ..H, ..S, ..Swap
export Qubit, controlled, graph_state, create_graph_state, X, Y, Z, Sx, Sy, Sz, S2, Sp, Sm, H, S, Swap

export graph_state, create_graph_state

end