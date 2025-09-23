export Qubits

"""
    type Qubit

A site type for representing qubit sites, that is a two level system.

# Example

    Qubit()

# States

- `"Up", "Z+", "↑", "0"`  : the up state
- `"Dn", "Z-", "↓", "1"`  : the down state
- `"+", "X+"`             : the + state (+1 eigenvector of X)
- `"-", "X-"`             : the - state (-1 eigenvector of X)
- `"i", "Y+"`             : the i state (+1 eigenvector of Y)
- `"-i", "Y-"`            : the -i state (-1 eigenvector of Y)

# Operators

- `X, Y, Z`          : the Pauli operators
- `Sp, Sm`           : the ``S^+`` and ``S^-`` operators
- `Sx, Sy, Sz, S2`   : the ``S_x``, ``S_y``, ``S_z`` operators (half the Pauli operators) and ``S^2``
- `H, S, T, Swap`    : the Hadamard, S, T and Swap gates
- `Phase(t)`         : the phase gate
- `controlled(gate)` : controlled gate
"""
struct Qubit <: AbstractSite end

dim(::Qubit) = 2

"""
    controlled(op)

the controlled gate constructor

# Examples
    CZ = controlled(Z)
    Toffoli = controlled(controlled(X))
"""
controlled(op::GenericOp{Pure, N}) where N = 
    Proj("Up") ⊗ Identity(op) + Proj("Dn") ⊗ op

@def_states(Qubit(),
[
    ["Up", "Z+", "↑"] => [1., 0.],
    ["Dn", "Z-", "↓"] => [0., 1.],
    ["+", "X+"] => [1., 1.] / √2,
    ["-", "X-"] => [1., -1.] / √2,
    ["i", "Y+"] => [1., im] / √2,
    ["-i", "Y-"] => [1., -im] / √2,
])

@def_operators(Qubit(),
[
    involution_op =>
    [
        X = [0. 1. ; 1. 0.],
        Y = [0. -im; im 0.],
        Z = [1. 0. ; 0. -1.],
        H = [1. 1. ; 1. -1] / √2,
        Swap = (Id ⊗ Id + X ⊗ X + Y ⊗ Y + Z ⊗ Z) / 2
    ],
    selfadjoint_op =>
    [
        Sx = X / 2,
        Sy = Y / 2,
        Sz = Z / 2,
        S2 = 0.75 * Id,
    ],
    plain_op =>
    [
        Sp = [0. 1. ; 0. 0.],
        Sm = dag(Sp),
        S = [1. 0. ; 0. im],
        T = [1. 0. ; 0. (1 + im)/√2],
    ]
])


"""
    Phase(t)

the phase gate for qubits
"""
Phase(t) = Operator("Phase($t)", [1. 0 ; 0 exp(im * t)], plain_op)

#= """
    graph_state(Pure()|Mixed(), graph::Vector{Tuple{Int, Int}}; limits)

create a graph state corresponding to the given graph

# Examples

    graph_state(Pure(), complete_graph(10); limits = Limits(maxdim = 10))
"""
function graph_state(tp::Repr, g::Vector{Tuple{Int, Int}}; limits::Limits=Limits(cutoff=1.e-16))
    n = graph_base_size(g)
    s = System(n, Qubit())
    state = State(tp, s, "+")
    CZ = controlled(Z)
    gates = prod(CZ(i, j) for (i, j) in g)
    state = apply(gates, state; limits)
    return state
end

"""
    create_graph_state(Pure()|Mixed(), graph::Vector{Tuple{Int, Int}}; limits)

create a phase for building a graph state to use in `SimData` and `runTMS`
"""
create_graph_state(tp::Repr, g::Vector{Tuple{Int, Int}}; kwargs...) = 
    [
        CreateState(
            name = "Creating initial state |++...++> for graph state",
            type = tp,
            system = System(graph_base_size(g), Qubit()),
            state = "+",
        ),
        Gates(;
            name = "Applying gates CZ for building graph state",
            gates = prod(controlled(Z)(i, j) for (i, j) in g),
            kwargs...
        )
    ]
 =#

module Qubits

import ..Qubit, ..controlled, ..graph_state, ..create_graph_state, ..X, ..Y, ..Z, ..Sx, ..Sy, ..Sz, ..S2, ..Sp, ..Sm, ..H, ..S, ..T, ..Swap, ..Phase
export Qubit, controlled, graph_state, create_graph_state, X, Y, Z, Sx, Sy, Sz, S2, Sp, Sm, H, S, T, Swap, Phase

export graph_state, create_graph_state

end