# The Qubit site type, its states and operators, the controlled gates and the graph states,
# gathered in the module Qubits.

export Qubits

"""
    Qubit(; conserve = ())

the site type of a qubit, a two level system, whose basis is `"Up"`, `"Dn"`.

# Examples

    Qubit()
    Qubit(conserve = N)      # the number of excitations
    Qubit(conserve = 2Sz)    # the same conservation, in spin language

# States

- `"Up", "Z+", "↑", "0"` : the up state
- `"Dn", "Z-", "↓", "1"` : the down state
- `"+", "X+"`            : the +1 eigenvector of `X`
- `"-", "X-"`            : the -1 eigenvector of `X`
- `"i", "Y+"`            : the +1 eigenvector of `Y`
- `"-i", "Y-"`           : the -1 eigenvector of `Y`

# Operators

- `X, Y, Z`          : the Pauli operators
- `Sp, Sm`           : the ``S^+`` and ``S^-`` operators
- `Sx, Sy, Sz, S2`   : the ``S_x``, ``S_y``, ``S_z`` operators, half the Pauli operators, and
                       ``S^2``
- `N`                : the number of excitations, ``(1 - Z)/2``: `"Dn"` is the occupied state
                       and `Sm` creates an excitation. `Qubit(conserve = N)` conserves what
                       `Qubit(conserve = 2Sz)` does, with charges 0 and 1 rather than ±1
- `H, S, T, Swap`    : the Hadamard, S, T and Swap gates
- `Phase(t)`         : the phase gate
- `controlled(gate)` : the controlled gate
"""
struct Qubit <: AbstractSite
    conserve::String
end

Qubit(; conserve = ()) = Qubit(conserve_string(Qubit(""), conserve))

dim(::Qubit) = 2

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
    ],
    selfadjoint_op =>
    [
        Sx = X / 2,
        Sy = Y / 2,
        Sz = Z / 2,
        S2 = 0.75 * Id,
        N = [0. 0. ; 0. 1.],
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
    controlled_name(op)

the default name of `controlled(op)`: `C` followed by the name of an `Operator`, as `CZ`, and
`controlled(op)` for an expression.
"""
controlled_name(a::Operator) = "C" * a.name
controlled_name(a::Op) = "controlled($a)"

"""
    controlled_type(op)

the default `OpType` of `controlled(op)`: that of an `Operator`, and `plain_op` for an
expression or a fermionic `Operator`, no fermionic `Operator` acting on several sites
"""
controlled_type(a::Operator) = a.type == fermionic_op ? plain_op : a.type
controlled_type(a::Op) = plain_op

"""
    controlled(op; name, type)

the gate applying `op` to the following sites when the first one, a qubit, is in the state
`"1"`, and nothing otherwise.

- `name`: its name (default `"C"` followed by the name of a named operator, `CZ` for
  `controlled(Z)`, and `"controlled(op)"` otherwise)
- `type`: its `OpType` (default that of a named operator other than fermionic, and `plain_op`
  otherwise)

# Examples

    CZ = controlled(Z)
    Toffoli = controlled(controlled(X))
    CZ(1, 2)
"""
controlled(op::GenericOp{Pure, N}; name::String = controlled_name(op), type = controlled_type(op)) where N =
    Operator{N+1}(name, Proj(0) ⊗ IdentityOp(op) + Proj(1) ⊗ op, type)

# factors of a definite charge, which X ⊗ X + Y ⊗ Y lacks, and dag(Sp) rather than Sm, which
# `measure` does not know to be its adjoint when it tests for a real value
"""
    Swap

the gate exchanging the states of two qubits.
"""
const Swap = Operator{2}("Swap", (Id ⊗ Id + Z ⊗ Z) / 2 + Sp ⊗ dag(Sp) + dag(Sp) ⊗ Sp, involution_op)

"""
    Phase(t)

the phase gate of a qubit, ``\\mathrm{diag}(1, e^{it})``.
"""
Phase(t) = Operator{1}("Phase($t)", [1. 0 ; 0 exp(im * t)], plain_op)

"""
    graph_state(graph; limits = Limits(cutoff = 1e-16))

the graph state of `graph`, a list of edges: the pure state of `graph_base_size(graph)` qubits
all in `"+"`, to which `controlled(Z)` is applied on every edge, truncated by `limits`.

# Examples

    graph_state(complete_graph(10); limits = Limits(maxdim = 10))
"""
function graph_state(g::Vector{Tuple{Int, Int}}; limits::Limits=Limits(cutoff=1.e-16))
    n = graph_base_size(g)
    s = System(n, Qubit())
    state = State{Pure}(s, "+")
    CZ = controlled(Z)
    gates = prod(CZ(i, j) for (i, j) in g)
    state = apply(gates, state; limits)
    return state
end

"""
    create_graph_state(graph; kwargs...)

the phases building the graph state of `graph`, see `graph_state`: a `CreateState` of qubits
all in `"+"`, then a `Gates` phase applying `controlled(Z)` on every edge, which receives the
keyword arguments. `SimData` accepts this list wherever a phase is expected.

# Examples

    create_graph_state(complete_graph(10); limits = Limits(cutoff = 1e-14))
"""
create_graph_state(g::Vector{Tuple{Int, Int}}; kwargs...) = 
    [
        CreateState(
            type = Pure(),
            name = "Creating initial state |++...++> for graph state",
            system = System(graph_base_size(g), Qubit()),
            state = "+",
        ),
        Gates(;
            name = "Applying gates CZ for building graph state",
            gates = prod(controlled(Z)(i, j) for (i, j) in g),
            kwargs...
        )
    ]

@create_site_module(Qubits, [Qubit, controlled, graph_state, create_graph_state, X, Y, Z, Sx, Sy, Sz, S2, N, Sp, Sm, H, S, T, Swap, Phase])