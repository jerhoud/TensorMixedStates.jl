# The Qubit site type, its states and operators and the controlled gates, gathered in the
# module Qubits.

export Qubits

"""
    Qubit(; conserve = ())

the site type of a qubit, a two level system, whose basis is `"Up"`, `"Dn"`. `Z` is a charge of
values ±1, so `conserve = Z` conserves the magnetization, as `2Sz` does; its parity is
`parity(N)`.

# Examples

    Qubit()
    Qubit(conserve = N)      # the number of excitations
    Qubit(conserve = 2Sz)    # the same conservation, in spin language
    Qubit(conserve = parity(N))

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
                       and `Sm` creates an excitation
- `H, S, T`          : the Hadamard, S and T gates
- `Phase(t)`         : the phase gate
- `controlled(gate)` : the controlled gate
"""
struct Qubit <: AbstractSite
    conserve::Conservation
end

Qubit(; conserve = ()) = conserving(Qubit; conserve)

dim(::Qubit) = 2

@def_states Qubit [
    ["Up", "Z+", "↑"] => [1., 0.],
    ["Dn", "Z-", "↓"] => [0., 1.],
    ["+", "X+"] => [1., 1.] / √2,
    ["-", "X-"] => [1., -1.] / √2,
    ["i", "Y+"] => [1., im] / √2,
    ["-i", "Y-"] => [1., -im] / √2,
]

@def_operators Qubit [
    X = [0. 1. ; 1. 0.],
    Y = [0. -im ; im 0.],
    Z = [1. 0. ; 0. -1.],
    H = [1. 1. ; 1. -1] / √2,
    Sz = Z / 2,
    N = [0. 0. ; 0. 1.],
    Sp = [0. 1. ; 0. 0.],
    S = [1. 0. ; 0. im],
    T = [1. 0. ; 0. (1 + im) / √2],
]

"""
    controlled_name(op)

the default name of `controlled(op)`: `C` followed by the name of an `Operator`, as `CZ`, and
`controlled(op)` for an expression
"""
controlled_name(a::Operator) = "C" * a.name
controlled_name(a::Op) = "controlled($a)"

"""
    controlled(op; name)

the gate applying `op` to the following sites when the first one, a qubit, is in the state
`"1"`, and nothing otherwise. Its name is `"C"` followed by the name of a named operator,
`CZ` for `controlled(Z)`, and `"controlled(op)"` otherwise.

# Examples

    CZ = controlled(Z)
    Toffoli = controlled(controlled(X))
    CZ(1, 2)
"""
controlled(op::GenericOp{IsOperator, N}; name::String = controlled_name(op)) where N =
    Operator{N+1}(name, Proj(Level => 0) ⊗ IdentityOp(op) + Proj(Level => 1) ⊗ op)

"""
    Phase(t)

the phase gate of a qubit, ``\\mathrm{diag}(1, e^{it})``.
"""
Phase(t) = Operator{1}("Phase($t)", [1. 0 ; 0 exp(im * t)])

"""
    graph_state(graph; limits = Limits())

the graph state of `graph`, a list of edges: the pure state of `graph_base_size(graph)` qubits
all in `"+"`, to which `controlled(Z)` is applied on every edge, truncated by `limits`.

# Examples

    graph_state(complete_graph(10); limits = Limits(maxdim = 10))
"""
function graph_state(g::Vector{Tuple{Int, Int}}; limits::Limits = Limits())
    if isempty(g)
        error("graph_state needs a graph with edges, which give its qubits")
    end
    state = State{Pure}(System(graph_base_size(g), Qubit()), "+")
    CZ = controlled(Z)
    return apply(prod(CZ(i, j) for (i, j) in g), state; limits)
end

"""
    thermal_relaxation_check(T1, T2)

refuse times `T1` and `T2` that are not positive, or a `T2` above `2T1`, which would ask for a
negative rate of dephasing
"""
function thermal_relaxation_check(T1::Real, T2::Real)
    if !(T1 > 0 && T2 > 0)
        error("the times T1 = $T1 and T2 = $T2 of a thermal relaxation must be positive")
    elseif T2 > 2T1
        error("T2 = $T2 is above 2T1 = $(2T1), which no relaxation reaches: the decay alone " *
              "takes the coherences in a time 2T1")
    end
end

"""
    thermal_relaxation_dissipator(T1, T2)

the generator of the relaxation of a qubit of times `T1` and `T2`: its decay to `"Up"`, of time
`T1`, and the dephasing that brings the decay of its coherences to the time `T2`,
`amplitude_damping_dissipator(1 / T1) + dephasing_dissipator(1 / T2 - 1 / (2T1))`. `T2` cannot
exceed `2T1`. Evolving under it for a time `t` is `thermal_relaxation_gate(T1, T2, t)`.

# Examples

    -im * H + sum(thermal_relaxation_dissipator(50., 30.)(i) for i in 1:10)
"""
function thermal_relaxation_dissipator(T1::Real, T2::Real)
    thermal_relaxation_check(T1, T2)
    return amplitude_damping_dissipator(1 / T1) + dephasing_dissipator(1 / T2 - 1 / (2T1))
end

"""
    thermal_relaxation_gate(T1, T2, t)

the channel of the relaxation of a qubit of times `T1` and `T2` over a time `t`,
`exp(t * thermal_relaxation_dissipator(T1, T2))`: the population of `"Dn"` decays by
``e^{-t/T_1}`` and the coherences by ``e^{-t/T_2}``. `T2` cannot exceed `2T1`.

# Examples

    prod(thermal_relaxation_gate(50., 30., 0.1)(i) for i in 1:10)
"""
function thermal_relaxation_gate(T1::Real, T2::Real, t::Real)
    if t < 0
        error("thermal_relaxation_gate takes a time of zero or above, not $t")
    end
    return exp(t * thermal_relaxation_dissipator(T1, T2))
end

@create_site_module(Qubits, [Qubit, controlled, graph_state, X, Y, Z, Sx, Sy, Sz, S2, N, Sp, Sm,
                             H, S, T, Phase, thermal_relaxation_dissipator,
                             thermal_relaxation_gate])
