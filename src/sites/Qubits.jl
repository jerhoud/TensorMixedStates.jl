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

@create_site_module(Qubits, [Qubit, controlled, X, Y, Z, Sx, Sy, Sz, S2, N, Sp, Sm, H, S, T,
                             Phase])
