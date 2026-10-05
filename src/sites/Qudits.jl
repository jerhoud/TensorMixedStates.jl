# The Qudit site type, a level system of any dimension with the generalized Pauli operators,
# with its states and operators and the gate Sumd, gathered in the module Qudits.

export Qudits

"""
    Qudit(dim; conserve = ())

the site type of a qudit, a `dim` level system carrying the generalized Pauli operators rather
than spin operators.

Writing ``d`` for the dimension and ``\\omega = e^{2i\\pi/d}``, the clock operator acts as
``Z_d|n\\rangle = \\omega^n |n\\rangle`` and the shift operator as
``X_d|n\\rangle = |n+1 \\mod d\\rangle``. Both are unitary of order ``d`` and obey the Weyl
relation ``Z_d X_d = \\omega X_d Z_d``.

`Qudit(2)` is the qubit: `Zd` and `Xd` are the Pauli `Z` and `X`, `Hd` and `S` the Hadamard and
phase gates. `Zd` and `Xd` are named apart because for ``d > 2`` they are not involutions.

`Hd`, `S` and `Sumd` generate the Clifford group, the unitaries mapping the operators generated
by `Zd` and `Xd` onto themselves by conjugation.

# Examples

    Qudit(3)
    Qudit(3, conserve = Zd)     # the Zd symmetry of the clock model

# States

`"0", "1", ..., "dim-1"` : the state of that level

# Operators

- `N`  : the level operator, `diag(0, 1, ..., d-1)`
- `Zd` : the clock operator
- `Xd` : the shift operator
- `Hd` : the Fourier operator, or generalized Hadamard gate
- `S`  : the phase operator, of order ``d`` for odd ``d`` and ``2d`` for even ``d``
- `Sumd(d)` : the generalized controlled not of two qudits, see `Sumd`
"""
struct Qudit <: AbstractSite
    dim::Int
    conserve::String
    Qudit(dim::Int, conserve::AbstractString) =
        if dim ≥ 1
            new(dim, conserve)
        else
            error("a Qudit has a dimension of 1 at least, not $dim")
        end
end

Qudit(dim::Int; conserve = ()) = Qudit(dim, conserve_string(Qudit(dim, ""), conserve))

dim(a::Qudit) = a.dim

@def_operators(Qudit(2),
[
    selfadjoint_op =>
    [
        N = s -> [ i == j ? i - 1. : 0. for i in 1:dim(s), j in 1:dim(s) ],
    ],
    plain_op =>
    [
        # a function of N rather than its matrix, so that as a conserved quantity it carries its
        # modulus: on a qubit the eigenvalues ±1 were read as integers, and their sum conserved
        Zd = s -> dim(s) == 1 ? ones(ComplexF64, 1, 1) : mod(N, dim(s)),
        Xd = s -> [ mod(i - j - 1, dim(s)) == 0 ? 1.0 + 0im : 0.0im
                    for i in 1:dim(s), j in 1:dim(s) ],
        Hd = s -> [ exp(2im * π * (i - 1) * (j - 1) / dim(s)) / √dim(s)
                    for i in 1:dim(s), j in 1:dim(s) ],
        S = s -> begin
            d = dim(s)
            τ = exp(1im * π * (d^2 + 1) / d)      # a square root of ω
            [ i == j ? τ^((i - 1)^2) : 0.0im for i in 1:d, j in 1:d ]
        end,
    ]
])

"""
    Sumd(d)

the generalized controlled not of two qudits of dimension `d`,
``|i, j\\rangle \\mapsto |i, i + j \\mod d\\rangle``; `Sumd(2)` is the CNOT. It is built as
the expression ``\\sum_i P_i \\otimes X_d^i`` rather than as a matrix, so that it can be used
in a hamiltonian as well as applied as a gate.

# Examples

    Sumd(3)(1, 2)
"""
Sumd(d::Int) =
    Operator{2}("Sumd($d)", sum(Proj(i) ⊗ Xd^i for i in 0:d-1), plain_op)

@create_site_module(Qudits, [Qudit, N, Zd, Xd, Hd, S, Sumd])
