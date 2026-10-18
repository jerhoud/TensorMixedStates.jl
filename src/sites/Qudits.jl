# The Qudit site type, a level system of any dimension with the generalized Pauli operators,
# with its operators, gathered in the module Qudits, and its gates H, S and Sum, which the module
# Circuits gathers.

export Qudits

"""
    Qudit(dim; conserve = ())

the site type of a qudit, a `dim` level system carrying the generalized Pauli operators.

Writing ``d`` for the dimension and ``\\omega = e^{2i\\pi/d}``, the clock operator acts as
``Z|n\\rangle = \\omega^n |n\\rangle`` and the shift operator as
``X|n\\rangle = |n+1 \\mod d\\rangle``. `H`, `S` and `Sum` generate the Clifford group. `X`,
`Z`, `H` and `S` are those of `Qubit` for ``d = 2``.

# Examples

    Qudit(3)
    Qudit(3, conserve = Z)      # the Z3 symmetry of the clock model

# States

`"0", "1", ..., "dim-1"` : the state of that level

# Operators

- `N`  : the level operator, `diag(0, 1, ..., d-1)`
- `Z`  : the clock operator, `mod(N, d)`
- `X`  : the shift operator, ``ZX = \\omega XZ``

In the module `Circuits`:

- `H`  : the Fourier operator, or generalized Hadamard gate
- `S`  : the phase operator, of order ``d`` for odd ``d`` and ``2d`` for even ``d``
- `Sum` : the generalized controlled not of two qudits, see `Sum`
"""
struct Qudit <: AbstractSite
    dim::Int
    conserve::Conservation
    Qudit(dim::Int, conserve::Conservation) =
        if dim ≥ 1
            new(dim, conserve)
        else
            error("a Qudit has a dimension of 1 at least, not $dim")
        end
end

Qudit(dim::Int; conserve = ()) = conserving(Qudit, dim; conserve)

dim(a::Qudit) = a.dim

@def_operators Qudit [
    N = s -> [ i == j ? i - 1. : 0. for i in 1:dim(s), j in 1:dim(s) ],
    # a function of N, so that as a conserved quantity it carries its modulus
    Z = s -> dim(s) == 1 ? ones(ComplexF64, 1, 1) : mod(N, dim(s)),
    X = s -> [ mod(i - j - 1, dim(s)) == 0 ? 1.0 + 0im : 0.0im
               for i in 1:dim(s), j in 1:dim(s) ],
    H = s -> [ exp(2im * π * (i - 1) * (j - 1) / dim(s)) / √dim(s)
               for i in 1:dim(s), j in 1:dim(s) ],
    S = s -> begin
        d = dim(s)
        τ = exp(1im * π * (d^2 + 1) / d)      # a square root of ω
        [ i == j ? τ^((i - 1)^2) : 0.0im for i in 1:d, j in 1:d ]
    end,
]

"""
    Sum

the generalized controlled not of two qudits, ``|i, j\\rangle \\mapsto |i, i + j \\mod d\\rangle``,
``d`` being the dimension of the second, read on the sites it is placed on; on two `Qudit(2)`
it is the CNOT.

# Examples

    Sum(1, 2)
"""
const Sum = Operator{2}("Sum", (a, _) -> sum(Proj(Level => i) ⊗ X^i for i in 0:dim(a)-1))

@create_site_module(Qudits, [Qudit, N, Z, X])
