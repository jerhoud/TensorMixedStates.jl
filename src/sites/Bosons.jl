# The Boson site type, a boson mode truncated to a given number of levels, with its states and
# operators, gathered in the module Bosons.

export Bosons

"""
    Boson(dim; conserve = ())

the site type of a boson mode truncated to `dim` levels, the occupation going from 0 to
`dim - 1`.

# Examples

    Boson(4)
    Boson(4, conserve = N)
    Boson(4, conserve = parity(N))

# States

`"0", "1", ..., "dim-1"` : the state of that occupation

# Operators

- `A` : the destruction operator
- `N` : the number of bosons
- `Q`, `P` : the position and momentum quadratures, ``(A + A^\\dagger)/\\sqrt{2}`` and
  ``(A - A^\\dagger)/(i\\sqrt{2})``
"""
struct Boson <: AbstractSite
    dim::Int
    conserve::String
end

Boson(dim::Int; conserve = ()) = Boson(dim, conserve_string(Boson(dim, ""), conserve))

dim(a::Boson) = a.dim


@def_operators(Boson(2),
[
    selfadjoint_op =>
    [
        N = s -> [ i==j ? i - 1. : 0. for i in 1:dim(s), j in 1:dim(s) ]
    ],
    plain_op =>
    [
        A = s -> [ i==j-1 ? sqrt(i) : 0. for i in 1:dim(s), j in 1:dim(s) ]
    ],
    selfadjoint_op =>
    [
        Q = (A + dag(A)) / √2,
        P = (A - dag(A)) / (im * √2),
    ]
])

@create_site_module(Bosons, [Boson, N, A, Q, P])
