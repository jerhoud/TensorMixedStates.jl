# The Qboson site type, a q-deformed boson mode truncated to a given number of levels, with its
# states and operators, gathered in the module Qbosons.

export Qbosons

"""
    Qboson(q, dim; conserve = ())

the site type of a q-boson mode truncated to `dim` levels, the occupation going from 0 to
`dim - 1`, on which
``a|n\\rangle = \\sqrt{1-q^n} |n-1\\rangle`` and
``a^\\dagger |n\\rangle = \\sqrt{1 - q^{n+1}} |n+1\\rangle``.

# Examples

    Qboson(0.1, 4)
    Qboson(0.1, 4, conserve = N)

# States

`"0", "1", ..., "dim-1"` : the state of that occupation

# Operators

- `A` : the destruction operator ``a``
- `N` : the number of q-bosons
- `Q`, `P` : the position and momentum quadratures, ``(a + a^\\dagger)/\\sqrt{2}`` and
  ``(a - a^\\dagger)/(i\\sqrt{2})``
"""
struct Qboson <: AbstractSite
    q::Float64
    dim::Int
    conserve::String
end

Qboson(q::Real, dim::Int; conserve = ()) =
    Qboson(q, dim, conserve_string(Qboson(q, dim, ""), conserve))

dim(a::Qboson) = a.dim

@def_operators(Qboson(1., 2),
[
    selfadjoint_op =>
    [
        N = s -> [ i==j ? i - 1. : 0. for i in 1:dim(s), j in 1:dim(s) ]
    ],
    plain_op =>
    [
        A = s -> [ i==j-1 ? sqrt(1-s.q^i) : 0. for i in 1:dim(s), j in 1:dim(s) ]
    ],
    selfadjoint_op =>
    [
        Q = (A + dag(A)) / √2,
        P = (A - dag(A)) / (im * √2),
    ]
])

@create_site_module(Qbosons, [Qboson, N, A, Q, P])
