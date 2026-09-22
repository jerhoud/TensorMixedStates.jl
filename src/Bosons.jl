export Bosons

"""
    Boson(dim)

a site type to represent boson sites, it is parametrised by the dimension of the Hilbert space
(maximum occupancy is `dim - 1`)

# Examples

    Boson(4)
    Boson(4, conserve = N)

# States

`"0", "1", ...`

# Operators

- `A` : the destruction operator
- `N` : the number of bosons operator
- `Q`, `P` : the position and momentum quadratures
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
