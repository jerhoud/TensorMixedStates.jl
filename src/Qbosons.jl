export Qbosons

"""
    Qboson(q, dim)

a site type to represent q-boson sites, it is parametred by `q` and the dimension of the Hilbert space
(maximum occupancy is `dim - 1`).

``a|n\\rangle = \\sqrt{1-q^n} |n-1\\rangle`` and ``a^\\dagger |n\\rangle = \\sqrt{1 - q^{n+1}} |n+1\\rangle``


# Examples

    Qboson(0.1, 4)

# States

`"0", "1", ...`

# Operators

- `A` : the destruction operator
- `N` : the number of q-bosons operator
- `Q`, `P` : the position and momentum quadratures
"""
struct Qboson <: AbstractSite
    q::Float64
    dim::Int
end

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
