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
"""
struct Qboson <: AbstractSite
    q::Float64
    dim::Int
end

dim(a::Qboson) = a.dim

@def_operators(Qboson(1., 2),
[
    F = Id,
], involution_op)

@def_operators(Qboson(1., 2),
[
    N = s -> [ i==j ? i - 1. : 0. for i in 1:dim(s), j in 1:dim(s) ]
], selfadjoint_op)

@def_operators(Qboson(1., 2),
[
    A = s -> [ i==j-1 ? sqrt(1-s.q^i) : 0. for i in 1:dim(s), j in 1:dim(s) ],
])

module Qbosons

import ..Qboson, ..N, ..A
export Qboson, N, A

end
