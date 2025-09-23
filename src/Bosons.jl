export Bosons

"""
    Boson(dim)

a site type to represent boson sites, it is parametred by the dimension of the Hilbert space
(maximum occupancy is `dim - 1`)

# Examples

    Boson(4)

# States

`"0", "1", ...`

# Operators

- `A` : the destruction operator
- `N` : the number of bosons operator
"""
struct Boson <: AbstractSite
    dim::Int
end

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
    ]
])

module Bosons

import ..Boson, ..N, ..A
export Boson, N, A

end
