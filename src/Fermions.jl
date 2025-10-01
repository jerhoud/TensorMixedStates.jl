export Fermions

"""
    Fermion()

a site type to represent fermion sites (dim is 2)

# Examples

    Fermion()

# States

- `"0", "Emp"` : empty state
- `"1", "Occ"` : occupied state

# Operators

- `C` : the destruction operator
- `A` : the Jordan-Wigner transform of C (...C = FFFA)
- `N` : the number of fermions operator
"""
struct Fermion <: AbstractSite end

dim(::Fermion) = 2

@def_states(Fermion(),
[
    "Emp" => [1., 0.],
    "Occ" => [0., 1.],
])

@def_operators(Fermion(),
[
    fermionic_op => 
    [
        C = [0. 1. ; 0. 0.],
    ],
    selfadjoint_op =>
    [
        N = dag(C) * C,
    ],
    plain_op =>
    [
        A = C,
    ],
    involution_op =>
    [
        F = Float64[1 0 ; 0 -1]
    ]
])

module Fermions

import ..Fermion, ..C, ..N, ..A
export Fermion, C, N, A

end
