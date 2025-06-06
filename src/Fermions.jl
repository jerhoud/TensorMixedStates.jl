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
    C = [0. 1. ; 0. 0.],
],
true)

@def_operators(Fermion(),
[
    F = [1. 0. ; 0. -1.],
    A = C,
    N = dag(C)*C
])

module Fermions

import ..Fermion, ..C, ..N, ..A
export Fermion, C, N, A

end
