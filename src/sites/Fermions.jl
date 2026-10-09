# The Fermion site type, a spinless fermion mode, with its states and operators, gathered in the
# module Fermions.

export Fermions

"""
    Fermion(; conserve = ())

the site type of a spinless fermion mode, of dimension 2.

# Examples

    Fermion()
    Fermion(conserve = N)
    Fermion(conserve = parity(N))
    Fermion(conserve = strong(N))

# States

- `"Emp", "0"` : the empty state
- `"Occ", "1"` : the occupied state

# Operators

- `C` : the annihilation operator
- `N` : the number of fermions
- `F` : the fermionic parity, ``(-1)^N``
"""
struct Fermion <: AbstractSite
    conserve::Conservation
end

Fermion(; conserve = ()) = conserving(Fermion; conserve)

dim(::Fermion) = 2

fermion_parities(::Fermion) = (0, 1)

@def_states Fermion [
    "Emp" => [1., 0.],
    "Occ" => [0., 1.],
]

@def_operators Fermion [
    C = [0. 1. ; 0. 0.],
    N = dag(C) * C,
]

fermion_species(::Fermion) = (C,)

@create_site_module(Fermions, [Fermion, C, N])
