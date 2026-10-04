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

- `C` : the destruction operator, fermionic
- `N` : the number of fermions
- `F` : the Jordan-Wigner operator, ``(-1)^N``
"""
struct Fermion <: AbstractSite
    conserve::String
end

Fermion(; conserve = ()) = Fermion(conserve_string(Fermion(""), conserve))

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
    involution_op =>
    [
        F = Float64[1 0 ; 0 -1]
    ]
])

fermion_species(::Fermion) = (C,)

@create_site_module(Fermions, [Fermion, C, N])
