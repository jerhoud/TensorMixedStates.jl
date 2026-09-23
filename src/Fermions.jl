export Fermions

"""
    Fermion()

a site type to represent fermion sites (dim is 2)

# Examples

    Fermion()
    Fermion(conserve = N)
    Fermion(conserve = parity(N))

# States

- `"0", "Emp"` : empty state
- `"1", "Occ"` : occupied state

# Operators

- `C` : the destruction operator
- `N` : the number of fermions operator
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

@create_site_module(Fermions, [Fermion, C, N])
