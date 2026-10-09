# The Tj site type, an electron orbital without double occupancy as in the t-J model, with its
# states and operators, gathered in the module Tjs.

export Tjs

"""
    Tj(; conserve = ())

the site type of the t-J model, an electron orbital without double occupancy, of dimension 3,
whose basis is `"Emp"`, `"Up"`, `"Dn"`.

# Examples

    Tj()
    Tj(conserve = (Ntot, 2Sz))

# States

- `"Emp", "0"`         : the empty state
- `"Up", "↑"`          : one electron of spin up
- `"Dn", "↓"`          : one electron of spin down
- `"MixedSpin", "↑|↓"` : one electron of fully mixed spin, the density matrix
  ``(|↑⟩⟨↑| + |↓⟩⟨↓|)/2``

# Operators

- `Cup, Cdn`              : the annihilation operators
- `Nup, Ndn, Ntot`        : the numbers of electrons of spin up, of spin down and in total
- `Sx, Sy, Sz, Sp, Sm`    : the spin operators
- `S2`                    : the total spin squared of the site, `3/4 Ntot`
- `F`                     : the fermionic parity, ``(-1)^{N_{tot}}``
"""
struct Tj <: AbstractSite
    conserve::Conservation
end

Tj(; conserve = ()) = conserving(Tj; conserve)

dim(::Tj) = 3

fermion_parities(::Tj) = (0, 1, 1)

string_state(::Tj, ::String) = nothing

@def_states Tj [
    ["Emp", "0"] => [1., 0., 0.],
    ["Up", "↑"] => [0., 1., 0.],
    ["Dn", "↓"] => [0., 0., 1.],
    ["MixedSpin", "↑|↓"] => [0. 0. 0. ; 0. 0.5 0. ; 0. 0. 0.5],
]

@def_operators Tj [
    Cup = Float64[
        0 1 0
        0 0 0
        0 0 0
    ],
    Cdn = Float64[
        0 0 1
        0 0 0
        0 0 0
    ],
    Sp = Float64[
        0  0  0
        0  0  1
        0  0  0
    ],
    Sm = dag(Sp),
    Nup = dag(Cup) * Cup,
    Ndn = dag(Cdn) * Cdn,
    Ntot = Nup + Ndn,
    Sz = 0.5 * [
        0  0  0
        0  1  0
        0  0 -1
    ],
    Sx = (Sp + Sm) / 2,
    Sy = (Sp - Sm) / (2im),
    S2 = Sx^2 + Sy^2 + Sz^2,
]

@create_site_module(Tjs, [Tj, Cup, Cdn, Nup, Ndn, Ntot, Sx, Sy, Sz, Sp, Sm, S2])
