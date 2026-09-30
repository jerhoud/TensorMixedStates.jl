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
  ``(|↑⟩⟨↑| + |↓⟩⟨↓|)/2``, which a strong conservation of `Ntot` allows where `"FullyMixed"`
  spreads over several numbers of electrons

# Operators

- `Cup, Cdn`              : the destruction operators, fermionic
- `Nup, Ndn, Ntot`        : the numbers of electrons of spin up, of spin down and in total
- `Sx, Sy, Sz, Sp, Sm`    : the spin operators
- `S2`                    : the total spin squared of the site, `3/4 Ntot`, three quarters of
                            the projector on the singly occupied states
- `F`                     : the Jordan-Wigner operator, ``(-1)^{N_{tot}}``
- `Fup, Fdn`              : the partial Jordan-Wigner operators, ``(-1)^{N_\\uparrow}`` and
                            ``(-1)^{N_\\downarrow}``
"""
struct Tj <: AbstractSite
    conserve::String
end

Tj(; conserve = ()) = Tj(conserve_string(Tj(""), conserve))

dim(::Tj) = 3

string_state(::Tj, ::String) = error("no generic state for Tj")

@def_states(Tj(),
[
    ["Emp", "0"] => [1., 0., 0.],
    ["Up", "↑"] => [0., 1., 0.],
    ["Dn", "↓"] => [0., 0., 1.],
    ["MixedSpin", "↑|↓"] => [0. 0. 0. ; 0. 0.5 0. ; 0. 0. 0.5],
])

@def_operators(Tj(),
[
    fermionic_op =>
    [
        Cup = Float64[
            0 1 0
            0 0 0
            0 0 0
        ],
        Cdn = Float64[
            0 0 1
            0 0 0
            0 0 0
        ]
    ],
    involution_op =>
    [
        F = Float64[
            1  0  0
            0 -1  0
            0  0 -1
        ],
        Fup = Float64[
            1  0  0
            0 -1  0
            0  0  1
        ],
        Fdn = Fup * F,
    ],
    plain_op =>
    [
        Sp = Float64[
            0  0  0
            0  0  1
            0  0  0
        ],
        Sm = dag(Sp),
    ],
    selfadjoint_op =>
    [
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
])

@create_site_module(Tjs, [Tj, Cup, Cdn, Fup, Fdn, Nup, Ndn, Ntot, Sx, Sy, Sz, Sp, Sm, S2])
