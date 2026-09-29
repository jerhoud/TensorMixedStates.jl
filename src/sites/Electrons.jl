# The Electron site type, a spinful fermion orbital, with its states and operators, gathered in
# the module Electrons.

export Electrons

"""
    Electron(; conserve = ())

the site type of an electron orbital, of dimension 4.

# Examples

    Electron()
    Electron(conserve = (Ntot, 2Sz))
    Electron(conserve = (strong(Ntot), 2Sz))

# States

- `"Emp", "0"`         : the empty state
- `"Up", "↑"`          : one electron of spin up
- `"Dn", "↓"`          : one electron of spin down
- `"UpDn", "↑↓"`       : two electrons
- `"MixedSpin", "↑|↓"` : one electron of fully mixed spin, the density matrix
  ``(|↑⟩⟨↑| + |↓⟩⟨↓|)/2``, which a strong conservation of `Ntot` allows where `"FullyMixed"`
  spreads over several numbers of electrons

# Operators

- `Cup, Cdn`              : the destruction operators, fermionic
- `Nup, Ndn, Nupdn, Ntot` : the numbers of electrons of spin up, of spin down, of pairs and in
                            total
- `Sx, Sy, Sz, Sp, Sm`    : the spin operators
- `S2`                    : the total spin squared of the site, `3/4 (Ntot - 2 Nupdn)`, three
                            quarters of the projector on the singly occupied states
- `F`                     : the Jordan-Wigner operator, ``(-1)^{N_{tot}}``
- `Fup, Fdn`              : the partial Jordan-Wigner operators, ``(-1)^{N_\\uparrow}`` and
                            ``(-1)^{N_\\downarrow}``
"""
struct Electron <: AbstractSite
    conserve::String
end

Electron(; conserve = ()) = Electron(conserve_string(Electron(""), conserve))

dim(::Electron) = 4

string_state(::Electron, ::String) = error("no generic state for Electron")

@def_states(Electron(),
[
    ["Emp", "0"] => [1., 0., 0., 0.],
    ["Up", "↑"] => [0., 1., 0., 0.],
    ["Dn", "↓"] => [0., 0., 1., 0.],
    ["UpDn", "↑↓"] => [0., 0., 0., 1.],
    ["MixedSpin", "↑|↓"] => [0. 0. 0. 0. ; 0. 0.5 0. 0. ; 0. 0. 0.5 0. ; 0. 0. 0. 0.],
])

@def_operators(Electron(),
[
    fermionic_op =>
    [
        Cup = Float64[
            0 1 0 0
            0 0 0 0
            0 0 0 1
            0 0 0 0
        ],
        Cdn = Float64[
            0 0 1 0
            0 0 0 -1
            0 0 0 0
            0 0 0 0
        ]
    ],
    involution_op =>
    [
        F = Float64[
            1  0  0  0
            0 -1  0  0
            0  0 -1  0
            0  0  0  1
        ],
        Fup = Float64[
            1  0  0  0
            0 -1  0  0
            0  0  1  0
            0  0  0 -1
        ],
        Fdn = Fup * F
    ],
    plain_op =>
    [
        Sp = Float64[
            0  0  0  0
            0  0  1  0
            0  0  0  0
            0  0  0  0
        ],
        Sm = dag(Sp)
    ],
    selfadjoint_op =>
    [
        Nup = dag(Cup) * Cup,
        Ndn = dag(Cdn) * Cdn,
        Nupdn = Float64[
            0 0 0 0
            0 0 0 0
            0 0 0 0
            0 0 0 1
        ],
        Ntot = Nup + Ndn,
        Sz = 0.5 * [
            0  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0  0
        ],
        Sx = (Sp + Sm) / 2,
        Sy = (Sp - Sm) / (2im),
        S2 = Sx^2 + Sy^2 + Sz^2,
    ]
])

@create_site_module(Electrons, [Electron, Cup, Cdn, Fup, Fdn, Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz, Sp, Sm, S2])
