# The Electron site type, a spinful fermion orbital, with its states and operators, gathered in
# the module Electrons.

export Electrons

"""
    Electron(; conserve = ())

the site type of an electron orbital, of dimension 4, whose basis is `"Emp"`, `"Up"`, `"Dn"`,
`"UpDn"`, the doubly occupied state being ``c^\\dagger_\\uparrow c^\\dagger_\\downarrow |0⟩``.

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

- `Cup, Cdn`              : the annihilation operators
- `Nup, Ndn, Nupdn, Ntot` : the numbers of electrons of spin up, of spin down, of pairs and in
                            total
- `Sx, Sy, Sz, Sp, Sm`    : the spin operators
- `S2`                    : the total spin squared of the site, `3/4 (Ntot - 2 Nupdn)`
- `F`                     : the fermionic parity, ``(-1)^{N_{tot}}``
- `Fup, Fdn`              : the partial parities, ``(-1)^{N_\\uparrow}`` and
                            ``(-1)^{N_\\downarrow}``
"""
struct Electron <: AbstractSite
    conserve::Conservation
end

Electron(; conserve = ()) = conserving(Electron; conserve)

dim(::Electron) = 4

fermion_parities(::Electron) = (0, 1, 1, 0)

string_state(::Electron, ::String) = nothing

@def_states Electron [
    ["Emp", "0"] => [1., 0., 0., 0.],
    ["Up", "↑"] => [0., 1., 0., 0.],
    ["Dn", "↓"] => [0., 0., 1., 0.],
    ["UpDn", "↑↓"] => [0., 0., 0., 1.],
    ["MixedSpin", "↑|↓"] => [0. 0. 0. 0. ; 0. 0.5 0. 0. ; 0. 0. 0.5 0. ; 0. 0. 0. 0.],
]

@def_operators Electron [
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
    ],
    Fup = Float64[
        1  0  0  0
        0 -1  0  0
        0  0  1  0
        0  0  0 -1
    ],
    Fdn = Fup * F,
    Sp = Float64[
        0  0  0  0
        0  0  1  0
        0  0  0  0
        0  0  0  0
    ],
    Sm = dag(Sp),
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

fermion_species(::Electron) = (Cup, Cdn)

@create_site_module(Electrons, [Electron, Cup, Cdn, Fup, Fdn, Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz,
                                Sp, Sm, S2])
