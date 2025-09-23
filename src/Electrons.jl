export Electrons

"""
    Eletron()

a site type to represent electron sites (dim is 4)

# Examples

    Electron()

# States

- `"0", "Emp"`   : empty state
- `"Up", "↑"`    : up state
- `"Dn", "↓"`    : down state
- `"UpDn", "↑↓"` : up and down state

# Operators

- `Cup, Cdn`              : the destruction operators
- `Aup, Adn`              : the Jordan-Wigner transforms of Cup and Cdn (...C = FFFA)
- `Nup, Ndn, Nupdn, Ntot` : the numbers operator for up, down, up and down, and total
- `Sx, Sy, Sz, Sp, Sm`    : spin operators
- `Fup, Fdn`              : partial Jordan-Wigner F operators
"""
struct Electron <: AbstractSite end

dim(::Electron) = 4

string_state(::Electron, ::String) = error("no generic state for Electron")

@def_states(Electron(),
[
    ["Emp", "0"] => [1., 0., 0., 0.],
    ["Up", "↑"] => [0., 1., 0., 0.],
    ["Dn", "↓"] => [0., 0., 1., 0.],
    ["UpDn", "↑↓"] => [0., 0., 0., 1.],
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
        Aup = Cup,
        Adn = F * Cdn,
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
        Sy = (Sp - Sm) / (2im)
    ]
])

module Electrons

import ..Electron, ..Cup, ..Cdn, ..Fup, ..Fdn, ..Aup, ..Adn, ..Nup, ..Ndn, ..Nupdn, ..Ntot, ..Sx, ..Sy, ..Sz, ..Sp, ..Sm
export Electron, Cup, Cdn, Fup, Fdn, Aup, Adn, Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz, Sp, Sm

end
