export Electrons

struct Electron <: AbstractSite end

dim(::Electron) = 4

generic_state(::Electron, ::String) = error("no generic state for Electron")

@def_states(Electron(),
[
    ["Emp", "0"] => [1., 0., 0., 0.],
    ["Up", "↑"] => [0., 1., 0., 0.],
    ["Dn", "↓"] => [0., 0., 1., 0.],
    ["UpDn", "↑↓"] => [0., 0., 0., 1.],
])

@def_operators(Electron(),
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
true)

@def_operators(Electron(),
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
    Fdn = Float64[
        1  0  0  0
        0  1  0  0
        0  0 -1  0
        0  0  0 -1
    ],
    Aup = Cup,
    Adn = Float64[
        0 0 1 0
        0 0 0 1
        0 0 0 0
        0 0 0 0
    ],
    Nup = dag(Cup)*Cup,
    Ndn = dag(Cdn)*Cdn,
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
    Sp = Float64[
        0  0  0  0
        0  0  1  0
        0  0  0  0
        0  0  0  0
    ],
    Sm = dag(Sp),
    Sx = (Sp + Sm) / 2,
    Sy = (Sp - Sm) / (2im)
])

module Electrons

import ..Electron, ..Cup, ..Cdn, ..Fup, ..Fdn, ..Aup, ..Adn, ..Nup, ..Ndn, ..Nupdn, ..Ntot, ..Sx, ..Sy, ..Sz, ..Sp, ..Sm
export Electron, Cup, Cdn, Fup, Fdn, Aup, Adn, Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz, Sp, Sm

end
