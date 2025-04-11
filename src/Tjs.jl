export Tjs

struct Tj <: AbstractSite end

dim(::Tj) = 3

generic_state(::Tj, ::String) = error("no generic state for Tj")

@def_states(Tj(),
[
    ["Emp", "0"] => [1., 0., 0.],
    ["Up", "↑"] => [0., 1., 0.],
    ["Dn", "↓"] => [0., 0., 1.],
])

@def_operators(Tj(),
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
    ],
],
true)

@def_operators(Tj(),
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
    Fdn = Float64[
        1  0  0
        0  1  0
        0  0 -1
    ],
    Aup = Cup,
    Adn = Cdn,
    Nup = dag(Cup)*Cup,
    Ndn = dag(Cdn)*Cdn,
    Ntot = Nup + Ndn,
    Sz = 0.5 * [
        0  0  0
        0  1  0
        0  0 -1
    ],
    Sp = Float64[
        0  0  0
        0  0  1
        0  0  0
    ],
    Sm = dag(Sp),
    Sx = (Sp + Sm) / 2,
    Sy = (Sp - Sm) / (2im)
])

module Tjs

import ..Tj, ..Cup, ..Cdn, ..Fup, ..Fdn, ..Aup, ..Adn, ..Nup, ..Ndn, ..Nupdn, ..Ntot, ..Sx, ..Sy, ..Sz, ..Sp, ..Sm
export Tj, Cup, Cdn, Fup, Fdn, Aup, Adn, Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz, Sp, Sm

end
