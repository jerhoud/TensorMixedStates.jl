export Fermions

struct Fermion <: AbstractSite end

dim(::Fermion) = 2

@def_states(Fermion(),
[
    "Emp" => [1., 0.],
    "Occ" => [0., 1.],
])


@def_operators(Fermion(),
[
    C = [0. 0. ; 1. 0.],
],
true)

@def_operators(Fermion(),
[
    F = [1. 0. ; 0. -1.],
    A = [0. 0. ; 1. 0.],
    N = [0. 0. ; 0. 1.],
])

module Fermions

import ..Fermion, ..C, ..N, ..A
export Fermion, C, N, A

end
