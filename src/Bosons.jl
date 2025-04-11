export Bosons

struct Boson <: AbstractSite
    dim::Int
end

dim(a::Boson) = a.dim

@def_operators(Boson(2),
[
    F = Id,
    A = s -> [ i==j-1 ? sqrt(i) : 0. for i in 1:dim(s), j in 1:dim(s) ],
    N = s -> [ i==j ? i - 1. : 0. for i in 1:dim(s), j in 1:dim(s) ]
])

module Bosons

import ..Boson, ..N, ..A
export Boson, N, A

end
