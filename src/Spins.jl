export Spins

struct Spin <: AbstractSite
    s::Float64
    Spin(s::Number) =
        if isinteger(2 * s)
            new(s)
        else
            error("Spin requires an half integer as argument")
        end
end

dim(a::Spin) = Int(2 * a.s + 1)
function generic_state(::AbstractSite, st::String)
    i = 1 + 2 * Int(parse(Float64, st) + a.s)
    v = zeros(Float64, dim(site))
    v[i] = 1.0
    return v
end

@def_operators(Spin(0),
[
    F = Id,
    Sz = s -> [ i==j ? s.s - i + 1 : 0. for i in 1:dim(s), j in 1:dim(s) ],
    Sp = s -> [ i==j-1 ? sqrt(s.s*(s.s+1) - (s.s - i + 1)*(s.s - i)) : 0. for i in 1:dim(s), j in 1:dim(s) ],
    Sm = dag(Sp),
    Sx = (Sp + Sm) / 2,
    Sy = (Sp - Sm) / (2im),
    S2 = s -> s.s * (s.s + 1) * identity_operator(s),
])

module Spins

import ..Spin, ..Sp, ..Sm, ..Sx, ..Sy, ..Sz, ..S2
export Spin, Sp, Sm, Sx, Sy, Sz, S2

end
