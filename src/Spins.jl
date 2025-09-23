export Spins

"""
    Spin(spin)

A site type for representing spin sites (dim is `2 spin + 1`)

# Example

    Spin(3/2)
    Spin(2)

# States

"0", "1", "-1"... for integer spins
"1/2", "-1/2", "3/2", "-3/2"... for half integer spins

"X0", "X1/2", "X-1/2", ... for eigenstate of `Sx`
"Y0", "Y1/2", "Y-1/2", ... for eigenstate of `Sy`
"Z0", "Z1/2", "Z-1/2", ... for eigenstate of `Sz` (same as "0", "1/2" ...)

# Operators

- `Sp, Sm`           : the ``S^+`` and ``S^-`` operators
- `Sx, Sy, Sz, S2`   : the ``S_x``, ``S_y``, ``S_z`` operators and ``S^2``
"""
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

function string_state(a::Spin, st::String)
    c = st[1]
    if c == 'X' || c == 'Y' || c == 'Z'
        st = st[2:end]
    else
        c = 'Z'
    end
    i = 1 + Int(a.s - eval(Meta.parse(st)))
    v = zeros(Float64, dim(a))
    v[i] = 1.0
    if c == 'Z'
        return v
    elseif c == 'X'
        return exp(-0.5 * im * pi * matrix(Sy, a)) * v
    else
        return exp(0.5 * im * pi * matrix(Sx, a)) * v
    end
end

@def_operators(Spin(0),
[
    plain_op =>
    [
        Sp = s -> [ i==j-1 ? sqrt(s.s*(s.s + 1) - (s.s - i + 1)*(s.s - i)) : 0. for i in 1:dim(s), j in 1:dim(s) ],
        Sm = dag(Sp),
    ],
    selfadjoint_op =>
    [
        Sz = s -> [ i==j ? s.s - i + 1 : 0. for i in 1:dim(s), j in 1:dim(s) ],
        Sx = (Sp + Sm) / 2,
        Sy = (Sp - Sm) / (2im),
        S2 = s -> s.s * (s.s + 1) * Id,
    ],
])

module Spins

import ..Spin, ..Sp, ..Sm, ..Sx, ..Sy, ..Sz, ..S2
export Spin, Sp, Sm, Sx, Sy, Sz, S2

end
