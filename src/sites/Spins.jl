# The Spin site type, a spin of any integer or half integer value, with its states and
# operators, gathered in the module Spins.

export Spins

"""
    Spin(s; conserve = ())

the site type of a spin `s`, integer or half integer, of dimension ``2s + 1``, its basis going
from ``S_z = s`` down to ``S_z = -s``.

# Examples

    Spin(3/2)
    Spin(2)
    Spin(1, conserve = Sz)       # integer eigenvalues
    Spin(1/2, conserve = 2Sz)    # half integer ones are written doubled
    Spin(1/2, conserve = N)      # or counted from the top, which is integer for any spin

# States

For `m` from `s` down to `-s`, written as an integer or a fraction, `"1"`, `"-1/2"`:

- `"m"`, `"Zm"` : the eigenstate of `Sz` of eigenvalue `m`
- `"Xm"`        : the eigenstate of `Sx` of eigenvalue `m`
- `"Ym"`        : the eigenstate of `Sy` of eigenvalue `m`

A name gives an eigenvalue, where a number gives a basis state counted from 0, as on any site:
on `Spin(1)`, `"1"` is ``S_z = 1`` and `1` is ``S_z = 0``, so that `Proj("1")` and `Proj(1)`
differ.

# Operators

- `Sp, Sm`           : the ``S^+`` and ``S^-`` operators
- `Sx, Sy, Sz, S2`   : the ``S_x``, ``S_y``, ``S_z`` operators and ``S^2``
- `N`                : the number of excitations above the state of maximal ``S_z``,
                       ``s - S_z``, as in the Holstein-Primakoff mapping. Its eigenvalues are
                       integers for every spin, so `Spin(3/2, conserve = N)` conserves what
                       `Spin(3/2, conserve = 2Sz)` does without the doubling. As for a qubit,
                       `Sm` raises it
"""
struct Spin <: AbstractSite
    s::Float64
    conserve::String
    Spin(s::Number, conserve::AbstractString) =
        if isinteger(2 * s)
            new(s, conserve)
        else
            error("Spin requires an half integer as argument")
        end
end

Spin(s::Number; conserve = ()) = Spin(s, conserve_string(Spin(s, ""), conserve))

dim(a::Spin) = Int(2 * a.s + 1)

function string_state(a::Spin, st::String)
    c = st[1]
    if c == 'X' || c == 'Y' || c == 'Z'
        st = st[2:end]
    else
        c = 'Z'
    end
    i = 1 + Int(a.s - parse(Rational{Int}, st))
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
        N = s -> [ i==j ? i - 1. : 0. for i in 1:dim(s), j in 1:dim(s) ],
        Sx = (Sp + Sm) / 2,
        Sy = (Sp - Sm) / (2im),
        S2 = s -> s.s * (s.s + 1) * Id,
    ],
])

@create_site_module(Spins, [Spin, Sp, Sm, Sx, Sy, Sz, S2, N])
