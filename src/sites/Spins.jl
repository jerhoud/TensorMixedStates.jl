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

"""
    aklt_state(system; left = "Up", right = "Up")

the AKLT state of a chain of spins one, `Spin(1)`, the exact ground state of
``\\sum_i \\vec S_i \\cdot \\vec S_{i+1} + \\frac{1}{3} (\\vec S_i \\cdot \\vec S_{i+1})^2``,
of energy ``-\\frac{2}{3}(n - 1)``: the MPS of bond dimension 2 of tensors
``A^{+1} = \\sqrt{2/3}\\, \\sigma^+``, ``A^0 = -\\sqrt{1/3}\\, \\sigma^z`` and
``A^{-1} = -\\sqrt{2/3}\\, \\sigma^-``. On an open chain the ground state is fourfold degenerate,
the virtual spins one half at the two ends being free: `left` and `right` set them, `"Up"` or
`"Dn"`, or a vector of two components, and the total ``S^z`` is that of `left` less that of
`right`, 0 by default.

# Examples

    aklt_state(System(20, Spin(1)))
    aklt_state(System(20, Spin(1, conserve = Sz)); right = "Dn")    # total Sz of 1
"""
function aklt_state(system::System; left = "Up", right = "Up")
    if !all(s -> s isa Spin && s.s == 1, system.sites)
        error("aklt_state takes a chain of spins one, Spin(1)")
    end
    edge(v) = v == "Up" ? [1., 0.] : v == "Dn" ? [0., 1.] : v
    l, r = edge(left), edge(right)
    # in the basis of Spin(1), m = 1, 0, -1, and of the virtual spins, up then down
    A = zeros(2, 3, 2)
    A[1, 1, 2] = sqrt(2 / 3)
    A[1, 2, 1], A[2, 2, 2] = -sqrt(1 / 3), sqrt(1 / 3)
    A[2, 3, 1] = -sqrt(2 / 3)
    n = length(system)
    tensors = map(1:n) do k
        t = A
        if k == 1
            t = reshape(sum(l[a] * t[a, :, :] for a in 1:2), 1, 3, 2)
        end
        if k == n
            t = reshape(sum(t[:, :, b] * r[b] for b in 1:2), size(t, 1), 3, 1)
        end
        t
    end
    return mps_state(system, tensors)
end

@create_site_module(Spins, [aklt_state, Spin, Sp, Sm, Sx, Sy, Sz, S2, N])
