# The Spin site type, a spin of any integer or half integer value, with its states and
# operators, and the AKLT state, gathered in the module Spins.

export Spins

"""
    Spin(s; conserve = ())

the site type of a spin `s`, integer or half integer, of dimension ``2s + 1``, its basis going
from ``S_z = s`` down to ``S_z = -s``.

# Examples

    Spin(3/2)
    Spin(1, conserve = Sz)       # integer eigenvalues
    Spin(1/2, conserve = 2Sz)    # half integer ones are written doubled
    Spin(1/2, conserve = N)      # or counted from the top, which is integer for any spin

# States

For `m` from `s` down to `-s`, written as an integer or a fraction, `"1"`, `"-1/2"`:

- `"m"`, `"Zm"` : the eigenstate of `Sz` of eigenvalue `m`
- `"Xm"`        : the eigenstate of `Sx` of eigenvalue `m`
- `"Ym"`        : the eigenstate of `Sy` of eigenvalue `m`

and `"Up"`, `"↑"`, `"Dn"`, `"↓"` for ``m = s`` and ``m = -s``.

A name gives an eigenvalue: on `Spin(1)`, `"1"` is ``S_z = 1``, and `Level => 1` the second
basis state, ``S_z = 0``.

# Operators

- `Sp, Sm`           : the ``S^+`` and ``S^-`` operators
- `Sx, Sy, Sz, S2`   : the ``S_x``, ``S_y``, ``S_z`` operators and ``S^2``
- `N`                : the number of excitations ``s - S_z``, integer for every spin
"""
struct Spin <: AbstractSite
    s::Float64
    conserve::Conservation
    Spin(s::Number, conserve::Conservation) =
        if isinteger(2 * s) && s ≥ 0
            new(s, conserve)
        else
            error("Spin requires a half integer of 0 or more as argument, not $s")
        end
end

Spin(s::Number; conserve = ()) = conserving(Spin, s; conserve)

dim(a::Spin) = Int(2 * a.s + 1)

@def_operators Spin [
    Sp = s -> [ i == j - 1 ? sqrt(s.s * (s.s + 1) - (s.s - i + 1) * (s.s - i)) : 0.
                for i in 1:dim(s), j in 1:dim(s) ],
    Sz = s -> [ i == j ? s.s - i + 1 : 0. for i in 1:dim(s), j in 1:dim(s) ],
    N = s -> [ i == j ? i - 1. : 0. for i in 1:dim(s), j in 1:dim(s) ],
    S2 = s -> s.s * (s.s + 1) * Id,
]

@def_states Spin [
    ["Up", "↑"] => s -> [ i == 1 ? 1. : 0. for i in 1:dim(s) ],
    ["Dn", "↓"] => s -> [ i == dim(s) ? 1. : 0. for i in 1:dim(s) ],
]

"""
    rational_value(s)

the number the string `s` writes as an integer or a fraction, `"1"`, `"-3/2"`, or `nothing`
"""
function rational_value(s::AbstractString)
    parts = split(s, '/')
    if length(parts) > 2
        return nothing
    end
    n = tryparse(Int, parts[1])
    d = length(parts) == 2 ? tryparse(Int, parts[2]) : 1
    return isnothing(n) || isnothing(d) || d == 0 ? nothing : n // d
end

# the eigenstates of Sx and Sy are rotations of those of Sz, exp(-iπ/2 Sy) and exp(iπ/2 Sx)
function string_state(a::Spin, st::String)
    axis, value = !isempty(st) && st[1] in ('X', 'Y', 'Z') ? (st[1], st[nextind(st, 1):end]) :
                                                              ('Z', st)
    m = rational_value(value)
    if isnothing(m) || !isinteger(a.s - m) || abs(m) > a.s
        return nothing
    end
    v = zeros(Float64, dim(a))
    v[1 + Int(a.s - m)] = 1.0
    if axis == 'Z'
        return v
    elseif axis == 'X'
        return exp(-0.5 * im * pi * matrix(Sy, a)) * v
    else
        return exp(0.5 * im * pi * matrix(Sx, a)) * v
    end
end

"""
    aklt_state(system; left = "Up", right = "Up")

the AKLT state of a chain of spins one, `Spin(1)`, the exact ground state of
``\\sum_i \\vec S_i \\cdot \\vec S_{i+1} + \\frac{1}{3} (\\vec S_i \\cdot \\vec S_{i+1})^2``,
as an MPS of bond dimension 2. `left` and `right` set the free virtual spins one half at the two
ends, pure states as `Spin(1/2)` takes them: `"Up"`, `"Dn"`, `"X1/2"` or a vector of two
components.

# Examples

    aklt_state(System(20, Spin(1)))
    aklt_state(System(20, Spin(1, conserve = Sz)); right = "Dn")    # total Sz of 1
"""
function aklt_state(system::System; left = "Up", right = "Up")
    if !all(s -> s isa Spin && s.s == 1, system.sites)
        error("aklt_state takes a chain of spins one, Spin(1)")
    end
    function edge(v)
        e = local_state(Spin(1/2), v)
        if !(e isa AbstractVector)
            error("aklt_state takes pure states for its free spins one half, not $(repr(v))")
        end
        return e
    end
    l, r = edge(left), edge(right)
    # in the basis of Spin(1), m = 1, 0, -1, and of the virtual spins, up then down, of the
    # element type of the edges for every tensor, which mps_state takes
    A = zeros(promote_type(Float64, eltype(l), eltype(r)), 2, 3, 2)
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

@create_site_module(Spins, [Spin, aklt_state, Sp, Sm, Sx, Sy, Sz, S2, N])
