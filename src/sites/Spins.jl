export Spins

"""
    Spin(spin)

A site type for representing spin sites (dim is `2 spin + 1`)

# Example

    Spin(3/2)
    Spin(2)
    Spin(1, conserve = Sz)       # integer eigenvalues
    Spin(1/2, conserve = 2Sz)    # half integer ones are written doubled
    Spin(1/2, conserve = N)      # or counted from the top, which is integer for any spin

# States

"0", "1", "-1"... for integer spins
"1/2", "-1/2", "3/2", "-3/2"... for half integer spins

"X0", "X1/2", "X-1/2", ... for eigenstate of `Sx`
"Y0", "Y1/2", "Y-1/2", ... for eigenstate of `Sy`
"Z0", "Z1/2", "Z-1/2", ... for eigenstate of `Sz` (same as "0", "1/2" ...)

# Operators

- `Sp, Sm`           : the ``S^+`` and ``S^-`` operators
- `Sx, Sy, Sz, S2`   : the ``S_x``, ``S_y``, ``S_z`` operators and ``S^2``
- `N`                : the number of excitations above the state of maximal ``S_z``, that is
                       ``s - S_z``, which is the counting of the Holstein-Primakoff mapping.
                       Its eigenvalues are integers for every spin, half integer ones
                       included, so `Spin(3/2, conserve = N)` conserves the same quantity as
                       `Spin(3/2, conserve = 2Sz)` without the doubling. As for a qubit, it is
                       `Sm` that raises it
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
