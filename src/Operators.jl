# The operator algebra: the types of operators, generic or placed on sites, acting on pure
# states or on density matrices, and how they are built, combined, compared and printed, before
# any site or state is involved.

export Representation, Pure, Mixed, GenericOp, IndexedOp, SimpleOp
export OpType, plain_op, fermionic_op, selfadjoint_op, involution_op
export Op, Operator, Id, F, Proj, Gate, Dissipator, Evolver, Left, Right, SetState
export named, parity
export dag, ⊗, isfermionic, hasfermionic, lindblad_terms, map_sites

############# Types ################

"""
    abstract type Representation

the supertype of the representations in which `CreateState` creates a state: `Pure` and
`Mixed`, and those an extension defines, whose states are `AbstractState`s.
"""
abstract type Representation end

"""
    abstract type PM <: Representation

the supertype of `Pure` and `Mixed`, the two representations of a state that the package
defines, which parametrize operators as well as states
"""
abstract type PM <: Representation end

"""
    struct Pure <: PM
    Pure()

the pure representation, in which a state is a wave function. It parametrizes states and
operators, `State{Pure}` and `Op{Pure}`, and `Pure()` selects it where a representation is
passed as a value.

# Examples

    State{Pure}(System(4, Qubit()), "Up")
    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up")
"""
struct Pure <: PM end


"""
    struct Mixed <: PM
    Mixed()

the mixed representation, in which a state is a density matrix. It parametrizes states and
operators, `State{Mixed}` and `Op{Mixed}`, and `Mixed()` selects it where a representation is
passed as a value. Superoperators such as `Gate`, `Dissipator` or `Left` act on it only.

# Examples

    State{Mixed}(System(3, Qubit()), ["FullyMixed", "Up", "Dn"])
    CreateState(type = Mixed(), system = System(3, Qubit()), state = ["Up", "Dn", "Up"])
"""
struct Mixed <: PM end

"""
    abstract type GI

the supertype of `Generic` and `Indexed`
"""
abstract type GI end

"""
    Generic

the parameter of `Op` for a generic operator, not yet placed on sites, as `X` or `X ⊗ Y`
"""
struct Generic <: GI end

"""
    Indexed

the parameter of `Op` for an operator placed on sites, as `X(1)` or `X(1) * Y(3)`
"""
struct Indexed <: GI end

"""
    Op{R <: PM, T <: GI, N}

the supertype of all operators.

# Type parameters

- `R`: `Pure` or `Mixed`, the representation of the states the operator acts on
- `T`: `Generic` for an operator not yet placed, `Indexed` for one placed on sites
- `N`: the number of sites a generic operator acts on, 1 for a placed one
"""
abstract type Op{R <: PM, T <: GI, N} end

"""
    GenericOp{R, N}

the generic operators, not yet placed on sites, that is `Op{R, Generic, N}`.
"""
const GenericOp{R, N} = Op{R, Generic, N}

"""
    IndexedOp{R}

the operators placed on sites, that is `Op{R, Indexed, 1}`.
"""
const IndexedOp{R} = Op{R, Indexed, 1}

"""
    SimpleOp

the generic operators of one site on pure states, that is `GenericOp{Pure, 1}`. It is an
abstract type: `X` is a `SimpleOp`, and so is any combination of one site such as `2X`,
`X * Y` or `exp(X)`.
"""
const SimpleOp = GenericOp{Pure, 1}

############## Showing ###############

"""
    show_func(io, name, args, kwargs = (;))

print `name(args; kwargs)`, `args` being a vector of arguments or a single one
"""
function show_func(io::IO, name, args, kwargs=(;))
    print(io, name, "(")
    if args isa Vector
        join(io, (repr(a) for a in args), ",")
    else
        show(io, args)
    end
    if !isempty(kwargs)
        print(io, ";")
        join(io, ("$sym=$(repr(val))" for (sym, val) in pairs(kwargs)), ",")
    end
    print(io, ")")
end

"""
    print_number(io, x)

print the real number `x`, a float of integer value as that integer, `2` rather than `2.0`: an
integer coefficient is stored as a float, see `float_coef`, and prints as it was written.
Beyond `maxintfloat`, where a float no longer stands for a single integer, it prints as a float.
"""
print_number(io::IO, x::AbstractFloat) =
    if isinteger(x) && abs(x) < maxintfloat(Float64)
        print(io, Int(x))
    else
        print(io, x)
    end
print_number(io::IO, x::Real) = print(io, x)

"""
    print_coef(io, a)

print the number `a` as the coefficient in front of an operator: nothing for 1, `-` for -1,
`2im*` for an imaginary number, a rational or a complex number in parentheses, the floats of
integer value as integers, see `print_number`
"""
print_coef(io::IO, a::Number) =
if a ≠ 1
    if a == -1
        print(io, "-")
    elseif isa(a, Complex)
        if imag(a) == 0
            print_coef(io, real(a))
        elseif real(a) == 0
            print_coef(io, imag(a))
            print(io, "im*")
        elseif isfinite(a)
            print(io, "(")
            print_number(io, real(a))
            print(io, imag(a) < 0 ? " - " : " + ")
            print_number(io, abs(imag(a)))
            print(io, "im)")
        else
            print(io, "(", a, ")")
        end
    elseif isa(a, Rational)
        # (1//2)X and not 1//2X, which reads 1//(2X)
        print(io, "(", a, ")")
    else
        print_number(io, a)
    end
end

"""
    infix(io, subs, op)

print `subs` separated by `op`, the operands after the first at a precedence one higher: `*`
and `⊗` are left associative, and `X ⊗ (Y*Z)` must not print as `X⊗Y*Z`, which reads
`(X⊗Y)*Z`.
"""
function infix(io::IO, subs, op::String)
    p = get(io, :precedence, 0)
    for (k, s) in enumerate(subs)
        if k > 1
            print(io, op)
        end
        print(k == 1 ? io : IOContext(io, :precedence => p + 1), s)
    end
end

"""
    paren(f, io, out_prec, in_prec = out_prec)

print through `f(io)` an expression of precedence `out_prec`, in parentheses when the
surrounding precedence is higher, its inside printed at precedence `in_prec`
"""
function paren(f, io::IO, out_prec::Int, in_prec::Int = out_prec)
    ext_prec = get(io, :precedence, 0)
    if out_prec < ext_prec
        print(io, "(")
    end
    f(IOContext(io, :precedence => in_prec))
    if out_prec < ext_prec
        print(io, ")")
    end
end


"""
    no_signed_zero(x)

`x` with its signed zeros made positive (`-0.0 + 1.0im` becomes `0.0 + 1.0im`), any other
value unchanged. Operators store their numbers through it: `-0.0 == 0.0`, but `isless` and
`hash` tell them apart, so `simplify` could leave two equal terms unmerged and `measure` could
compute the same operator twice.
"""
no_signed_zero(x::Union{AbstractFloat, Complex{<:AbstractFloat}}) = x + zero(x)
no_signed_zero(x::AbstractArray) = map(no_signed_zero, x)
no_signed_zero(x) = x

"""
    state_key(state)

the key a `Proj` or a `SetState` is ordered by: the kind of its state (an index, a name or an
array), then its value, an array by its size and by the real and imaginary parts of its
elements. The order is total, complex numbers included, and two keys tie exactly when the
states are `==`, so that `simplify` merges the projectors on `[1, 0]` and on `[1.0, 0.0]`.
"""
state_key(x::Int) = (1, x)
state_key(x::AbstractString) = (2, x)
state_key(x::AbstractArray) = (3, size(x), [ (real(y), imag(y)) for y in vec(x) ])

############### Operator ###############

"""
    @enum OpType

the possible types of an `Operator`.

# Enumeration values

- `plain_op`: no particular property
- `fermionic_op`: fermionic, the Jordan-Wigner transform applies to it
- `selfadjoint_op`: invariant under `dag`
- `involution_op`: invariant under `dag`, and its square is the identity

`simplify` relies on the type: `dag(A)` is `A` for a self-adjoint `A`, `A * A` is `Id` for an
involution, and the `F` of a site crosses a fermionic operator with a sign, any other without.
The type is therefore checked against the matrix whenever the operator is placed on a site,
and a wrong type is refused.
"""
@enum OpType plain_op fermionic_op selfadjoint_op involution_op

"""
    plain_op

the `OpType` of an operator with no particular properties.
"""
plain_op

"""
    fermionic_op

the `OpType` of a fermionic operator, for which the Jordan-Wigner transform must be used.
"""
fermionic_op

"""
    selfadjoint_op

the `OpType` of an operator invariant under `dag`.
"""
selfadjoint_op

"""
    involution_op

the `OpType` of an operator invariant under `dag` and whose square is the identity.
"""
involution_op

"""
    struct Operator{N} <: GenericOp{Pure, N}
    Operator{N}(name, expr, type)

a named operator of `N` sites, as `X`, `Swap` or `C`, of the given `OpType`. `expr` defines
it: a matrix, a function of the sites, an expression of other operators, or `nothing` for an
operator whose matrix each site type gives through `@def_operators`. An operator of several
sites cannot be `fermionic_op`. To define an operator of one's own, see `named`.

# Examples
    Operator{1}("X", nothing, involution_op)  # value predefined by the sites
    Operator{1}("Z", [1 0 ; 0 -1], involution_op)
    Operator{2}("Swap", [ 1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1], involution_op)
    Operator{1}("Sx", (Sp + Sm) / 2, selfadjoint_op)
    Operator{1}("C", [0 1 ; 0 0], fermionic_op)
"""
struct Operator{N} <: GenericOp{Pure, N}
    name::String
    expr::Union{Nothing, Matrix, Function, GenericOp{Pure, N}}
    type::OpType
    Operator{N}(name::String, expr::Union{Nothing, Matrix, Function, GenericOp{Pure, N}}, type::OpType) where N =
        if N > 1 && type == fermionic_op
            error("cannot deal with a fermionic multi site operator $name")
        else
            new{N}(name, no_signed_zero(expr), type)
        end
end

show(io::IO, op::Operator) =
    print(io, op.name)

isless(a::Operator, b::Operator) = isless(a.name, b.name)


############ Identity ############

"""
    struct IdentityOp{R, T, N}

the identity, one value for each kind of operator (pure or mixed, generic on `N` sites or
placed). Every way of writing an identity gives it, `Id ⊗ Id`, `Left(Id)`, `Right(Id)`,
`Gate(Id)` and `Id(3)` included, so an identity is recognized by its type alone. Placed, it is
the identity of the whole system and has no site: the code that needs its tensor on a site
builds it there.
"""
struct IdentityOp{R, T, N} <: Op{R, T, N} end

IdentityOp(::Op{R, T, N}) where {R, T, N} = IdentityOp{R, T, N}()

"""
    Id

the identity operator, defined on every site type.
"""
const Id = IdentityOp{Pure, Generic, 1}()

show(io::IO, ::IdentityOp{Pure, Generic, N}) where N = print(io, join(fill("Id", N), "⊗"))
show(io::IO, ::IdentityOp{Mixed, Generic, N}) where N = print(io, "Gate(", join(fill("Id", N), "⊗"), ")")
show(io::IO, ::IdentityOp{Pure, Indexed}) = print(io, "Id")
show(io::IO, ::IdentityOp{Mixed, Indexed}) = print(io, "Gate(Id)")

isless(::IdentityOp, ::IdentityOp) = false


################ Sums ##############

"""
    struct SumOp{R, T, N} <: Op{R, T, N}

a sum of operators, nested sums flattened and terms of coefficient zero left out. An empty
sum is `0 * Id`, a sum of one term that term.
"""
struct SumOp{R, T, N} <: Op{R, T, N}
    subs::Vector{<:Op{R, T, N}}
    function SumOp(subs::Vector{<:Op{R, T, N}}) where {R, T, N}
        # gathered in one pass: a sum built term by term, as `sum` does from a generator, copies
        # what it holds at every term, and concatenating and filtering it again made that cost
        # grow much faster, 3.8 s for 20000 terms
        s = Op{R, T, N}[]
        for x in subs
            if x isa SumOp
                # its terms are already those of coefficient other than zero
                append!(s, x.subs)
            elseif scalarcoef(x) ≠ 0
                # a term of coefficient zero is left out: the zero operator has every parity,
                # and kept, it made `0C + dag(C)` a sum of fermionic and non fermionic operators
                push!(s, x)
            end
        end
        if isempty(s)
            return 0 * IdentityOp{R, T, N}()
        elseif length(s) == 1
            return s[1]
        else
            return new{R, T, N}(s)
        end
    end
end

"""
    sumsubs(a)

the terms of a sum, an operator that is not a sum being its only term
"""
sumsubs(a::SumOp) = a.subs
sumsubs(a::Op) = [a]

(a::Op{R, T, N} + b::Op{R, T, N}) where {R, T, N} = SumOp([a, b])
(a::Op - b::Op) = a + (-b)

show(io::IO, a::SumOp) =
    paren(io, Base.operator_precedence(:+)) do io
        n = length(a.subs)
        subs = if n > 6 && get(io, :compact, false)
            dots = startswith(sprint(print, a.subs[4]; context = io), '-') ? "-..." : "..."
            [a.subs[1:3]; dots; a.subs[n-1:n]]
        else
            a.subs
        end
        for (k, x) in enumerate(subs)
            s = sprint(print, x; context = io)
            if k > 1 && s[1] ≠ '-'
                print(io, "+")
            end
            print(io, s)
        end
    end

isless(a::SumOp, b::SumOp) = isless(a.subs, b.subs)


################ Product by a number #############

"""
    float_coef(x)

the number `x` as an operator stores it as its coefficient: an integer of fixed width, or a
complex number of them, as a float, whose products cannot wrap around as those of such integers
silently do, `prod(2Sz(i) for i in 1:63)` having had the coefficient `-2^63`. Any other number
is kept as it is, a rational refusing to overflow. The float prints as the integer, see
`print_number`.
"""
float_coef(x::Union{Bool, Base.BitInteger, Complex{<:Union{Bool, Base.BitInteger}}}) = float(x)
float_coef(x::Number) = x

"""
    struct ScalarOp{R, T, N} <: Op{R, T, N}

a number times an operator, the operator never a `ScalarOp` itself: `2 * (3X)` is `6X`. A
coefficient of 1 gives the operator, of 0 gives `0 * Id`, and a number times a sum multiplies
each term. An integer coefficient is stored as a float, see `float_coef`.
"""
struct ScalarOp{R, T, N} <: Op{R, T, N}
    coef::Number
    arg::Op{R, T, N}
    ScalarOp(coef::Number, arg::Op{R, T, N}) where {R, T, N} =
        if coef == 0
            new{R, T, N}(0., IdentityOp{R, T, N}())
        elseif coef == 1
            arg
        elseif arg isa SumOp
            SumOp(map(x -> coef * x, arg.subs))
        else
            new{R, T, N}(no_signed_zero(float_coef(coef * scalarcoef(arg))), scalararg(arg))
        end
end

"""
    scalarcoef(a)

the coefficient of an operator, 1 unless it is a `ScalarOp`
"""
scalarcoef(a::ScalarOp) = a.coef
scalarcoef(::Op) = 1

"""
    scalararg(a)

the operator without its coefficient, see `scalarcoef`
"""
scalararg(a::ScalarOp) = a.arg
scalararg(a::Op) = a

(a::Number * b::Op{R, T, N}) where {R, T, N} = 
    ScalarOp(a * scalarcoef(b), scalararg(b))
(a::Op * b::Number) = b * a
(a::Op / b::Number) = inv(b) * a
-(a::Op) = -1 * a

show(io::IO, a::ScalarOp) =
    paren(io, Base.operator_precedence(:*)) do io
        print_coef(io, a.coef)
        print(io, a.arg)
    end

isless(a::ScalarOp, b::ScalarOp) =
    isless((a.arg, real(a.coef), imag(a.coef)), (b.arg, real(b.coef), imag(b.coef)))


############### Products ###############

"""
    struct ProdOp{R, T, N} <: Op{R, T, N}

a product of operators, nested products flattened and the coefficients of the factors
gathered in front of it. An empty product is `Id`, a product of one factor that factor.
"""
struct ProdOp{R, T, N} <: Op{R, T, N}
    subs::Vector{<:Op{R, T, N}}
    function ProdOp(subs::Vector{<:Op{R, T, N}}) where {R, T, N}
        c = prod(scalarcoef.(subs))
        if isempty(subs)
            return c * IdentityOp{R, T, N}()
        end
        s = reduce(vcat, prodsubs.(subs))
        if length(s) == 1
            c * s[1]
        else
            c * new{R, T, N}(s)
        end
    end
end

"""
    prodsubs(a)

the factors of a product, without their coefficients, an operator that is not a product
being its only factor
"""
prodsubs(a::ProdOp) = a.subs
prodsubs(a::ScalarOp) = prodsubs(a.arg)
prodsubs(a::Op) = [a]

(a::Op{R, T, N} * b::Op{R, T, N}) where {R, T, N} =
    ProdOp([a, b])

show(io::IO, a::ProdOp) =
    paren(io, Base.operator_precedence(:*)) do io
        n = length(a.subs)
        if n > 6 && get(io, :compact, false)
            infix(io, [a.subs[1:3]; "..."; a.subs[n-1:n]], "*")
        else
            infix(io, a.subs, "*")
        end
    end

isless(a::ProdOp, b::ProdOp) = isless(a.subs, b.subs)


############### Tensor products ############

"""
    struct TensorOp{N} <: GenericOp{Pure, N}

a tensor product of generic operators on pure states, acting on `N` sites, the sum of theirs.
The coefficients of the factors are gathered in front of it, and a product of identities is
the identity of `N` sites.
"""
struct TensorOp{N} <: GenericOp{Pure, N}
    subs::Vector{<:GenericOp{Pure}}
    TensorOp{N}(subs::Vector{<:GenericOp{Pure}}) where N =
        if length(subs) == 1
            subs[1]
        else
            c = prod(scalarcoef.(subs))
            s = scalararg.(subs)
            if all(x -> x isa IdentityOp, s)
                return c * IdentityOp{Pure, Generic, N}()
            end
            c * new{N}(s)
        end
end

"""
    tensorsubs(a)

the factors of a tensor product, an operator that is not one being its only factor
"""
tensorsubs(a::TensorOp) = a.subs
tensorsubs(a::GenericOp) = [a]

"""
    nsites(a)

the number of sites the generic operator `a` acts on
"""
nsites(::GenericOp{R, N}) where {R, N} = N

"""
    factor_sites(a::TensorOp)

the positions, among the sites of the tensor product `a`, that each of its factors acts on, as
ranges in the order of the factors: the first factor takes the first sites, the next one the
following ones, and so on
"""
function factor_sites(a::TensorOp)
    stop = 0
    return map(a.subs) do o
        start = stop + 1
        stop += nsites(o)
        start:stop
    end
end

"""
    op1 ⊗ op2

the tensor product of generic operators on pure states, acting on the sites of `op1` followed
by those of `op2`: `(A ⊗ B)(i, j)` is `A(i) * B(j)`. `⊗` is typed `\\otimes`, and
`tensor(op1, op2, ...)` is the same product.

# Examples

    mycontrolled(op) = Proj("Up") ⊗ Id + Proj("Dn") ⊗ op     # op of one site
    Rxy(t) = exp(-im * t * (X ⊗ X + Y ⊗ Y) / 4)

"""
(a::GenericOp{Pure, N} ⊗ b::GenericOp{Pure, M}) where {N, M} =
    TensorOp{N + M}([tensorsubs(a) ; tensorsubs(b)])
tensor(a::GenericOp, b::GenericOp) = a ⊗ b
tensor(a::GenericOp, b::GenericOp, c::GenericOp, d::GenericOp...) = tensor(a ⊗ b, c, d...)

show(io::IO, a::TensorOp) =
    paren(io, Base.operator_precedence(:⊗)) do io
        infix(io, a.subs, "⊗")
    end

isless(a::TensorOp, b::TensorOp) = isless(a.subs, b.subs)


############# Jordan_Wigner transformation ##############

"""
    struct JW <: SimpleOp

a fermionic operator of one site once `simplify` has put its Jordan-Wigner string in front
of it: `C(5)` becomes `Multi_F{Pure}(1, 4, false, false) * JW(C)(5)`. It has the matrix of
the operator and anticommutes with `F`, but it is not fermionic, its string being in place.
"""
struct JW <: SimpleOp
    arg::Operator{1}
end

isless(a::JW, b::JW) = isless(a.arg, b.arg)

"""
    struct JW_F

the type of `F`
"""
struct JW_F <: SimpleOp end

"""
    F

the Jordan-Wigner factor of a site, defined on every site type: the identity on a site that
declares none, which is not fermionic. It has to be an involution.
"""
const F = JW_F()

show(io::IO, ::JW_F) =
    print(io, "F")

isless(::JW_F, ::JW_F) = false

"""
    Multi_F{R}(start, stop, left, right)

a Jordan-Wigner string, `F` on each site from `start` to `stop`: `C(5)` becomes
`Multi_F{Pure}(1, 4, false, false) * JW(C)(5)`. On mixed states, `left` and `right` tell
whether it is `Left(F)`, `Right(F)` or both on each site. An empty string, or a mixed one on
neither side, is the identity, and a string of one site is the `F` of that site.
"""
struct Multi_F{R} <: IndexedOp{R}
    start::Int
    stop::Int
    left::Bool
    right::Bool
    Multi_F{R}(start::Int, stop::Int, left::Bool, right::Bool) where R =
        if start > stop || (R == Mixed && !left && !right)
            IdentityOp{R, Indexed, 1}()
        elseif start < stop
            new{R}(start, stop, left, right)
        elseif R == Pure
            F(start)
        elseif !left
            Right(F)(start)
        elseif right
            (Left(F)*Right(F))(start)
        else
            Left(F)(start)
        end
end

isless(a::Multi_F, b::Multi_F) =
    isless((a.start, a.stop, a.left, a.right), (b.start, b.stop, b.left, b.right))

############ Proj ################

"""
    Proj(state)

the projector ``|s\\rangle\\langle s|`` on a pure state of one site, given by its name, by
its vector in the basis of the site, taken as it is without normalization, or by the number
of a basis state, counted from 0.

# Examples

    Proj("Up")
    Proj([1, 0])
    Proj(1)       # the second basis state
"""
struct Proj <: SimpleOp
    state::Union{Int, String, Vector}
    Proj(state::Union{Int, String, Vector}) = new(no_signed_zero(state))
end

isless(a::Proj, b::Proj) = isless(state_key(a.state), state_key(b.state))


############ AtIndex ################

"""
    struct AtIndex{R, N} <: IndexedOp{R}

an operator placed on sites, as `X(1)` or `Swap(2, 4)`. The sites must be distinct:
`Swap(1, 1)` and `(X ⊗ Y)(1, 1)` are refused, the product on one site is `X(1) * Y(1)`.
"""
struct AtIndex{R, N} <: IndexedOp{R}
    op::GenericOp{R, N}
    index::NTuple{N, Int}
    function AtIndex(op::GenericOp{R, N}, index::NTuple{N, Int}) where {R, N}
        if !allunique(index)
            error("$op acts on $N sites and cannot be placed on $index, which repeats a site")
        end
        if scalararg(op) isa IdentityOp
            return scalarcoef(op) * IdentityOp{R, Indexed, 1}()
        end
        return scalarcoef(op) * new{R, N}(scalararg(op), index)
    end
end

(op::GenericOp{R, N})(index::Vararg{Int, N}) where {R, N} =
    AtIndex(op, index)

show(io::IO, ind::AtIndex) =
    show_func(io, repr(ind.op; context=:precedence=>500), collect(ind.index))

isless(a::AtIndex, b::AtIndex) =
    isless((a.index, a.op), (b.index, b.op))


############ ComOp ################

"""
    struct ComOp{R} <: IndexedOp{R}

a block of channels of an MPO, the form `compact` gives to the terms of several sites of an
operator: ``\\sum_{i<j} C_i \\left(\\prod_{i<l<j} A_l\\right) B_j``, where on each site the row
``C`` opens channels, the matrix ``A`` carries them and the column ``B`` closes them, each
entry an operator of one site. `pieces` holds the entries of each site from `start` on, as
`(l, r, op)`: `l` a channel of the link on the left of the site, or 0 for the terms not yet
begun, `r` a channel of the link on its right, or 0 for the terms finished. `linkdims` gives
the number of channels of each link, from the right of the first site to the left of the last.
It prints as `com(sites,linkdims)`.
"""
struct ComOp{R} <: IndexedOp{R}
    start::Int
    linkdims::Vector{Int}
    pieces::Vector{Vector{Tuple{Int, Int, GenericOp{R, 1}}}}
end

"""
    com_sites(a)

the sites the com `a` acts on, from the first to the last.
"""
com_sites(a::ComOp) = a.start:a.start + length(a.pieces) - 1

"""
    map_pieces(f, R, a)

the com `a` with `f` applied to each of its pieces, a com of the representation `R`
"""
map_pieces(f, ::Type{R}, a::ComOp) where R =
    ComOp{R}(a.start, a.linkdims, [ [ (l, r, f(o)) for (l, r, o) in p ] for p in a.pieces ])

show(io::IO, a::ComOp) = print(io, "com(", com_sites(a), ",[", join(a.linkdims, ","), "])")

isless(a::ComOp, b::ComOp) =
    isless((a.start, a.linkdims, a.pieces), (b.start, b.linkdims, b.pieces))


############## Mixers ###############

# Dissipator

"""
    Dissipator(L)

the Lindblad dissipator of the jump operator `L`,
``\\rho \\mapsto L\\rho L^\\dagger - \\frac{1}{2}\\{L^\\dagger L, \\rho\\}``, to be added to
the evolver of a mixed state. `Dissipator(c * L)` is `abs2(c) * Dissipator(L)`. `L` is a
generic operator, placed on its sites afterwards, `Dissipator(C)(3)`: `Dissipator(C(3))` is
refused.

# Examples
    Dissipator(Sp)
    Dissipator(0.1 * C)
    Dissipator(Z ⊗ Z)
    Dissipator(0.9X⊗Sm⊗X + 0.1Y⊗Sm⊗Y)
"""
struct Dissipator{N} <: GenericOp{Mixed, N}
    arg::GenericOp{Pure, N}
    Dissipator(arg::GenericOp{Pure, N}) where N =
        abs2(scalarcoef(arg)) * new{N}(scalararg(arg))
end

show(io::IO, a::Dissipator) =
    paren(io, 1000, 0) do io
        show_func(io, "Dissipator", a.arg)
    end

isless(a::Dissipator, b::Dissipator) = isless(a.arg, b.arg)

# Evolver

"""
    Evolver(A)

the superoperator ``\\rho \\mapsto A\\rho + \\rho A^\\dagger`` of a placed operator `A` on
pure states, which for `A = -im * H` is the hamiltonian part ``-i[H, \\rho]`` of an evolver.
It is rarely written: a placed pure operator added to a mixed one, or turned into an MPO for a
mixed state, becomes its `Evolver`.
"""
struct Evolver <: IndexedOp{Mixed}
    arg::IndexedOp{Pure}
end

# the converse of Gate, which takes an operator before it is placed
Evolver(a::GenericOp) =
    error("Evolver takes an operator placed on sites: write Evolver(X(1)) rather than " *
          "Evolver(X)(1)")

show(io::IO, a::Evolver) =
    paren(io, 1000, 0) do io
        show_func(io, "Evolver", a.arg)
    end

isless(a::Evolver, b::Evolver) = isless(a.arg, b.arg)

(a::IndexedOp{Mixed} + b::IndexedOp{Pure}) = a + Evolver(b)
(a::IndexedOp{Pure} + b::IndexedOp{Mixed}) = Evolver(a) + b

"""
    map_sites(f, op)

the operator placed on sites `op` with each of its factors moved from its sites to their
images by `f`, a function from a site to a site: `map_sites(i -> 2i - 1, X(1) * Y(2))` is
`X(1) * Y(3)`. The factors keep their order, and so does the product of fermionic operators
they make. A representation of one's own whose tensors lie on another system, a purification
interleaving an ancilla with each site for instance, places with it on that system the
operators it is given.
"""
map_sites(f, a::AtIndex) = AtIndex(a.op, map(f, a.index))
map_sites(f, a::SumOp) = SumOp(map(x -> map_sites(f, x), a.subs))
map_sites(f, a::ProdOp) = ProdOp(map(x -> map_sites(f, x), a.subs))
map_sites(f, a::ScalarOp) = a.coef * map_sites(f, a.arg)
map_sites(f, a::Evolver) = Evolver(map_sites(f, a.arg))
map_sites(_, a::IdentityOp) = a

# Left

"""
    Left(A)

the superoperator ``\\rho \\mapsto A\\rho``, acting on the left of the density matrix. `A` is a
generic operator, placed on its sites afterwards, `Left(X)(1)`, or an operator already placed,
a sum or a product, `Left(X(1) * Z(3))` being `Left(X)(1) * Left(Z)(3)`, the strings of its
fermionic factors included.

`Gate`, `Dissipator` and `Evolver` act on both sides at once, as ``A\\rho A^\\dagger`` or
``-i[H, \\rho]``. `Left` and `Right` act on a single side, which is not trace preserving. They
can be used in an evolution or applied as gates, but they are not observables.

# Examples

    apply(Left(X)(1), rho)          # rho -> X rho
    make_mpo(rho, Left(X)(1))
    make_mpo(rho, Left(dag(C)(1) * C(3) + N(2)))
"""
struct Left{N} <: GenericOp{Mixed, N}
    arg::GenericOp{Pure, N}
    Left(arg::GenericOp{Pure, N}) where N =
        scalarcoef(arg) *
        (scalararg(arg) isa IdentityOp ? IdentityOp{Mixed, Generic, N}() : new{N}(scalararg(arg)))
end

show(io::IO, a::Left) =
    paren(io, 1000, 0) do io
        show_func(io, "Left", a.arg)
    end

isless(a::Left, b::Left) =
    isless(a.arg, b.arg)

# Right

"""
    Right(A)

the superoperator ``\\rho \\mapsto \\rho A^\\dagger``, acting on the right of the density
matrix. `Right(c * A)` is `conj(c) * Right(A)`. See `Left`, of which it is the mirror, and
which also takes a generic operator or an operator already placed.

# Examples

    apply(Right(X)(1), rho)         # rho -> rho X†
"""
struct Right{N} <: GenericOp{Mixed, N}
    arg::GenericOp{Pure, N}
    Right(arg::GenericOp{Pure, N}) where N =
        conj(scalarcoef(arg)) *
        (scalararg(arg) isa IdentityOp ? IdentityOp{Mixed, Generic, N}() : new{N}(scalararg(arg)))
end

show(io::IO, a::Right) =
    paren(io, 1000, 0) do io
        show_func(io, "Right", a.arg)
    end

isless(a::Right, b::Right) =
    isless(a.arg, b.arg)

"""
    sided(S, a)

the superoperator `S` (`Left` or `Right`) of an operator `a` on pure states. A generic
operator gives `S(a)`. A placed one is built factor by factor, the strings `Multi_F` included:
`Left` is linear and `Right` conjugates the coefficients, and both are multiplicative in the
order of the factors, `Right(A) Right(B) ρ` being `ρ B† A† = Right(A B) ρ`. There is no
fallback for placed forms: one with no method raises rather than being dropped.
"""
sided(S, a::GenericOp{Pure}) = S(a)
sided(S, a::AtIndex{Pure}) = AtIndex(S(a.op), a.index)
sided(S, a::SumOp{Pure, Indexed, 1}) = SumOp(map(x -> sided(S, x), a.subs))
sided(S, a::ProdOp{Pure, Indexed, 1}) = ProdOp(map(x -> sided(S, x), a.subs))
sided(::Type{Left}, a::ScalarOp{Pure, Indexed, 1}) = a.coef * sided(Left, a.arg)
sided(::Type{Right}, a::ScalarOp{Pure, Indexed, 1}) = conj(a.coef) * sided(Right, a.arg)
sided(_, ::IdentityOp{Pure, Indexed, 1}) = IdentityOp{Mixed, Indexed, 1}()
sided(::Type{Left}, a::Multi_F{Pure}) = Multi_F{Mixed}(a.start, a.stop, true, false)
sided(::Type{Right}, a::Multi_F{Pure}) = Multi_F{Mixed}(a.start, a.stop, false, true)
sided(S, a::ComOp{Pure}) = map_pieces(S, Mixed, a)

# an operator already placed, a sum or a product with its coefficients, through sided
Left(a::IndexedOp{Pure}) = sided(Left, a)
Right(a::IndexedOp{Pure}) = sided(Right, a)

# Gate

"""
    Gate(A)

the superoperator ``\\rho \\mapsto A\\rho A^\\dagger`` of an operator `A` on pure states: `A`
applied as a gate to a density matrix. Gates combine linearly, which is how a noisy gate is
written, and `Gate(c * A)` is `abs2(c) * Gate(A)`. `A` is not placed on sites, the gate is:
`Gate(X)(1)` rather than `Gate(X(1))`, and `Gate(X ⊗ Z)(1, 2)` for several sites. A placed
operator on pure states applied to a mixed state, or multiplied by an operator on mixed states,
is turned into its gate.

# Examples

    G = 0.9 * Gate(Id) + 0.1 * Gate(X)      # X with probability 0.1
"""
struct Gate{N} <: GenericOp{Mixed, N}
    arg::GenericOp{Pure, N}
    Gate(arg::GenericOp{Pure, N}) where N =
        abs2(scalarcoef(arg)) *
        (scalararg(arg) isa IdentityOp ? IdentityOp{Mixed, Generic, N}() : new{N}(scalararg(arg)))
end

Gate(a::IndexedOp) =
    error("cannot take the gate of $a, which is placed on sites: write Gate(X)(1) rather than Gate(X(1))")

"""
    build_gate(a)

the gate of the operator `a`, placed on sites, which `apply` takes on a mixed state and a
product of `a` by an operator on mixed states takes
"""
build_gate(a::ProdOp{Pure, Indexed, 1}) = ProdOp(build_gate.(a.subs))
build_gate(ind::AtIndex{Pure}) = AtIndex(Gate(ind.op), ind.index)
build_gate(::IdentityOp{Pure, Indexed, 1}) = IdentityOp{Mixed, Indexed, 1}()
build_gate(a::ScalarOp{Pure, Indexed, 1}) = abs2(a.coef) * build_gate(a.arg)

# a gate built from a sum does not distribute: (A + B) rho (A + B)' has cross terms. On one
# site the sum is an operator of that site, whose gate is placed whole; otherwise the gate of K
# is Left(K) Right(K), which the placed factors of K give one by one
function build_gate(a::SumOp{Pure, Indexed, 1})
    s = scalararg.(a.subs)
    if all(x -> x isa AtIndex && length(x.index) == 1, s) && allequal(x -> x.index, s)
        return Gate(sum(scalarcoef(x) * scalararg(x).op for x in a.subs))(only(first(s).index))
    end
    return ProdOp([sided(Left, a), sided(Right, a)])
end

# a com is a sum, see above: what asked for its gate refuses it with its own message, a product
# as compacted and apply as a sum
build_gate(a::ComOp{Pure}) = ProdOp([sided(Left, a), sided(Right, a)])

(a::IndexedOp{Mixed} * b::IndexedOp{Pure}) = a * build_gate(b)
(a::IndexedOp{Pure} * b::IndexedOp{Mixed}) = build_gate(a) * b

show(io::IO, a::Gate) =
    paren(io, 1000, 0) do io
        show_func(io, "Gate", a.arg)
    end

isless(a::Gate, b::Gate) = isless(a.arg, b.arg)

"""
    lindblad_terms(evolver)

the hamiltonian and the jump operators of an evolver written `-im * H + Σ Dissipator(L)(sites)`,
as `(H, jumps)`: `H` an operator placed on sites, zero when there is none, and `jumps` a vector
of pairs `L => sites`, `L` a generic operator with its rate taken in, so that
`Dissipator(L)(sites...)` is the term again: `γ * Dissipator(Sm)(2)` gives
`sqrt(γ) * Sm => (2,)`. An evolver on pure states is a hamiltonian part alone.

What it reads is what a representation of a state needs to unravel the evolution, quantum
trajectories for instance, from the evolver of an `Evolve` phase. A term of another form, as
`Left(X)(1)`, or a dissipator whose rate is negative or complex, which no jump operator gives, is
refused. The terms of a time dependent evolver are read one by one: a rate multiplied by a
function ``f(t)`` gives the jump ``\\sqrt{f(t)}\\,L``.

# Examples

    H, jumps = lindblad_terms(-im * sum(Z(i) for i in 1:4) + 0.5 * sum(Dissipator(Sm)(i) for i in 1:4))
"""
function lindblad_terms(evolver::IndexedOp)
    hs = IndexedOp{Pure}[]
    jumps = Pair{GenericOp{Pure}, Tuple}[]
    c0 = scalarcoef(evolver)
    for t in sumsubs(scalararg(evolver))
        c, a = c0 * scalarcoef(t), scalararg(t)
        if a isa IndexedOp{Pure}
            push!(hs, im * c * a)
        elseif a isa Evolver && isreal(c)
            push!(hs, im * real(c) * a.arg)
        elseif a isa AtIndex{Mixed}
            for u in sumsubs(a.op)
                r, d = c * scalarcoef(u), scalararg(u)
                if !(d isa Dissipator)
                    error("lindblad_terms reads -im * H and dissipators, and $evolver holds $(d)")
                elseif !isreal(r) || real(r) < 0
                    error("a dissipator of rate $r has no jump operator: the rate must be real " *
                          "and not negative")
                end
                push!(jumps, sqrt(real(r)) * d.arg => a.index)
            end
        else
            error("lindblad_terms reads -im * H and dissipators, and $evolver holds $(c * a)")
        end
    end
    h = isempty(hs) ? 0 * IdentityOp{Pure, Indexed, 1}() : sum(hs)
    return (h, jumps)
end


# SetState

"""
    SetState(state)

the superoperator resetting a site to `state`, given by its name, its vector or its density
matrix: ``\\rho \\mapsto \\sigma \\otimes \\mathrm{tr}_i \\rho``, where ``\\sigma`` is
the density matrix of `state` and ``\\mathrm{tr}_i`` the trace over the site. A vector or a
matrix is taken as it is, without normalization, as for `Proj`. It acts on mixed states only,
and is refused on a site conserving something strongly, which it does not preserve: `weaken`
the state first.

# Examples

    SetState("Up")(3)
    SetState([0.3 0. ; 0. 0.7])(3)
"""
struct SetState <: GenericOp{Mixed, 1}
    state::Union{String, Vector, Matrix}
    SetState(state::Union{String, Vector, Matrix}) = new(no_signed_zero(state))
end

isless(a::SetState, b::SetState) = isless(state_key(a.state), state_key(b.state))


############## Operator functions ###########

# powers

"""
    is_involution(op)

whether an operator is its own inverse by its type alone: the identity, `F`, or an `Operator`
declared `involution_op`. Its integer powers then reduce to it or to the identity.
"""
is_involution(::IdentityOp) = true
is_involution(::JW_F) = true
is_involution(a::Operator) = a.type == involution_op
is_involution(::Op) = false

"""
    is_natural(p)

whether the exponent `p` is an integer of zero or above, whatever its type (`2.0` and
`2 + 0im` are), which makes the power a product
"""
is_natural(p::Number) = isreal(p) && isinteger(real(p)) && real(p) ≥ 0

"""
    struct IntPowOp{R, N} <: GenericOp{R, N}
    struct GenPowOp{R, N} <: GenericOp{R, N}

internal types for the powers of an operator.

- `IntPowOp`, an integer exponent of zero or above: `A^3` is the product `A * A * A`, kept
  whole only to print that way, with the adjoint, parity and Jordan-Wigner strings of that
  product. A power of an involution is reduced at once, `X^10000` is `Id`.
- `GenPowOp`, any other exponent (not an integer, negative or complex): a function of the
  operator, the principal power taken through its logarithm, as `exp` is. It does not split
  over factors, and of a coefficient only the modulus comes out, the phase staying inside.

Both merge, `A^p * A^q` being `A^(p + q)`.
"""
struct IntPowOp{R, N} <: GenericOp{R, N}
    arg::GenericOp{R, N}
    expo::Int
    function IntPowOp(arg::GenericOp{R, N}, n::Int) where {R, N}
        c = scalarcoef(arg)
        a = scalararg(arg)
        if n == 0 || a isa IdentityOp
            return c^n * IdentityOp{R, Generic, N}()
        elseif is_involution(a)
            return c^n * (iseven(n) ? IdentityOp{R, Generic, N}() : a)
        elseif n == 1
            return arg
        elseif a isa IntPowOp
            return c^n * IntPowOp(a.arg, a.expo * n)
        end
        return c^n * new{R, N}(a, n)
    end
end

struct GenPowOp{R, N} <: GenericOp{R, N}
    arg::GenericOp{R, N}
    expo::Number
    function GenPowOp(arg::GenericOp{R, N}, p::Number) where {R, N}
        # a superoperator holding a fermionic operator takes, once placed, strings on both
        # sides of the density matrix, which its power does not commute with: Gate(C +
        # dag(C))^0.5 on site 3 was the power of the gate without them. Refused whatever its
        # sites, rather than on all but the first. `hasfermionic` comes further on, reading the
        # parity of this type
        if R === Mixed && hasfermionic(arg)
            error("$arg has no non integer power: the Jordan-Wigner strings of its fermionic " *
                  "operators do not commute with it")
        end
        c = scalarcoef(arg)
        a = scalararg(arg)
        m = abs(c)
        if a isa IdentityOp
            return complex(c)^p * a
        elseif !iszero(m) && !isone(m)
            # the modulus commutes with everything and leaves the logarithm exactly, the
            # phase would turn the eigenvalues across its cut and stays inside
            return m^p * GenPowOp((c / m) * a, p)
        elseif is_involution(a) && isreal(p) && isone(c)
            # A² = 1 leaves the power modulo 2, which the principal logarithm agrees with
            q = mod(real(p), 2)
            return is_natural(q) ? IntPowOp(a, Int(q)) : new{R, N}(a, no_signed_zero(q))
        end
        return new{R, N}(arg, no_signed_zero(p))
    end
end

"""
    power(a, p)

`a^p`: an `IntPowOp` when `p` is natural, see `is_natural`, a `GenPowOp` otherwise. It is
also what two powers of the same operator merge into, the sum of their exponents choosing
the kind.
"""
power(a::GenericOp, p::Number) = is_natural(p) ? IntPowOp(a, Int(real(p))) : GenPowOp(a, p)

(a::GenericOp ^ p::Number) = power(a, p)

# a literal exponent goes through literal_pow, which for a negative one calls inv, which an
# operator has no method of: X^-1 failed where p = -1; X^p did not
Base.literal_pow(::typeof(^), a::Op, ::Val{p}) where p = a ^ p

# a placed operator is its operator on its sites, which takes any power, and an integer power of
# anything placed is its product; a function of a placed sum or product has nowhere to go
(a::AtIndex ^ p::Number) = power(a.op, p)(a.index...)

# the identity to any power, the principal one of 1 being 1
(a::IdentityOp{R, Indexed} ^ ::Number) where R = a

(a::IndexedOp ^ p::Number) =
    if !is_natural(p)
        error("$a has no power $p once placed: take the power of the operator before placing it")
    elseif iszero(p)
        IdentityOp(a)
    else
        prod(fill(a, Int(real(p))))
    end

"""
    sqrt(::GenericOp)

the principal square root of a generic operator, that is `op^0.5`
"""
sqrt(a::GenericOp) = a ^ 0.5

show(io::IO, a::Union{IntPowOp, GenPowOp}) =
    paren(io, Base.operator_precedence(:^)) do io
        # the base a precedence higher, ^ being right associative: (X^0.5)^0.5 printed
        # X^0.5^0.5, which reads X^(0.5^0.5). A complex or rational exponent in parentheses,
        # X^(0.0 + 0.5im) and X^(1//2), not X^0.0 + 0.5im and X^1//2
        print(IOContext(io, :precedence => Base.operator_precedence(:^) + 1), a.arg)
        print(io, "^", a.expo isa Union{Complex, Rational} ? "($(a.expo))" : a.expo)
    end

isless(a::IntPowOp, b::IntPowOp) = isless((a.arg, a.expo), (b.arg, b.expo))
# an exponent may be complex, which has no order: its two parts are compared
isless(a::GenPowOp, b::GenPowOp) =
    isless((a.arg, real(a.expo), imag(a.expo)), (b.arg, real(b.expo), imag(b.expo)))

# ExpOp

"""
    struct ExpOp{N} <: GenericOp{Pure, N}

the exponential of a generic operator on pure states, see `exp`.
"""
struct ExpOp{N} <: GenericOp{Pure, N}
    arg::GenericOp{Pure, N}
end

"""
    exp(::GenericOp{Pure})

the exponential ``e^A`` of a generic operator on pure states. That of a fermionic operator has
no definite parity, and `isfermionic` refuses it.
"""
exp(a::GenericOp{Pure}) = ExpOp(a)

show(io::IO, a::ExpOp) =
    paren(io, 1000, 0) do io
        show_func(io, "exp", a.arg)
    end

isless(a::ExpOp, b::ExpOp) = isless(a.arg, b.arg)

# DagOp

"""
    struct DagOp{N} <: GenericOp{Pure, N}

the adjoint of a generic operator on pure states, see `dag`. `dag(dag(A))` is `A`, and
`dag(c * A)` is `conj(c) * dag(A)`.
"""
struct DagOp{N} <: GenericOp{Pure, N}
    arg::GenericOp{Pure, N}
    DagOp(arg::DagOp) = arg.arg          # dag is an involution
    DagOp(arg::ScalarOp{Pure, Generic}) =
        conj(arg.coef) * DagOp(arg.arg)
    DagOp(arg::GenericOp{Pure, N}) where N =
        new{N}(arg)
end

"""
    dag(::GenericOp{Pure})

the adjoint ``A^\\dagger`` of a generic operator on pure states.

# Examples

    dag(C) * C
    Sp ⊗ dag(Sp) + dag(Sp) ⊗ Sp
"""
dag(a::GenericOp{Pure}) = DagOp(a)

show(io::IO, a::DagOp) =
    paren(io, 1000, 0) do io
        show_func(io, "dag", a.arg)
    end

isless(a::DagOp, b::DagOp) = isless(a.arg, b.arg)

# ModOp

"""
    struct ModOp{N} <: GenericOp{Pure, N}

an operator `A` taken modulo an integer `m` of at least 2, ``e^{2i\\pi A/m}``, see `mod`. A
conserved quantity written this way is conserved modulo `m`.
"""
struct ModOp{N} <: GenericOp{Pure, N}
    arg::GenericOp{Pure, N}
    modulus::Int
    ModOp(arg::GenericOp{Pure, N}, modulus::Int) where N =
        if modulus < 2
            error("a modulus is at least 2, got $modulus")
        else
            new{N}(arg, modulus)
        end
end

"""
    mod(A, m)

the operator ``e^{2i\\pi A/m}``, of eigenvalues ``e^{2i\\pi a/m}`` for the eigenvalues ``a``
of `A`. It can be measured like any other, and it is how a conserved quantity of
``\\mathbb{Z}_m`` is written.

The modulus is carried by the operator rather than read from its eigenvalues, which cannot
tell it for ``m = 2``: ``\\pm 1`` could as well be integer charges.

# Examples

    mod(N, 3)                       # a Z3 charge
    Boson(6, conserve = mod(N, 3))
"""
mod(a::GenericOp{Pure}, m::Int) = ModOp(a, m)

"""
    parity(A)

the parity operator ``(-1)^A``, that is `mod(A, 2)`: measured, it gives the usual parity, and
conserved, a charge of ``\\mathbb{Z}_2``.

# Examples

    parity(N)
    Fermion(conserve = parity(N))
"""
parity(a::GenericOp{Pure}) = ModOp(a, 2)

show(io::IO, a::ModOp) =
    paren(io, 1000, 0) do io
        if a.modulus == 2
            show_func(io, "parity", a.arg)
        else
            show_func(io, "mod", [a.arg, a.modulus])
        end
    end

isless(a::ModOp, b::ModOp) = isless((a.arg, a.modulus), (b.arg, b.modulus))

"""
    fermion_parity(a, strung)

the parity of an operator: `0` if even, `1` if odd, `nothing` if it has none. An operator of
several sites has that of the product of its pieces placed: a tensor product the sum of the
parities of its factors, an operator defined by an expression that of its expression. It
serves `jw_parity`, and `isfermionic` and `hasfermionic`, which ask two different questions,
told apart by `strung`:

- `strung = true`, for `jw_parity`: how the operator behaves when the `F` of its site crosses
  it, the Jordan-Wigner strings being in place. A `JW` transform is then odd, and a projector
  on a vector is taken to have no parity, since the vector may mix even and odd states.
- `strung = false`, for `isfermionic` and `hasfermionic`: whether the operator holds a fermionic factor whose
  string `simplify` has not inserted yet. A `JW` transform then counts as even, its string
  being already in place, and so does a projector, which never takes a string.

`strung` changes nothing else: a composite operator passes it on to its pieces.
"""
fermion_parity(::Op, ::Bool) = 0
fermion_parity(::JW, strung::Bool) = strung ? 1 : 0
# an operator of several sites cannot be declared fermionic: one defined by an expression has
# its parity, one defined by a matrix or a function is taken as even, its matrix being laid
# with no Jordan-Wigner string
fermion_parity(a::Operator{N}, strung::Bool) where N =
    if a.type == fermionic_op
        1
    elseif N > 1 && a.expr isa Op
        fermion_parity(a.expr, strung)
    else
        0
    end
fermion_parity(a::Union{ScalarOp, DagOp}, strung::Bool) = fermion_parity(a.arg, strung)
# a projector on a basis state, given by its index or by a name, is even: the named states of
# the fermionic sites are all basis states, which a site defined outside the package is taken
# to follow
fermion_parity(a::Proj, strung::Bool) = strung && a.state isa Vector ? nothing : 0

# a tensor product placed is the product of its factors placed, see ⊗
function fermion_parity(a::Union{ProdOp, TensorOp}, strung::Bool)
    ps = map(x -> fermion_parity(x, strung), a.subs)
    return any(isnothing, ps) ? nothing : mod(sum(ps), 2)
end

function fermion_parity(a::SumOp, strung::Bool)
    ps = unique(map(x -> fermion_parity(x, strung), a.subs))
    return length(ps) == 1 ? only(ps) : nothing
end

function fermion_parity(a::IntPowOp, strung::Bool)
    p = fermion_parity(a.arg, strung)
    return isnothing(p) ? nothing : mod(p * a.expo, 2)
end

# a function of an odd operator, its exponential or a non integer power, mixes the two parities
fermion_parity(a::Union{ExpOp, ModOp, GenPowOp}, strung::Bool) =
    fermion_parity(a.arg, strung) == 0 ? 0 : nothing

################## isfermionic #################

"""
    isfermionic(::SimpleOp)

whether an operator of one site on pure states is odd under the fermion parity. An operator
of no definite parity raises an error: a sum of fermionic and non fermionic terms, or the
exponential, `mod` or non integer power of a fermionic operator. For operators of several
sites and superoperators, see `hasfermionic`.
"""
function isfermionic(a::SimpleOp)
    p = fermion_parity(a, false)
    if isnothing(p)
        error("$a has no definite fermionic parity: it sums fermionic and non fermionic " *
              "operators, or is a function of a fermionic one")
    end
    return p == 1
end

# named

"""
    named_type(op)

the `OpType` `named` gives by default: the type of `op` for an `Operator`, `fermionic_op` for
a fermionic expression of one site, `plain_op` otherwise.

A renamed operator of one site keeps its name through `simplify`, so its type is the only thing
that gives it a Jordan-Wigner string. `controlled_type` does the opposite for a fermionic
target: a controlled operator acts on several sites, and `simplify` replaces it by its
expression.
"""
named_type(a::Operator) = a.type
named_type(a::Op) = a isa SimpleOp && isfermionic(a) ? fermionic_op : plain_op

"""
    named(def, name[, sites...]; type)

the operator `name` defined by `def`, an expression, a matrix or a function of the sites: the
way to define an operator of one's own.

Its sites:
- none given: an expression acts on as many sites as it does, a matrix or a function on one;
- given: the operator is computed on them once and for all, which lets an operator of several
  sites into a hamiltonian, see `Operator{N}(name, def, type, sites...)`. A single site
  stands for as many identical ones as the size of a matrix asks for. A function then serves
  these sites only: to keep it for every site, give it without sites and with its `type`.

Its type, see `OpType`, unless `type` gives it:
- an expression without sites: the type of the operator when it is a single one, else
  `fermionic_op` when it is fermionic of one site, `plain_op` otherwise;
- a matrix, or anything given with its sites: the strongest type its matrix satisfies,
  `involution_op`, `selfadjoint_op` or `plain_op`, or `fermionic_op` when it acts on a single
  site, given with it, and anticommutes with its `F`. A matrix given without sites is never
  taken for fermionic: give it `type = fermionic_op`;
- a function without sites: `plain_op`.

A matrix given without sites is checked against its type at once, except for its parity,
which needs the `F` of a site and is checked when the matrix is laid on one. A type given to a
function or an expression without sites is taken as it is: `simplify` relies on it, squaring an
involution to the identity for instance, so a wrong one gives a wrong result.

The name also identifies a conserved quantity: sites conserving operators of the same name
share one charge, and renaming keeps two quantities apart.

# Examples

    named([1 1 ; 1 -1] / √2, "MyH")             # involution_op
    named(Sp + Sm, "Sx2", Spin(1))               # selfadjoint_op, read on Spin(1)
    named(swap_matrix, "MySwap", Qubit())        # two qubits, from the size of the matrix
    named(s -> ..., "K"; type = selfadjoint_op)
    Fermion(conserve = named(N, "Nf"))    # these two numbers are
    Boson(4, conserve = named(N, "Nb"))   # conserved separately
"""
named(op::GenericOp{Pure, N}, name::String; type::OpType = named_type(op)) where N =
    Operator{N}(name, op, type)


"""
    hasfermionic(::Op)

whether an operator still holds an odd operator of one site whose Jordan-Wigner string
`simplify` has not inserted, wherever it sits: in a tensor product, in the expression of an
operator of several sites, or inside a superoperator. The tensor of such an operator has
legs on the sites it acts on only, so that it cannot hold the string on the sites in between:
`(C ⊗ dag(C))(1, 3)`, which is `C(1) * dag(C)(3)`, would miss the string on site 2.

After `simplify` the factors are `JW` transforms, which are not fermionic, and the answer is
false unless a factor was kept whole, as an exponential of several sites is. A factor of no
definite parity, as `C + N`, answers true instead of raising.
"""
hasfermionic(a::GenericOp{Pure, 1}) = fermion_parity(a, false) ≠ 0

function hasfermionic(a::Op)
    for f in fieldnames(typeof(a))
        x = getfield(a, f)
        if x isa Op && hasfermionic(x)
            return true
        elseif x isa Vector && any(y -> y isa Op && hasfermionic(y), x)
            return true
        end
    end
    return false
end


"""
    jw_parity(a)

how a factor of a product behaves when the `F` of its sites cross it: `0` if it commutes with
them, `1` if it anticommutes, `nothing` if neither.

A fermionic operator and its Jordan-Wigner transform are odd, any other operator is taken as
even, the convention the strings rest on, and a composite factor gets the parity of its pieces.
Sums count too, since `simplify` gathers the terms of one site into one factor: taking
`(C + dag(C))(1)` as even gave `C(3) * (C + dag(C))(1)` the wrong sign.
"""
jw_parity(a) = fermion_parity(a, true)

################## Equality #################

# An operator is compared by what it is made of, never by the identity of the object that
# holds it. The default `==` of an immutable struct falls back to `===`, which walks the
# fields but compares a `Vector` field by identity, so the types carrying `subs` need a
# definition of their own — and so does every type that may hold one of them, since a
# `ScalarOp` wrapping a `ProdOp` is `===` only to itself. Picking those types one by one
# is what let `2X(1)*Y(2)`, `dag(X*Y)`, `Left(X*Y)`, `Phase(0.3)`, `controlled(Z)` and
# half the hierarchy fall through, so it is read from the type instead, once, the way
# `phase_hash` reads a phase.
#
# `Set` and `Dict` pick their bucket by `hash` and only then compare, so the two have to
# be defined together or two equal operators land apart. `Measure` relies on this to ask
# for a measurement shared by two observables only once.

(a::Op == b::Op) =
    typeof(a) == typeof(b) &&
    all(f -> getfield(a, f) == getfield(b, f), fieldnames(typeof(a)))

hash(a::Op, h::UInt) =
    foldl((h, f) -> hash(getfield(a, f), h), fieldnames(typeof(a)); init = hash(typeof(a), h))


################## Global Ordering ###############

"""
    ranking(a)

the rank of the type of an operator in the order of all operators: operators of different
types are ordered by it, those of one type by their own `isless`. A type with no rank raises.
"""
ranking(a) = error("ranking not defined for ($a)")

isless(a::Op, b::Op) = isless((ranking(a), a), (ranking(b), b))

ranking(::IdentityOp) = 1
ranking(::JW_F) = 2
ranking(::Operator) = 3
ranking(::JW) = 4
ranking(::Proj) = 5
ranking(::Multi_F) = 6

ranking(::AtIndex) = 10
ranking(::Gate) = 11
ranking(::Dissipator) = 12
ranking(::Evolver) = 13
ranking(::Left) = 14
ranking(::Right) = 15
ranking(::SetState) = 16

ranking(::ScalarOp) = 20
ranking(::ProdOp) = 21
ranking(::SumOp) = 22
ranking(::TensorOp) = 23
ranking(::ComOp) = 24

ranking(::IntPowOp) = 30
ranking(::GenPowOp) = 34
ranking(::ExpOp) = 31
ranking(::DagOp) = 32
ranking(::ModOp) = 33

"""
    obs_name(op)

the default name of a measurement: the compact printed form of the operator, in which long
sums and products are abbreviated. Only an `Operator` has a `name` field; `X * Y`, `2X` or
`X + Y` have no other description than their printed form.
"""
obs_name(op) = sprint(print, op; context = :compact => true)
