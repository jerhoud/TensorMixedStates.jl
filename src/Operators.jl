# The operator algebra: the types of operators, generic or placed on sites, acting on pure
# states or on density matrices, and how they are built, combined, compared and printed, before
# any site or state is involved.

export Representation, Pure, Mixed, GenericOp, IndexedOp, SimpleOp
export OpType, plain_op, fermionic_op, selfadjoint_op, involution_op
export Op, Operator, Id, F, Proj, Basis, Gate, Dissipator, Evolver, Left, Right, SetState, Dephase
export relaxing_dissipator, relaxing_gate, depolarizing_dissipator, depolarizing_gate
export dephasing_dissipator, dephasing_gate
export named, parity
export dag, ⊗, isfermionic, hasfermionic, map_sites

############# Types ################

"""
    abstract type Representation

the supertype of the representations of a state: `Pure`, `Mixed`, and those an extension
defines, whose states are `AbstractState`s.
"""
abstract type Representation end

"""
    abstract type PM <: Representation

the supertype of `Pure` and `Mixed`, which parametrize operators as well as states
"""
abstract type PM <: Representation end

"""
    struct Pure <: PM
    Pure()

the representation of a state as a wave function, as in `State{Pure}` and `Op{Pure}`;
`Pure()` selects it where a representation is passed as a value.

# Examples

    State{Pure}(System(4, Qubit()), "Up")
    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up")
"""
struct Pure <: PM end


"""
    struct Mixed <: PM
    Mixed()

the representation of a state as a density matrix, as in `State{Mixed}` and `Op{Mixed}`;
`Mixed()` selects it where a representation is passed as a value. Superoperators, as `Gate`,
`Dissipator` or `Left`, act on it only.

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

print the real number `x`, a float of integer value as that integer, `2` rather than `2.0`,
below `maxintfloat`
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
`2im*` for an imaginary number, a rational or a complex number in parentheses
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

`x` with its signed zeros made positive, as operators store their numbers: `-0.0 == 0.0`,
but `isless` and `hash` tell them apart, which would keep equal terms apart in `simplify` and
`measure`
"""
no_signed_zero(x::Union{AbstractFloat, Complex{<:AbstractFloat}}) = x + zero(x)
no_signed_zero(x::AbstractArray) = map(no_signed_zero, x)
no_signed_zero(x) = x

"""
    state_key(state)

the key a `Proj` or a `SetState` is ordered by, total, complex numbers included, and tying
exactly when the states are `==`
"""
state_key(x::Int) = (1, x)
state_key(x::AbstractString) = (2, x)
state_key(x::AbstractArray) = (3, size(x), [ (real(y), imag(y)) for y in vec(x) ])
state_key(x::Pair) = (4, x.first, x.second)

############### Operator ###############

"""
    @enum OpType

the possible types of an `Operator`.

# Enumeration values

- `plain_op`: no particular property
- `fermionic_op`: fermionic, the Jordan-Wigner transform applies to it
- `selfadjoint_op`: invariant under `dag`
- `involution_op`: invariant under `dag`, and its square is the identity

`simplify` relies on the type, which is checked against the matrix when the operator is
placed on a site.
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

the identity, one value for each kind of operator, so that every way of writing it, as
`Id ⊗ Id`, `Gate(Id)` or `Id(3)`, is recognized by its type. Placed, it has no site.
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
        # in one pass: a sum built term by term, as `sum` of a generator does, copies its terms
        # at each one, and filtering them again would make that cost grow faster still
        s = Op{R, T, N}[]
        for x in subs
            if x isa SumOp
                # its terms are already those of coefficient other than zero
                append!(s, x.subs)
            elseif scalarcoef(x) ≠ 0
                # left out: the zero operator has every parity, and `0C + dag(C)` would have none
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

the number `x` as the coefficient of an operator: an integer of fixed width, or a complex
number of them, as a float, whose products cannot silently wrap around
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
    struct TensorOp{R, N} <: GenericOp{R, N}

a tensor product of generic operators of representation `R`, acting on `N` sites, the sum of
theirs.
The coefficients of the factors are gathered in front of it, and a product of identities is
the identity of `N` sites.
"""
struct TensorOp{R, N} <: GenericOp{R, N}
    subs::Vector{<:GenericOp{R}}
    TensorOp{R, N}(subs::Vector{<:GenericOp{R}}) where {R, N} =
        if length(subs) == 1
            subs[1]
        else
            c = prod(scalarcoef.(subs))
            s = scalararg.(subs)
            if all(x -> x isa IdentityOp, s)
                return c * IdentityOp{R, Generic, N}()
            end
            c * new{R, N}(s)
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

the ranges of the sites of the tensor product `a` that each of its factors acts on, in order
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

the tensor product of generic operators, acting on the sites of `op1` followed by those of
`op2`: `(A ⊗ B)(i, j)` is `A(i) * B(j)`. Beside a superoperator, an operator on pure states is
its gate, as in a product. `⊗` is typed `\\otimes`, and `tensor(op1, op2, ...)` is the same
product.

# Examples

    mycontrolled(op) = Proj("Up") ⊗ Id + Proj("Dn") ⊗ op     # op of one site
    Rxy(t) = exp(-im * t * (X ⊗ X + Y ⊗ Y) / 4)
    FM = SetState("FullyMixed")
    depolarizing2(p) = (1 - p) * Gate(Id ⊗ Id) + p * (FM ⊗ FM)

"""
(a::GenericOp{R, N} ⊗ b::GenericOp{R, M}) where {R, N, M} =
    TensorOp{R, N + M}([tensorsubs(a) ; tensorsubs(b)])
(a::GenericOp{Pure} ⊗ b::GenericOp{Mixed}) = Gate(a) ⊗ b
(a::GenericOp{Mixed} ⊗ b::GenericOp{Pure}) = a ⊗ Gate(b)
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
    Proj(A => λ)

the projector ``|s\\rangle\\langle s|`` on a pure state of one site, given by its name, by
its vector in the basis of the site, taken as it is without normalization, or by the number
of a basis state, counted from 0; or the projector on the eigenspace of eigenvalue `λ` of the
Hermitian operator `A`, of one site or several, `λ` being refused if it is not an eigenvalue.
On several sites it is applied as a gate, and `named(Proj(A => λ), name, sites...)` writes it
as a sum of products, for a Hamiltonian or a measurement. On a fermionic site, the state must
have a definite parity and `A` must commute with the parity.

# Examples

    Proj("Up")
    Proj([1, 0])
    Proj(1)            # the second basis state
    Proj(Sz => 0)      # on a Spin(1), the state of Sz = 0
    Proj(Ntot => 1)    # on an Electron, Up and Dn
    Proj(Swap => -1)   # on two qubits, the singlet
    named(Proj(Sx⊗Sx + Sy⊗Sy + Sz⊗Sz => 1), "P2", Spin(1))   # spin 2 of two spins 1
"""
struct Proj{N} <: GenericOp{Pure, N}
    state::Union{Int, String, Vector, Pair{<:GenericOp{Pure}, <:Real}}
    Proj(state::Union{Int, String, Vector}) = new{1}(no_signed_zero(state))
    Proj(p::Pair{<:GenericOp{Pure, N}, <:Real}) where N = new{N}(first(p) => no_signed_zero(last(p)))
end

Proj(p::Pair{<:IndexedOp, <:Real}) =
    error("cannot project on an eigenspace of $(first(p)), which is placed on sites: write " *
          "Proj(N => 1)(3) rather than Proj(N(3) => 1)")

# the default would print the number of sites, Proj{1}("Up"), in the names of measurements
show(io::IO, a::Proj) = show_func(io, "Proj", [a.state])

isless(a::Proj, b::Proj) = isless(state_key(a.state), state_key(b.state))


############ Basis ################

"""
    Basis

the operator numbering the basis states of a site, ``\\mathrm{diag}(0, 1, \\dots, d - 1)``,
defined on every site type: its eigenbasis is the basis of the site. On a `Qubit`, a `Fermion`,
a `Boson`, a `Qboson`, a `Qudit` or a `Spin`, it is `N`.

# Examples

    Dephase(Basis)(3)
    collapse(state, Basis(3))
"""
const Basis = Operator{1}("Basis", s -> diagm(Float64.(0:dim(s) - 1)), selfadjoint_op)


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

Dissipator(a::IndexedOp) =
    error("cannot take the dissipator of $a, which is placed on sites: write Dissipator(Sm)(1) " *
          "rather than Dissipator(Sm(1))")

show(io::IO, a::Dissipator) =
    paren(io, 1000, 0) do io
        show_func(io, "Dissipator", a.arg)
    end

isless(a::Dissipator, b::Dissipator) = isless(a.arg, b.arg)

# Evolver

"""
    Evolver(A)

the superoperator ``\\rho \\mapsto A\\rho + \\rho A^\\dagger`` of a placed operator `A` on
pure states, which for `A = -im * H` is the Hamiltonian part ``-i[H, \\rho]`` of an evolver.
It is rarely written: a placed pure operator added to a mixed one, or turned into an MPO for a
mixed state, becomes its `Evolver`.
"""
struct Evolver <: IndexedOp{Mixed}
    arg::IndexedOp{Pure}
end

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
`X(1) * Y(3)`, the factors keeping their order. A representation of one's own whose tensors
lie on another system, a purification for instance, places with it the operators it is given.
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

Acting on a single side, `Left` and `Right` do not preserve the trace. They can be used in an
evolution or applied as gates, but are not observables.

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

the superoperator `S`, `Left` or `Right`, of an operator `a` on pure states, built factor by
factor for a placed one, strings included: both keep the order of the factors,
`Right(A) Right(B) ρ` being `ρ B† A† = Right(A B) ρ`, and `Right` conjugates the coefficients
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

Left(a::IndexedOp{Pure}) = sided(Left, a)
Right(a::IndexedOp{Pure}) = sided(Right, a)

# Gate

"""
    Gate(A)

the superoperator ``\\rho \\mapsto A\\rho A^\\dagger`` of an operator `A` on pure states.
Gates combine linearly, which is how a noisy gate is written, and `Gate(c * A)` is
`abs2(c) * Gate(A)`. The gate is placed, not `A`: `Gate(X)(1)`, `Gate(X ⊗ Z)(1, 2)`. A placed
operator on pure states applied to a mixed state, or multiplied by a mixed one, becomes its
gate.

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

the gate of the placed operator `a`, see `Gate`
"""
build_gate(a::ProdOp{Pure, Indexed, 1}) = ProdOp(build_gate.(a.subs))
build_gate(ind::AtIndex{Pure}) = AtIndex(Gate(ind.op), ind.index)
build_gate(::IdentityOp{Pure, Indexed, 1}) = IdentityOp{Mixed, Indexed, 1}()
build_gate(a::ScalarOp{Pure, Indexed, 1}) = abs2(a.coef) * build_gate(a.arg)

# the gate of a sum does not distribute, (A + B) ρ (A + B)† having cross terms: on one site it
# is the gate of the sum, otherwise Left(K) Right(K)
function build_gate(a::SumOp{Pure, Indexed, 1})
    s = scalararg.(a.subs)
    if all(x -> x isa AtIndex && length(x.index) == 1, s) && allequal(x -> x.index, s)
        return Gate(sum(scalarcoef(x) * scalararg(x).op for x in a.subs))(only(first(s).index))
    end
    return ProdOp([sided(Left, a), sided(Right, a)])
end

# a com is a sum, see above
build_gate(a::ComOp{Pure}) = ProdOp([sided(Left, a), sided(Right, a)])

(a::IndexedOp{Mixed} * b::IndexedOp{Pure}) = a * build_gate(b)
(a::IndexedOp{Pure} * b::IndexedOp{Mixed}) = build_gate(a) * b

show(io::IO, a::Gate) =
    paren(io, 1000, 0) do io
        show_func(io, "Gate", a.arg)
    end

isless(a::Gate, b::Gate) = isless(a.arg, b.arg)


# SetState

"""
    SetState(state)

the superoperator resetting a site to `state`, given by its name, the number of a basis state
counted from 0, its vector or its density matrix: ``\\rho \\mapsto \\sigma \\otimes \\mathrm{tr}_i \\rho``, where ``\\sigma`` is
the density matrix of `state` and ``\\mathrm{tr}_i`` the trace over the site. A vector or a
matrix is taken as it is, without normalization, as for `Proj`. It acts on mixed states only,
and is refused on a site conserving something strongly, which it does not preserve: `weaken`
the state first.

# Examples

    SetState("Up")(3)
    SetState(0)(3)
    SetState([0.3 0. ; 0. 0.7])(3)
"""
struct SetState <: GenericOp{Mixed, 1}
    state::Union{Int, String, Vector, Matrix}
    SetState(state::Union{Int, String, Vector, Matrix}) = new(no_signed_zero(state))
end

isless(a::SetState, b::SetState) = isless(state_key(a.state), state_key(b.state))


# Dephase

"""
    Dephase()
    Dephase(A)

the superoperator dephasing a site in the basis of the eigenspaces of the Hermitian operator
`A`, ``\\rho \\mapsto \\sum_\\lambda P_\\lambda \\rho P_\\lambda``, ``P_\\lambda`` being the
projector on the eigenspace of `A` of eigenvalue ``\\lambda``: the measurement of `A` whose
result is not read. It erases the coherences between eigenspaces and keeps those within one.
Without `A`, it is `Dephase(Basis)`, which dephases in the basis of the site. On a fermionic
site, `A` must commute with the parity.

# Examples

    Dephase()(3)
    Dephase(Ntot)(2)        # on an Electron, keeps the coherence of Up and Dn
"""
struct Dephase <: GenericOp{Mixed, 1}
    arg::GenericOp{Pure, 1}
    Dephase(arg::GenericOp{Pure, 1} = Basis) = new(arg)
end

Dephase(a::IndexedOp) =
    error("cannot dephase in the eigenbasis of $a, which is placed on sites: write " *
          "Dephase(N)(1) rather than Dephase(N(1))")

show(io::IO, a::Dephase) =
    paren(io, 1000, 0) do io
        show_func(io, "Dephase", a.arg)
    end

isless(a::Dephase, b::Dephase) = isless(a.arg, b.arg)


# Relaxation

"""
    tensor_power(a, n)

the tensor product of `n` copies of the generic operator `a`, `a` itself for `n = 1`
"""
tensor_power(a::GenericOp, n::Int) = reduce(⊗, fill(a, n))

"""
    resets(n, state)
    resets(states)

the number of sites and the superoperator resetting them, to `state` each or to a state per
site, a vector of numbers being the amplitudes of a single state, as for `State`
"""
resets(n::Int, state) = (n, tensor_power(SetState(state), n))
resets(state) = resets(1, state)
resets(states::Vector) = (length(states), reduce(⊗, [ SetState(s) for s in states ]))
resets(states::Vector{<:Number}) = resets(1, states)

"""
    relaxing_dissipator(γ, state)
    relaxing_dissipator(γ, n, state)
    relaxing_dissipator(γ, [state1, state2, ...])

the Lindblad generator relaxing a site, `n` sites, or as many sites as states, at rate `γ`
towards `state`, or towards a state per site, for an evolver, on sites of any type:
``\\gamma\\,(\\sigma \\otimes \\mathrm{tr}_S(\\rho) - \\rho)``, see `SetState`, the sites
relaxing together. A state is given as for `SetState`, and a vector of numbers is the amplitudes
of a single state. Evolving under it for a time `t` is `relaxing_gate(1 - exp(-γt), ...)`. It is
refused on sites conserving something strongly.

# Examples

    evolver = -im * H + sum(relaxing_dissipator(0.1, "Up")(i) for i in 1:10)
"""
function relaxing_dissipator(γ::Real, s...)
    n, r = resets(s...)
    return γ * (r - Gate(tensor_power(Id, n)))
end

"""
    relaxing_gate(p, state)
    relaxing_gate(p, n, state)
    relaxing_gate(p, [state1, state2, ...])

the channel relaxing a site, `n` sites, or as many sites as states, with probability `p` towards
`state`, or towards a state per site, for a gate, on sites of any type, the sites relaxing
together: ``(1 - p)\\,\\rho + p\\,\\sigma \\otimes \\mathrm{tr}_S(\\rho)``, the reset of
`SetState` taking place with probability `p`. It is `Gate(Id) + relaxing_dissipator(p, ...)`.

# Examples

    relaxing_gate(0.1, "Up")(3)
    relaxing_gate(0.1, ["Up", "Dn"])(1, 2)
"""
function relaxing_gate(p::Real, s...)
    n, r = resets(s...)
    return (1 - p) * Gate(tensor_power(Id, n)) + p * r
end

"""
    depolarizing_dissipator(γ, n = 1)

the Lindblad generator depolarizing `n` sites together at rate `γ`, for an evolver:
`relaxing_dissipator(γ, n, "FullyMixed")`, ``\\rho \\mapsto \\gamma\\,(\\mathrm{tr}_S(\\rho)
\\otimes I/d^n - \\rho)``. Evolving under it for a time `t` is `depolarizing_gate(1 - exp(-γt),
n)`.

# Examples

    evolver = -im * H + sum(depolarizing_dissipator(0.1)(i) for i in 1:10)
"""
depolarizing_dissipator(γ::Real, n::Int = 1) = relaxing_dissipator(γ, n, "FullyMixed")

"""
    depolarizing_gate(p, n = 1)

the channel depolarizing `n` sites together with probability `p`, for a gate: `relaxing_gate(p,
n, "FullyMixed")`, ``\\rho \\mapsto (1 - p)\\,\\rho + p\\,\\mathrm{tr}_S(\\rho) \\otimes
I/d^n``, on a qubit `(1 - 3p/4) * Gate(Id) + p/4 * (Gate(X) + Gate(Y) + Gate(Z))`. Two sites
depolarized together are not two sites depolarized each on its own, `depolarizing_gate(p) ⊗
depolarizing_gate(p)`.

# Examples

    Gates(gates = prod(depolarizing_gate(0.01)(i) for i in 1:10))
    apply(depolarizing_gate(0.02, 2)(1, 2) * controlled(Z)(1, 2), ρ)
"""
depolarizing_gate(p::Real, n::Int = 1) = relaxing_gate(p, n, "FullyMixed")


# Dephasing

"""
    dephasing_dissipator(γ, A = Basis)

the Lindblad generator dephasing a site at rate `γ` in the eigenbasis of `A`, by default the
basis of the site, for an evolver: ``\\gamma\\,(\\Delta - \\mathrm{Id})``, `Δ` being
`Dephase(A)`. Every coherence between two eigenspaces decays at the same rate `γ`, where
`Dissipator(A)` damps it at ``(\\lambda_a - \\lambda_b)^2 / 2``, faster for distant
eigenvalues: the two agree when `A` has two eigenvalues, `Dissipator(Z)` being
`dephasing_dissipator(2, Z)`. Evolving under it for a time `t` is
`dephasing_gate(1 - exp(-γt), A)`.

# Examples

    evolver = -im * H + sum(dephasing_dissipator(0.1)(i) for i in 1:10)
"""
dephasing_dissipator(γ::Real, a::GenericOp{Pure, 1} = Basis) = γ * (Dephase(a) - Gate(Id))

"""
    dephasing_gate(p, A = Basis)

the channel dephasing a site with probability `p` in the eigenbasis of `A`, by default the
basis of the site, for a gate: ``(1 - p)\\,\\rho + p\\,\\Delta(\\rho)``, `Δ` being
`Dephase(A)`. It is `Gate(Id) + dephasing_dissipator(p, A)`.

# Examples

    Gates(gates = prod(dephasing_gate(0.01)(i) for i in 1:10))
"""
dephasing_gate(p::Real, a::GenericOp{Pure, 1} = Basis) = (1 - p) * Gate(Id) + p * Dephase(a)


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

- `IntPowOp`, an integer exponent of zero or above: the product, kept whole to print that way
- `GenPowOp`, any other exponent: the principal power, a function of the operator, from which
  only the modulus of a coefficient comes out
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
        # placed, a superoperator holding a fermionic operator takes strings on both sides of
        # the density matrix, which its power does not commute with. `hasfermionic` is defined
        # further on
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

`a^p`: an `IntPowOp` when `p` is natural, see `is_natural`, a `GenPowOp` otherwise
"""
power(a::GenericOp, p::Number) = is_natural(p) ? IntPowOp(a, Int(real(p))) : GenPowOp(a, p)

(a::GenericOp ^ p::Number) = power(a, p)

# a literal negative exponent would go through `inv`, which an operator has no method of
Base.literal_pow(::typeof(^), a::Op, ::Val{p}) where p = a ^ p

# the power of a placed operator is that of its operator; anything else placed only takes an
# integer power, its product
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
        # the base a precedence higher, ^ being right associative: (X^0.5)^0.5, not X^0.5^0.5
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
of `A`. Measured, `mod(N, 3)` gives the cube roots of unity; conserved, `conserve = mod(N, 3)`
gives the charges ``n \\bmod 3``, 0, 1 or 2, added modulo 3 over the sites.

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

the parity of an operator: `0` if even, `1` if odd, `nothing` if it has none, an operator of
several sites taking that of the product of its pieces placed. `strung` tells two questions
apart:

- `true`, for `jw_parity`: how the operator behaves when the `F` of its site crosses it, the
  strings in place. A `JW` transform is odd.
- `false`, for `isfermionic` and `hasfermionic`: whether it holds a fermionic factor whose
  string `simplify` has not inserted yet. A `JW` transform is even.
"""
fermion_parity(::Op, ::Bool) = 0
fermion_parity(::JW, strung::Bool) = strung ? 1 : 0
# an operator of several sites defined by a matrix or a function is even, its matrix being laid
# with no string
fermion_parity(a::Operator{N}, strung::Bool) where N =
    if a.type == fermionic_op
        1
    elseif N > 1 && a.expr isa Op
        fermion_parity(a.expr, strung)
    else
        0
    end
fermion_parity(a::Union{ScalarOp, DagOp}, strung::Bool) = fermion_parity(a.arg, strung)
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
a fermionic expression of one site, whose string its type alone gives once renamed, and
`plain_op` otherwise
"""
named_type(a::Operator) = a.type
named_type(a::Op) = a isa SimpleOp && isfermionic(a) ? fermionic_op : plain_op

"""
    named(def, name[, sites...]; type)

the operator `name` defined by `def`, an expression, a matrix or a function of the sites: the
way to define an operator of one's own.

Without sites, an expression acts on as many sites as it does, a matrix or a function on one.
With sites, the operator is computed on them once and for all, a single site standing for as
many as the size of a matrix asks for, and a function then serves these sites only: give it
without sites, with its `type`, to keep it for every site.

Its type, see `OpType`, unless `type` gives it:
- an expression without sites: that of the operator when it is a single one, else
  `fermionic_op` when it is fermionic of one site, `plain_op` otherwise;
- a matrix, or anything given with its sites: the strongest type its matrix satisfies, and
  `fermionic_op` only on a single site given with it; a matrix without sites is never taken
  for fermionic;
- a function without sites: `plain_op`.

A type given to a function or an expression without sites is taken as it is, and a wrong one
gives wrong results. Sites conserving operators of the same name share one charge, and
renaming keeps two quantities apart.

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

whether an operator holds an odd operator of one site whose Jordan-Wigner string `simplify`
has not inserted, wherever it sits: in a tensor product, in the expression of an operator of
several sites, or inside a superoperator. Such an operator needs `simplify` before it is laid,
its tensor holding no string on the sites in between. A factor of no definite parity, as
`C + N`, answers true.
"""
hasfermionic(a::GenericOp{Pure, 1}) = fermion_parity(a, false) ≠ 0

# a projector on an eigenspace of fermionic operators of several sites, whose matrix lacks the
# strings of the sites in between, see `strung_function`
hasfermionic(a::Proj) = a.state isa Pair && hasfermionic(first(a.state))
hasfermionic(::Proj{1}) = false

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

A fermionic operator and its Jordan-Wigner transform are odd, any other operator is even, and
a composite factor, a sum included, gets the parity of its pieces.
"""
jw_parity(a) = fermion_parity(a, true)

################## Equality #################

# An operator is compared by what it is made of: the default `==` of an immutable struct
# compares a `Vector` field by identity, so equality is read from the fields of every type at
# once. `hash` goes with it, or two equal operators land apart in a `Set` or a `Dict`.

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
ranking(::Dephase) = 17

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

the default name of a measurement: the compact printed form of the operator
"""
obs_name(op) = sprint(print, op; context = :compact => true)
