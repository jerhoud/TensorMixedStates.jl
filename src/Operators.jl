export Pure, Mixed, GenericOp, IndexedOp, SimpleOp
export OpType, plain_op, fermionic_op, selfadjoint_op, involution_op
export Op, Operator, Id, F, Proj, Gate, Dissipator, Evolver, Left, Right, SetState
export named, parity
export dag, ⊗, isfermionic, has_fermionic

############# Types ################

"""
    abstract type PM

`PM` is a supertype for `Pure` and `Mixed` (i.e. `Pure()` and `Mixed()` are of type `PM`)
"""
abstract type PM end

"""
    type Pure <: PM
    Pure()

correspond to pure quantum representation
"""
struct Pure <: PM end


"""
    type Mixed <: PM
    Mixed()

correspond to mixed quantum representation
"""
struct Mixed <: PM end

"""
    abstract type GI end
    
`GI` is a supertype for `Generic` and `Indexed`
"""
abstract type GI end

"""
    Generic

a type representing generic operators (without site indices) to parametrize `Op`
"""
struct Generic <: GI end

"""
    Indexed

a type representing indexed operators (with site indices) to parametrize `Op`
"""
struct Indexed <: GI end

"""
    Op{R <: PM, T <: GI, N}

the type of all operators.

# Type parameters

- `R`: is `Pure` or `Mixed` (the type of the representations on which the operator may be applied)
- `T`: is `Generic` or `Indexed`
- `N`: the number of sites on which the operator must be applied (set to 1 for indexed operators)
"""
abstract type Op{R <: PM, T <: GI, N} end

"""
    GenericOp{R, N}

the type of generic operators (without site indices), that is `Op{R, Generic, N}`
"""
const GenericOp{R, N} = Op{R, Generic, N}

"""
    IndexedOp{R}

the type of indexed operators (with site indices), that is `Op{R, Indexed, 1}`
"""
const IndexedOp{R} = Op{R, Indexed, 1}

"""
    SimpleOp

the type of generic pure operators acting on one site, that is `GenericOp{Pure, 1}`

Note that it is abstract: `X` is a `SimpleOp`, but so is any one site combination such as
`2X`, `X * Y` or `exp(X)`.
"""
const SimpleOp = GenericOp{Pure, 1}

############## Showing ###############

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
        else
            print(io, "(", a, ")")
        end
    elseif isa(a, Rational)
        # (1//2)X and not 1//2X, which reads 1//(2X)
        print(io, "(", a, ")")
    else
        print(io, a)
    end
end

# the operands after the first a precedence higher, * and ⊗ being left associative: X ⊗ (Y*Z)
# printed X⊗Y*Z, which reads (X⊗Y)*Z
function infix(io::IO, subs, op::String)
    p = get(io, :precedence, 0)
    for (k, s) in enumerate(subs)
        if k > 1
            print(io, op)
        end
        print(k == 1 ? io : IOContext(io, :precedence => p + 1), s)
    end
end

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


############### Operator ###############

"""
    @enum OpType

the possible operator types for `Operator`

# Enumeration values

- `plain_op`: an operator with no particular properties
- `fermionic_op`: a fermionic operator for which Jordan-Wigner transform must be used
- `selfadjoint_op`: an operator invariant under `dag`
- `involution_op`: an operator invariant under `dag` and whose square is the identity 

`simplify` reasons with the type, taking `dag(X)` to be `X` or `X * X` to be `Id`, and moving
the `F` of a site across an operator with a sign for a fermionic one and without for any other.
The type is therefore checked against the matrix each time the operator is placed on a site,
and a type the matrix belies is refused rather than giving a wrong result.
"""
@enum OpType plain_op fermionic_op selfadjoint_op involution_op

"""
    plain_op

the `OpType` of an operator with no particular properties
"""
plain_op

"""
    fermionic_op

the `OpType` of a fermionic operator, for which the Jordan-Wigner transform must be used
"""
fermionic_op

"""
    selfadjoint_op

the `OpType` of an operator invariant under `dag`
"""
selfadjoint_op

"""
    involution_op

the `OpType` of an operator invariant under `dag` and whose square is the identity
"""
involution_op

"""
    type Operator{N} <: GenericOp{Pure, N}

the type of base operators (like `X`, `Swap`, `C` ...),
`N` is the number of sites on which it may be applied.

# Example
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
            new{N}(name, expr, type)
        end
end

show(io::IO, op::Operator) =
    print(io, op.name)

isless(a::Operator, b::Operator) = isless(a.name, b.name)


############ Identity ############

"""
    type IdentityOp{R, T, N}

the identity, a single value for each kind of operator: pure or on a density matrix, generic
on `N` sites, or placed. It is what every construction of an identity gives, `Id ⊗ Id`,
`Left(Id)`, `Right(Id)`, `Gate(Id)` and `Id(3)` included, so that an identity is told by its
type alone. Placed, it has no site: it is the identity of the whole system, and a tensor of it
on a given site is laid by the code that needs one, which knows the site.
"""
struct IdentityOp{R, T, N} <: Op{R, T, N} end

IdentityOp(::Op{R, T, N}) where {R, T, N} = IdentityOp{R, T, N}()

"""
    Id

the identity operator defined for all site types
"""
const Id = IdentityOp{Pure, Generic, 1}()

show(io::IO, ::IdentityOp{Pure, Generic, N}) where N = print(io, join(fill("Id", N), "⊗"))
show(io::IO, ::IdentityOp{Mixed, Generic, N}) where N = print(io, "Left(", join(fill("Id", N), "⊗"), ")")
show(io::IO, ::IdentityOp{Pure, Indexed}) = print(io, "Id")
show(io::IO, ::IdentityOp{Mixed, Indexed}) = print(io, "Left(Id)")

isless(::IdentityOp, ::IdentityOp) = false


################ Sums ##############

"""
    type SumOp{R, T, N} <: Op{R, T, N}

internal type for sum of operators
"""
struct SumOp{R, T, N} <: Op{R, T, N}
    subs::Vector{<:Op{R, T, N}}
    function SumOp(subs::Vector{<:Op{R, T, N}}) where {R, T, N}
        # a term of coefficient zero is left out: the zero operator has every parity, and
        # kept, it made `0C + dag(C)` a sum of fermionic and non fermionic operators
        s = filter(x -> scalarcoef(x) ≠ 0, reduce(vcat, sumsubs.(subs); init = Op{R, T, N}[]))
        if isempty(s)
            return 0 * IdentityOp{R, T, N}()
        elseif length(s) == 1
            return s[1]
        else
            return new{R, T, N}(s)
        end
    end
end

sumsubs(a::SumOp) = a.subs
sumsubs(a::Op) = [a]

(a::Op{R, T, N} + b::Op{R, T, N}) where {R, T, N} = SumOp([a, b])
(a::Op - b::Op) = a + (-b)

show(io::IO, a::SumOp) =
    paren(io, Base.operator_precedence(:+)) do io
        n = length(a.subs)
        compact = n > 6 && get(io, :compact, false)
        i = 1
        while (i <= n)
            s = sprint(show, a.subs[i]; context = io)
            if i == 1
                print(io, s)
            else
                if s[1] ≠ '-'
                    print(io, "+")
                elseif compact && i == 4
                    print(io, "-")
                end
                if compact && i == 4
                    print(io, "...")
                    i = n - 1
                    continue
                else
                    print(io, s)
                end
            end
            i = i + 1
        end
    end

isless(a::SumOp, b::SumOp) = isless(a.subs, b.subs)


################ Product by a number #############

"""
    type ScalarOp{R, T, N} <: Op{R, T, N}

internal type for product of number and operators
"""
struct ScalarOp{R, T, N} <: Op{R, T, N}
    coef::Number
    arg::Op{R, T, N}
    ScalarOp(coef::Number, arg::Op{R, T, N}) where {R, T, N} =
        if coef == 0
            new{R, T, N}(0, IdentityOp{R, T, N}())
        elseif coef == 1
            arg
        elseif arg isa SumOp
            SumOp(map(x -> coef * x, arg.subs))
        else
            # a signed zero taken out, -0.0 + 1.0im becoming 0.0 + 1.0im: `==` holds them
            # equal but `isless` and `hash` do not, and simplify, which sorts the terms before
            # merging the equal ones, left an interleaved pair unmerged
            c = coef * scalarcoef(arg)
            new{R, T, N}(c + zero(c), scalararg(arg))
        end
end

scalarcoef(a::ScalarOp) = a.coef
scalarcoef(::Op) = 1

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
    type ProdOp{R, T, N} <: Op{R, T, N}

internal type for product of operators
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

prodsubs(a::ProdOp) = a.subs
prodsubs(a::ScalarOp) = prodsubs(a.arg)
prodsubs(a::Op) = [a]

(a::Op{R, T, N} * b::Op{R, T, N}) where {R, T, N} =
    ProdOp([a, b])

show(io::IO, a::ProdOp) =
    paren(io, Base.operator_precedence(:*)) do io
        infix(io, a.subs, "*")
    end

isless(a::ProdOp, b::ProdOp) = isless(a.subs, b.subs)


############### Tensor products ############

"""
    type TensorOp{N} <: GenericOp{Pure, N}

internal type for tensor product of generic operators
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

tensorsubs(a::TensorOp) = a.subs
tensorsubs(a::GenericOp) = [a]

"""
    op1 ⊗ op2

tensor product for generic operators, alternative syntax: tensor(op1, op2)

# Examples

    controlled(op) = Proj("Up") ⊗ Id + Proj("Dn") ⊗ op
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
    type JW <: SimpleOp

type for operators transformed by the Jordan-Wigner transform.
for example C(5) is transformed into Multi_F(1,4)JW(C). 
JW operators anticommute with F
"""
struct JW <: SimpleOp
    arg::Operator{1}
end

isless(a::JW, b::JW) = isless(a.arg, b.arg)

"""
    type JW_F

the type of the F operator
"""
struct JW_F <: SimpleOp end

"""
    F

the Jordan Wigner F factor. Defined for all site types
"""
const F = JW_F()

show(io::IO, ::JW_F) =
    print(io, "F")

isless(::JW_F, ::JW_F) = false

"""
    Multi_F

a type for representing the F factors in the Jordan-Wigner transform.
for example C(5) is transformed into Multi_F(1,4)JW(C). 
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

an operator to project on the given state

# Examples

    Proj("Up")
    Proj([1, 0])
    Proj(1)   # project on the nth state (starting at 0)
"""
struct Proj <: SimpleOp
    state::Union{Int, String, Vector}
end

# the state is a number, a name or a vector, which have no order between them, and a vector
# may hold complex numbers, which have none either: what is compared is how they print, as
# for SetState, the global ordering of operators needing a total order, not a meaningful one
isless(a::Proj, b::Proj) = isless(repr(a.state), repr(b.state))


############ AtIndex ################

"""
    struct AtIndex{R, N} <: IndexedOp{R}

represent an indexed operator (like `X(1)` or `Swap(2, 4)`). An operator of several sites is
placed on distinct sites: `Swap(1, 1)` or `(X ⊗ Y)(1, 1)` is refused, write `X(1) * Y(1)`
for the product on one site.
"""
struct AtIndex{R, N} <: IndexedOp{R}
    op::GenericOp{R, N}
    index::NTuple{N, Int}
    # the matrix of an operator of several sites acts on distinct sites and says nothing of one
    # site taken twice: a gate then put one index into its tensor twice, and the adjoint of a
    # tensor product, developed factor by factor, relies on its sites being distinct
    function AtIndex(op::GenericOp{R, N}, index::NTuple{N, Int}) where {R, N}
        if !allunique(index)
            error("$op acts on $N sites and cannot be placed on $index, which repeats a site")
        end
        # the identity of the whole system, whatever sites it was placed on
        if scalararg(op) isa IdentityOp
            return scalarcoef(op) * IdentityOp{R, Indexed, 1}()
        end
        return scalarcoef(op) * new{R, N}(scalararg(op), index)
    end
end

(op::GenericOp{R, N})(index::Vararg{Int, N}) where {R, N} =
    AtIndex(op, index)

show(io::IO, ind::AtIndex) =
    if ind.op isa Operator
        show_func(io, ind.op.name, collect(ind.index))
    else
        show_func(io, repr(ind.op; context=:precedence=>500), collect(ind.index))
    end

isless(a::AtIndex, b::AtIndex) =
    isless((a.index, a.op), (b.index, b.op))


############## Mixers ###############

# Gate

"""
    Gate(op)

a generic operator acting as a gate on states in mixed representation. Useful for building noisy gates

# Examples

    G = 0.9 * Gate(Id) + 0.1 * Gate(X)
"""
struct Gate{N} <: GenericOp{Mixed, N}
    arg::GenericOp{Pure, N}
    Gate(arg::GenericOp{Pure, N}) where N =
        abs2(scalarcoef(arg)) *
        (scalararg(arg) isa IdentityOp ? IdentityOp{Mixed, Generic, N}() : new{N}(scalararg(arg)))
end

(a::IndexedOp{Mixed} * b::IndexedOp{Pure}) = a * Gate(b)
(a::IndexedOp{Pure} * b::IndexedOp{Mixed}) = Gate(a) * b

Gate(a::ProdOp{Pure, Indexed, 1}) = ProdOp(Gate.(a.subs))
Gate(ind::AtIndex{Pure}) = AtIndex(Gate(ind.op), ind.index)
Gate(::IdentityOp{Pure, Indexed, 1}) = IdentityOp{Mixed, Indexed, 1}()
# the same hoisting the inner constructor does, for an indexed operator: a gate built
# from c*A is rho -> (c A) rho (c A)' , that is abs2(c) times the gate built from A
Gate(a::ScalarOp{Pure, Indexed, 1}) = abs2(a.coef) * Gate(a.arg)

show(io::IO, a::Gate) =
    paren(io, 1000, 0) do io
        show_func(io, "Gate", a.arg)
    end

isless(a::Gate, b::Gate) = isless(a.arg, b.arg)

# Dissipator

"""
    Dissipator(op)

a Lindbladian dissipator based on `op` to be used in evolver for time evolution

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
    Evolver(op)

a Hamiltonian based on `op` to be used on mixed representation.
`op` should be of the form -im * hamiltonian
"""
struct Evolver <: IndexedOp{Mixed}
    arg::IndexedOp{Pure}
end

show(io::IO, a::Evolver) =
    paren(io, 1000, 0) do io
        show_func(io, "Evolver", a.arg)
    end

isless(a::Evolver, b::Evolver) = isless(a.arg, b.arg)

(a::IndexedOp{Mixed} + b::IndexedOp{Pure}) = a + Evolver(b)
(a::IndexedOp{Pure} + b::IndexedOp{Mixed}) = Evolver(a) + b

# Left

"""
    Left(op)

the superoperator acting on the left of the density matrix, ``\\rho \\mapsto A\\rho``.

`Gate`, `Dissipator` and `Evolver` are the usual ways of acting on a mixed representation,
and they all act on both sides at once: ``A\\rho A^\\dagger`` and ``-i[H, \\rho]``. This
one, and `Right`, are what is left for a term acting on a single side, which is not trace
preserving. They can be evolved with and applied as gates, but they are not observables:
`expect` has nothing to say about them.

# Examples

    apply(Left(X)(1), rho)          # rho -> X rho
    make_mpo(rho, Left(X)(1))
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
    Right(op)

the superoperator acting on the right of the density matrix, ``\\rho \\mapsto \\rho A†``.
See `Left`, of which this is the mirror.

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

# a gate built from a sum does not distribute: (A + B) rho (A + B)' has cross terms. On one
# site the sum is an operator of that site, whose gate is placed whole; otherwise the gate of K
# is Left(K) Right(K), which the placed factors of K give one by one
function Gate(a::SumOp{Pure, Indexed, 1})
    s = scalararg.(a.subs)
    if all(x -> x isa AtIndex && length(x.index) == 1, s) && allequal(x -> x.index, s)
        return Gate(sum(scalarcoef(x) * scalararg(x).op for x in a.subs))(only(first(s).index))
    end
    return ProdOp([sided(Left, a), sided(Right, a)])
end

# Left is linear and Right conjugates the coefficients, and both are multiplicative in the order
# of the factors, Right(A) Right(B) ρ being ρ B† A† = Right(A B) ρ. No fallback: a placed form
# with no method raises rather than being dropped
sided(S, a::AtIndex{Pure}) = AtIndex(S(a.op), a.index)
sided(S, a::SumOp{Pure, Indexed, 1}) = SumOp(map(x -> sided(S, x), a.subs))
sided(S, a::ProdOp{Pure, Indexed, 1}) = ProdOp(map(x -> sided(S, x), a.subs))
sided(::Type{Left}, a::ScalarOp{Pure, Indexed, 1}) = a.coef * sided(Left, a.arg)
sided(::Type{Right}, a::ScalarOp{Pure, Indexed, 1}) = conj(a.coef) * sided(Right, a.arg)
sided(_, ::IdentityOp{Pure, Indexed, 1}) = IdentityOp{Mixed, Indexed, 1}()

isless(a::Right, b::Right) =
    isless(a.arg, b.arg)


# SetState

"""
    SetState(state)

an operator to Set the local state to the one given, can only be used on mixed representations

# Examples

    SetState("Up")(3)
    SetState([0.3 0. ; 0. 0.7])(3)
"""
struct SetState <: GenericOp{Mixed, 1}
    state::Union{String, Vector, Matrix}
end

# the state is a name, a vector or a matrix. There is no order between those, and none
# at all between two matrices, so what is compared is how they print: the global ordering
# of operators needs a total order, not a meaningful one
isless(a::SetState, b::SetState) = isless(repr(a.state), repr(b.state))


############## Operator functions ###########

# powers

"""
    is_involution(op)

whether an operator is its own inverse, which its integer powers reduce to
"""
is_involution(::IdentityOp) = true
is_involution(::JW_F) = true
is_involution(a::Operator) = a.type == involution_op
is_involution(::Op) = false

# an exponent that makes a power a product: an integer of zero or above, whatever its type
is_natural(p::Number) = isreal(p) && isinteger(real(p)) && real(p) ≥ 0

"""
    type IntPowOp{R, N} <: GenericOp{R, N}
    type GenPowOp{R, N} <: GenericOp{R, N}

internal types for the powers of operators. An integer power of zero or above is a product,
`A^3` being `A * A * A`, kept whole only to be written that way: it has the adjoint, the
parity and the Jordan-Wigner strings of that product, and an involution is reduced at once,
`X^10000` being `Id`. Any other power, of an exponent that is not an integer, negative or
complex, is a function of the operator, the principal power taken through its logarithm, as
`exp` is: it cannot be split between factors, and of a coefficient of the operator only the
modulus comes out, which the logarithm takes apart exactly, the phase staying inside. Both merge,
`A^p * A^q` being `A^(p + q)`, the two powers being functions of the same logarithm.
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
            return is_natural(q) ? IntPowOp(a, Int(q)) : new{R, N}(a, q)
        end
        # a signed zero taken out, as in the coefficient of a ScalarOp
        return new{R, N}(arg, p + zero(p))
    end
end

# the power a merge of two powers of the same operator gives, either kind
power(a::GenericOp, p::Number) = is_natural(p) ? IntPowOp(a, Int(real(p))) : GenPowOp(a, p)

(a::GenericOp ^ p::Number) = power(a, p)

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

square root for generic operators
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
    type ExpOp{N} <: GenericOp{Pure, N}

an internal type to represent exponential of operators
"""
struct ExpOp{N} <: GenericOp{Pure, N}
    arg::GenericOp{Pure, N}
end

"""
    exp(::GenericOp{Pure})

exponential of a generic operator on pure states
"""
exp(a::GenericOp{Pure}) = ExpOp(a)

show(io::IO, a::ExpOp) =
    paren(io, 1000, 0) do io
        show_func(io, "exp", a.arg)
    end

isless(a::ExpOp, b::ExpOp) = isless(a.arg, b.arg)

# DagOp

"""
    type DagOp{N} <: GenericOp{Pure, N}

internal type to represent the `dag` operator
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

adjoint of a generic operator on pure states
"""
dag(a::GenericOp{Pure}) = DagOp(a)

show(io::IO, a::DagOp) =
    paren(io, 1000, 0) do io
        show_func(io, "dag", a.arg)
    end

isless(a::DagOp, b::DagOp) = isless(a.arg, b.arg)

# ModOp

"""
    type ModOp{N} <: GenericOp{Pure, N}

internal type to represent an operator taken modulo an integer, that is
``e^{2i\\pi A/m}``. Its eigenvalues are the `m`-th roots of unity of those of `A`, so a
conserved quantity written this way is conserved modulo `m` rather than as an integer.
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
    mod(op, m)

the operator ``e^{2i\\pi A/m}``, whose eigenvalues are the `m`-th roots of unity of those
of `A`. It is a genuine operator, which may be measured like any other, and it is what a
conserved quantity of ``\\mathbb{Z}_m`` is written with.

The modulus cannot be read back from the eigenvalues when it is 2, since ``\\pm 1`` is as
much a pair of integers as a pair of square roots of unity, and the two readings are
different conservations. That is why it is carried here rather than rediscovered.

# Examples

    mod(N, 3)                       # a Z3 charge
    Boson(6, conserve = mod(N, 3))
"""
mod(a::GenericOp{Pure}, m::Int) = ModOp(a, m)

"""
    parity(op)

the parity operator ``(-1)^A``, that is `mod(op, 2)`. Its expectation value is the usual
one, and as a conserved quantity it gives a charge of ``\\mathbb{Z}_2``.

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

# named

"""
    named_type(op)

the `OpType` a renamed operator keeps: its own when it has one, `fermionic_op` for a fermionic
expression, and nothing assumed otherwise.

A renamed operator of one site keeps its name through `simplify`, so its type is all that
tells it to take a Jordan-Wigner string, and `named(2C, "C2")` lost it. `controlled_type`
makes the opposite choice for a fermionic target: a controlled operator acts on several
sites, which no fermionic `Operator` can, and `simplify` replaces it by its expression.
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
  `fermionic_op` when it is fermionic and `plain_op` otherwise;
- a matrix, or anything given with its sites, takes the strongest type its matrix satisfies:
  `involution_op`, `selfadjoint_op` or `plain_op`, or `fermionic_op` when, on one site, it
  anticommutes with `F`;
- a function without sites is `plain_op`.

The type is checked against the matrix each time the operator is placed on a site.

Renaming also keeps two conserved quantities apart: two sites declaring the same operator
name the same charge, and another name makes another charge.

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


################## isfermionic #################

"""
    isfermionic(::SimpleOp)

whether an operator of one site on pure states is odd under the fermion parity. A sum mixing
fermionic and non fermionic operators, and an exponential, a `mod` or a non integer power of
a fermionic operator, have no such parity and raise an error. For an operator of several
sites or a superoperator, see `has_fermionic`.
"""
isfermionic(a::SimpleOp) = false
isfermionic(a::Operator{1}) = a.type == fermionic_op
isfermionic(a::ScalarOp{Pure}) = isfermionic(a.arg)
isfermionic(a::DagOp) = isfermionic(a.arg)
isfermionic(a::ProdOp{Pure, Generic, 1}) = isodd(count(isfermionic, a.subs))
isfermionic(a::ExpOp) =
    if isfermionic(a.arg)
        error("cannot take the exponential of the fermionic operator $(a.arg)")
    else
        false
    end

isfermionic(a::ModOp) =
    if isfermionic(a.arg)
        error("cannot compute $a, which exponentiates the fermionic operator $(a.arg)")
    else
        false
    end
function isfermionic(a::SumOp{Pure, Generic, 1})
    n = length(a.subs)
    nf = count(isfermionic, a.subs)
    if n == nf
        return true
    elseif nf == 0
        return false
    else
        error("cannot sum fermionic and non fermionic operators ($a)")
    end
end

isfermionic(a::IntPowOp) = isfermionic(a.arg) && isodd(a.expo)

isfermionic(a::GenPowOp) =
    if isfermionic(a.arg)
        error("cannot take the non integer power $(a.expo) of the fermionic operator $(a.arg)")
    else
        false
    end

"""
    has_fermionic(::Op)

whether an operator still holds a factor whose Jordan-Wigner string `simplify` has not
inserted: an odd operator of one site, wherever it sits. A tensor product, an operator of
several sites defined by an expression and the argument of a superoperator are all looked
into, since their tensor carries the strings between consecutive sites only: placed on sites
apart, it misses those of the sites in between, which `(C ⊗ dag(C))(1, 3) = C(1) * dag(C)(3)`
asks for.

The structure is read from the fields, the way `==` is, so that no wrapper can be forgotten.
Once `simplify` has run, the factors are `JW` transforms, which are not fermionic, and the
answer is false unless a factor was left whole, as an exponential of several sites is. A
factor of no definite parity, as `C + N`, holds an odd part and answers true, where it raised.
"""
has_fermionic(a::GenericOp{Pure, 1}) = fermion_parity(a, false) ≠ 0

function has_fermionic(a::Op)
    for f in fieldnames(typeof(a))
        x = getfield(a, f)
        if x isa Op && has_fermionic(x)
            return true
        elseif x isa Vector && any(y -> y isa Op && has_fermionic(y), x)
            return true
        end
    end
    return false
end


"""
    jw_parity(a)

how a factor of a one site product behaves when the `F` of that site crosses it: `0` when it
commutes with `F`, `1` when it anticommutes, `nothing` when it does neither.

A fermionic operator, and the Jordan-Wigner transform `simplify` makes of it, is odd, and any
other operator is taken to be even, the convention the strings themselves rest on. A composite
factor has the parity its pieces give it. Sums have to be read as well, since `simplify`
gathers the terms of one site into a single factor: `(C + dag(C))(1)` was taken to be even,
and `C(3) * (C + dag(C))(1)` came out with the wrong sign.
"""
jw_parity(a) = fermion_parity(a, true)

# the reading shared by jw_parity and has_fermionic. `strung` is whether the Jordan-Wigner
# transforms simplify inserts count: odd for jw_parity, which carries F across them, even for
# has_fermionic, which asks for a string not yet inserted. A projector on a vector may mix the
# two parities, which matters to the F crossing it and not to a string, which it never takes
fermion_parity(::Op, ::Bool) = 0
fermion_parity(::JW, strung::Bool) = strung ? 1 : 0
fermion_parity(a::Operator, ::Bool) = a.type == fermionic_op ? 1 : 0
fermion_parity(a::Union{ScalarOp, DagOp}, strung::Bool) = fermion_parity(a.arg, strung)
# a projector on a basis state, given by its index or by a name, is even: the named states of
# the fermionic sites are all basis states, which a site defined outside the package is taken
# to follow
fermion_parity(a::Proj, strung::Bool) = strung && a.state isa Vector ? nothing : 0

function fermion_parity(a::ProdOp, strung::Bool)
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

ranking(::IntPowOp) = 30
ranking(::GenPowOp) = 34
ranking(::ExpOp) = 31
ranking(::DagOp) = 32
ranking(::ModOp) = 33

"""
    obs_name(op)

the name a measurement is given when the caller did not choose one: how the operator
prints, compactly. Only `Operator` carries a name of its own; everything else a
measurement may be asked for is a composition, whose printed form is its only description.
`SimpleOp` is abstract, so reading a `name` field would work for a bare `X` and fail for
`X * Y`, `2X` or `X + Y`.
"""
obs_name(op) = sprint(print, op; context = :compact => true)
