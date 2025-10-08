export PM, Pure, Mixed, GI, Generic, Indexed, GenericOp, IndexedOp, SimpleOp
export OpType, plain_op, fermionic_op, selfadjoint_op, involution_op
export Operator, Identity, Id, JW, JW_F, F, Proj, Gate, Dissipator, Evolver, Left, Right, Multi_F
export dag, tensor, ⊗, isfermionic

############# Types ################

"""
    abstract type PM

PM is a supertype for Pure and Mixed (i.e. Pure() and Mixed() are of type PM)
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

abstract type GI end

struct Generic <: GI end
struct Indexed <: GI end

abstract type Op{R <: PM, T <: GI, N} end

GenericOp{R, N} = Op{R, Generic, N}
IndexedOp{R} = Op{R, Indexed, 1}
SimpleOp = GenericOp{Pure, 1}

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
                print(io, real(a))
            elseif real(a) == 0
                print_coef(io, imag(a))
                print(io, "im*")
            else
                print(io, "(", a, ")")
            end
        else
            print(io, a)
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

@enum OpType plain_op fermionic_op selfadjoint_op involution_op

"""
    type Operator{N} <: GenericOp{N}

the type of base operators (like `X`, `Swap`, `C` ...), `N` is an Int.

# Example

    Operator("X")                      a base operator whose value is predefined by the sites
    Operator("Z", [1 0 ; 0 -1])
    Operator{2}("Swap", [ 1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1])
    Operator("CX", controlled(X))
    Operator("Sx", (Sp + Sm) / 2)
    Operator("C", [0 1 ; 0 0], fermionic_op)   a fermionic operator
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

struct Identity <: SimpleOp end

"""
    Id

the identity operator defined for all site types
"""
const Id = Identity()

struct MakeIdentity{R, T, N} <: Op{R, T, N}
    MakeIdentity{Pure, Generic, N}() where N = TensorOp{N}(fill(Id, N))
    MakeIdentity{Mixed, Generic, N}() where N = Left(TensorOp{N}(fill(Id, N)))
    MakeIdentity{Pure, Indexed, 1}() = Id(1)
    MakeIdentity{Mixed, Indexed, 1}() = Left(Id)(1)
    MakeIdentity(::Op{R, T, N}) where {R, T, N} = MakeIdentity{R, T, N}()
end

show(io::IO, ::Identity) =
    print(io, "Id")

isless(::Identity, ::Identity) = false


################ Sums ##############

struct SumOp{R, T, N} <: Op{R, T, N}
    subs::Vector{<:Op{R, T, N}}
    function SumOp(subs::Vector{<:Op{R, T, N}}) where {R, T, N}
        if isempty(subs)
            return 0 * MakeIdentity{R, T, N}()
        end
        s = reduce(vcat, sumsubs.(subs))
        if length(s) == 1
            s[1]
        else
            new{R, T, N}(s)
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
(a::SumOp == b::SumOp) = a.subs == b.subs


################ Product by a number #############

struct ScalarOp{R, T, N} <: Op{R, T, N}
    coef::Number
    arg::Op{R, T, N}
    ScalarOp(coef::Number, arg::Op{R, T, N}) where {R, T, N} =
        if coef == 0
            new{R, T, N}(0, MakeIdentity{R, T, N}())
        elseif coef == 1
            arg
        elseif arg isa SumOp
            SumOp(map(x -> coef * x, arg.subs))
        else
            new{R, T, N}(coef * scalarcoef(arg), scalararg(arg))
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

isless(a::ScalarOp, b::ScalarOp) = isless((a.coef, a.arg), (b.coef, b.arg))


############### Products ###############

struct ProdOp{R, T, N} <: Op{R, T, N}
    subs::Vector{<:Op{R, T, N}}
    function ProdOp(subs::Vector{<:Op{R, T, N}}) where {R, T, N}
        c = prod(scalarcoef.(subs))
        if isempty(subs)
            return c * MakeIdentity{R, T, N}()
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
        join(io, a.subs, "*")
    end

isless(a::ProdOp, b::ProdOp) = isless(a.subs, b.subs)
(a::ProdOp == b::ProdOp) = a.subs == b.subs


############### Tensor products ############

struct TensorOp{N} <: GenericOp{Pure, N}
    subs::Vector{<:GenericOp{Pure}}
    TensorOp{N}(subs::Vector{<:GenericOp{Pure}}) where N =
        if length(subs) == 1
            subs[1]
        else
            c = prod(scalarcoef.(subs))
            s = scalararg.(subs)
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
    Rxy(t) = exp(im * t * (X ⊗ X + Y ⊗ Y) / 4)

"""
(a::GenericOp{Pure, N} ⊗ b::GenericOp{Pure, M}) where {N, M} =
    TensorOp{N + M}([tensorsubs(a) ; tensorsubs(b)])
tensor(a::GenericOp, b::GenericOp) = a ⊗ b
tensor(a::GenericOp, b::GenericOp, c::GenericOp, d::GenericOp...) = tensor(a ⊗ b, c, d...)

show(io::IO, a::TensorOp) =
    paren(io, Base.operator_precedence(:⊗)) do io
        join(io, a.subs, "⊗")
    end

isless(a::TensorOp, b::TensorOp) = isless(a.subs, b.subs)
(a::TensorOp == b::TensorOp) = a.subs == b.subs


############# Jordan_Wigner transformation ##############

struct JW <: SimpleOp
    arg::Operator{1}
end

isless(a::JW, b::JW) = isless(a.arg, b.arg)

struct JW_F <: SimpleOp end

"""
    The Jordan Wigner F factor. Defined for all site types
"""
const F = JW_F()

show(io::IO, ::JW_F) =
    print(io, "F")

isless(::JW_F, ::JW_F) = false

struct Multi_F{R} <: IndexedOp{R}
    start::Int
    stop::Int
    left::Bool
    right::Bool
    Multi_F{R}(start::Int, stop::Int, left::Bool, right::Bool) where R =
        if start > stop || (R == Mixed && !left && !right)
            MakeIdentity{R, Indexed, 1}()
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
"""
struct Proj <: SimpleOp
    state::Union{String, Vector}
end

isless(a::Proj, b::Proj) = isless(a.state, b.state)


############ AtIndex ################

"""
    struct AtIndex{R, N} <: IndexedOp{R}

represent a base indexed operator (like `X(1)` or `Swap(2, 4)`)
"""
struct AtIndex{R, N} <: IndexedOp{R}
    op::GenericOp{R, N}
    index::NTuple{N, Int}
    AtIndex(op::GenericOp{R, N}, index::NTuple{N, Int}) where {R, N} =
        scalarcoef(op) * new{R, N}(scalararg(op), index)
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
(a::AtIndex == b::AtIndex) = a.op == b.op && a.index == b.index


############## Mixers ###############

# Gate

"""
    Gate(op)

a generic operator acting as a gate on states in mixed representation. Usefull for building noisy gates

# Examples

    G = 0.9 * Gate(Id) + 0.1 Gate(X)
"""
struct Gate{N} <: GenericOp{Mixed, N}
    arg::GenericOp{Pure, N}
    Gate(arg::GenericOp{Pure, N}) where N =
        abs2(scalarcoef(arg)) * new{N}(scalararg(arg))
end

(a::IndexedOp{Mixed} * b::IndexedOp{Pure}) = a * Gate(b)
(a::IndexedOp{Pure} * b::IndexedOp{Mixed}) = Gate(a) * b

Gate(a::ProdOp{Pure, Indexed, 1}) = ProdOp(Gate.(a.subs))
Gate(ind::AtIndex{Pure}) = AtIndex(Gate(ind.op), ind.index)

show(io::IO, a::Gate) =
    paren(io, 1000, 0) do io
        show_func(io, "Gate", a.arg)
    end

isless(a::Gate, b::Gate) = isless(a.arg, b.arg)

# Dissipator

"""
    Dissipator(op)

a Lindbladian dissipator based on `op` to be used in evolver for time evolution
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

(a::IndexedOp{Mixed} + b::IndexedOp{Pure}) = a + Evolver(b)
(a::IndexedOp{Pure} + b::IndexedOp{Mixed}) = Evolver(a) + b

# Left

struct Left{N} <: GenericOp{Mixed, N}
    arg::GenericOp{Pure, N}
    Left(arg::GenericOp{Pure, N}) where N =
        scalarcoef(arg) * new{N}(scalararg(arg))
end

show(io::IO, a::Left) =
    paren(io, 1000, 0) do io
        show_func(io, "Left", a.arg)
    end

isless(a::Left, b::Left) =
    isless(a.arg, b.arg)

# Right

struct Right{N} <: GenericOp{Mixed, N}
    arg::GenericOp{Pure, N}
    Right(arg::GenericOp{Pure, N}) where N =
        conj(scalarcoef(arg)) * new{N}(scalararg(arg))
end

show(io::IO, a::Right) =
    paren(io, 1000, 0) do io
        show_func(io, "Right", a.arg)
    end

isless(a::Right, b::Right) =
    isless(a.arg, b.arg)


############## Operator functions ###########

# PowOp

struct PowOp{R, N} <: GenericOp{R, N}
    arg::GenericOp{R, N}
    expo::Number
    PowOp(arg::GenericOp{R, N}, expo::Number) where {R, N} =
        if expo == 0
            MakeIdentity{R, Generic, N}()
        elseif expo == 1
            arg
        elseif expo < 0
            error("cannot take negative powers of operators")
        else
            c = scalarcoef(arg)
            a = scalararg(arg)
            if c >= 0 || isinteger(expo)
                c^expo * new{R, N}(a, expo)
            else
                (-c)^expo * new{R, N}(-a, expo)
            end
        end
end

(a::GenericOp ^ b::Number) = PowOp(a, b)

"""
    sqrt(::GenericOp)

square root for generic tensors
"""
sqrt(a::GenericOp) = a ^ 0.5 

show(io::IO, a::PowOp) =
    paren(io, Base.operator_precedence(:^)) do io
        print(io, a.arg, "^", a.expo)
    end

isless(a::PowOp, b::PowOp) =
    isless((a.arg, a.expo), (b.arg, b.expo))

# ExpOp

struct ExpOp{N} <: GenericOp{Pure, N}
    arg::GenericOp{Pure, N}
end

"""
    exp(::GenericOp)

exponential for generic operators
"""
exp(a::GenericOp{Pure}) = ExpOp(a)

show(io::IO, a::ExpOp) =
    paren(io, 1000, 0) do io
        show_func(io, "exp", a.arg)
    end

isless(a::ExpOp, b::ExpOp) = isless(a.arg, b.arg)

# DagOp

struct DagOp{N} <: GenericOp{Pure, N}
    arg::GenericOp{Pure, N}
    DagOp(arg::DagOp) = arg
    DagOp(arg::ScalarOp{Pure, Generic}) =
        conj(a.coef) * DagOp(arg.arg)
    DagOp(arg::GenericOp{Pure, N}) where N =
        new{N}(arg)
end

"""
    dag(::GenericOp)

adjoint for generic operators
"""
dag(a::GenericOp{Pure}) = DagOp(a)

show(io::IO, a::DagOp) =
    paren(io, 1000, 0) do io
        show_func(io, "dag", a.arg)
    end

isless(a::DagOp, b::DagOp) = isless(a.arg, b.arg)


################## isfermionic #################

isfermionic(a::SimpleOp) = false
isfermionic(a::Operator{1}) = a.type == fermionic_op
isfermionic(a::Union{ScalarOp{Pure}, DagOp{Pure}}) = isfermionic(a.arg)
isfermionic(a::ProdOp{Pure, Generic, 1}) = isodd(count(isfermionic, a.subs))
isfermionic(a::Union{ExpOp}) =
    if isfermionic(a.arg)
        error("cannot take exp of fermionic operators ($a)")
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

isfermionic(a::PowOp) =
    if !isfermionic(a.arg)
        false
    elseif isinteger(a.expo)
        isodd(a.expo)
    else
        error("cannot determine fermionic nature of $a")
    end


################## Global Ordering ###############

ranking(a) = error("ranking not defined for ($a)")

isless(a::Op, b::Op) = isless((ranking(a), a), (ranking(b), b))

ranking(::Identity) = 1
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

ranking(::ScalarOp) = 20
ranking(::ProdOp) = 21
ranking(::SumOp) = 22
ranking(::TensorOp) = 23

ranking(::PowOp) = 30
ranking(::ExpOp) = 31
ranking(::DagOp) = 32
