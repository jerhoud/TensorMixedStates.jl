export Repr, Pure, Mixed, AnyOp, GenericOp, IndexedOp, SimpleOp
export OpType, plain_op, fermionic_op, selfadjoint_op, involution_op
export Operator, Proj, Gate, Dissipator, Evolver, LeftRight
export tensor, ⊗

############# Types ################

"""
    abstract type Repr

Repr is a supertype for Pure and Mixed (i.e. Pure() and Mixed() are of type Repr)
"""
abstract type Repr end

"""
    type Pure <: Repr
    Pure()

correspond to pure quantum representation
"""
struct Pure <: Repr end


"""
    type Mixed <: Repr
    Mixed()

correspond to mixed quantum representation
"""
struct Mixed <: Repr end

abstract type AnyOp end

abstract type GenericOp{R <: Repr, N} <: AnyOp end
abstract type IndexedOp{R <: Repr} <: AnyOp end

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
    if out_prec <= ext_prec
        print(io, "(")
    end
    f(IOContext(io, :precedence => in_prec))
    if out_prec <= ext_prec
        print(io, ")")
    end
end


############### Operator ###############

@enum OpType plain_op fermionic_op selfadjoint_op involution_op

"""
    type Operator{R, N} <: GenericOp{R, N}

the type of base operators (like `X`, `Swap`, `C` ...), `R` is Pure or Mixed and `N` is an Int.

# Example

    Operator("X")                      a base operator whose value is predefined by the sites
    Operator("Z", [1 0 ; 0 -1])
    Operator{Pure, 2}("Swap", [ 1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1])
    Operator("CX", controlled(X))
    Operator("Sx", (Sp + Sm) / 2)
    Operator("C", [0 1 ; 0 0], fermionic)   a fermionic operator
"""
struct Operator{R, N} <: GenericOp{R, N}
    name::String
    expr::Union{Nothing, Matrix, Function, GenericOp{R, N}}
    type::OpType
    Operator{R, N}(name, expr = nothing, type = plain_op) where {R, N} =
        new(name, expr, type)
end

Operator(args...) = Operator{Pure, 1}(args...)
Operator(name, expr::GenericOp{R, N}, type=plain_op) where {R, N} = Operator{R, N}(name, expr, plain_op)

show(io::IO, op::Operator) =
    print(io, op.name)

isless(a::Operator, b::Operator) = isless(a.name, b.name)

############ proj ################

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

############ Indexed ################

"""
    struct Indexed{R, N} <: IndexedOp{R}

represent a base indexed operator (like `X(1)` or `Swap(2, 4)`)
"""
struct Indexed{R, N} <: IndexedOp{R}
    op::GenericOp{R, N}
    index::NTuple{N, Int}
end

(op::GenericOp{R, N})(index::Vararg{Int, N}) where {R, N} =
    Indexed{R, N}(op, index)

show(io::IO, ind::Indexed) =
    if ind.op isa Operator
        show_func(io, ind.op.name, collect(ind.index))
    else
        show_func(io, repr(ind.op; context=:precedence=>500), collect(ind.index))
    end

isless(a::Indexed, b::Indexed) =
    isless((a.index, fermionic(a.op), a.op), (b.index, fermionic(b.op), b.op))

# Why sort on fermionic ?

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
end

Gate(ind::Indexed{Pure}) = Indexed(Gate(ind.op), ind.index)

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

# LeftRight

struct LeftRight <: IndexedOp{Mixed}
    left::IndexedOp{Pure}
    right::IndexedOp{Pure}
end

show(io::IO, a::LeftRight) =
    paren(io, 1000, 0) do io
        show_func(io, "LeftRight", [a.left, a.right])
    end

isless(a::LeftRight, b::LeftRight) =
    isless((a.left, a.right), (b.left, b.right))

############### Operator Products ###############

# ProdGenOp    product of generic operators
struct ProdGenOp{R, N} <: GenericOp{R, N}
    coef::Number
    subs::Vector{<:GenericOp{R, N}}
end

# ProdIndOp    product of indexed operators
struct ProdIndOp{R} <: IndexedOp{R}
    coef::Number
    subs::Vector{<:IndexedOp{R}}
end

AnyProd = Union{ProdGenOp, ProdIndOp}

prodsubs(a::AnyProd) = a.subs
prodsubs(a::AnyOp) = [a]

prodcoef(a::AnyProd) = a.coef
prodcoef(a::AnyOp) = 1

(a::GenericOp{R, N} * b::GenericOp{R, N}) where {R, N} =
    ProdGenOp{R, N}(prodcoef(a) * prodcoef(b), [prodsubs(a) ; prodsubs(b)])
(a::Number * b::GenericOp{R, N}) where {R, N} = 
    ProdGenOp{R, N}(a * prodcoef(b), prodsubs(b))
(a::IndexedOp{R} * b::IndexedOp{R}) where R =
    ProdIndOp{R, N}(prodcoef(a) * prodcoef(b), [prodsubs(a) ; prodsubs(b)])
(a::Number * b::IndexedOp{R}) where R = 
    ProdIndOp{R}(a * prodcoef(b), prodsubs(b))

(a::AnyOp * b::Number) = b * a
(a::AnyOp / b::Number) = inv(b) * a
-(a::AnyOp) = -1 * a

(a::IndexedOp{Mixed} * b::IndexedOp{Pure}) = a * Gate(b)
(a::IndexedOp{Pure} * b::IndexedOp{Mixed}) = Gate(a) * b

Gate(a::ProdIndOp{Pure}) = ProdIndOp{Mixed}(abs2(a.coef), Gate.(a.subs))

show(io::IO, a::AnyProd) =
    paren(io, Base.operator_precedence(:*)) do io
        print_coef(io, a.coef)
        join(io, a.subs, "*")
    end

isless(a::ProdGenOp, b::ProdGenOp) = isless(a.subs, b.subs)
isless(a::ProdIndOp, b::ProdIndOp) = isless(a.subs, b.subs)


################ Operator Sums ##############

# SumGenOp  sum of generic operators
struct SumGenOp{R, N} <: GenericOp{R, N}
    subs::Vector{<:GenericOp{R, N}}
end

# SumIndOp  sum of indexed operators
struct SumIndOp{R} <: IndexedOp{R}
    subs::Vector{<:IndexedOp{R}}
end

AnySum = Union{SumGenOp, SumIndOp}

sumsubs(a::AnySum) = a.subs
sumsubs(a::AnyOp) = [a]

(a::GenericOp{R, N} + b::GenericOp{R, N}) where {R, N} =
    SumGenOp{R, N}([sumsubs(a) ; sumsubs(b)])
(a::IndexedOp{R} + b::IndexedOp{R}) where R =
    SumIndOp{R}([sumsubs(a) ; sumsubs(b)])
(a::AnyOp - b::AnyOp) = a + (-b)

(a::IndexedOp{Mixed} + b::IndexedOp{Pure}) = a + Evolver(b)
(a::IndexedOp{Pure} + b::IndexedOp{Mixed}) = Evolver(a) + b

show(io::IO, a::AnyOp) =
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

isless(a::SumGenOp, b::SumGenOp) = isless(a.subs, b.subs)
isless(a::SumIndOp, b::SumIndOp) = isless(a.subs, b.subs)


############### Tensor products ############

struct TensorOp{R, N} <: GenericOp{R, N}
    subs::Vector{GenericOp{R}}
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
(a::GenericOp{R, N} ⊗ b::GenericOp{R, M}) where {R, N, M} =
    TensorOp{R, N + M}([tensorsubs(a) ; tensorsubs(b)])
tensor(a::GenericOp, b::GenericOp) = a ⊗ b
tensor(a::GenericOp, b::GenericOp, c::GenericOp, d::GenericOp...) = tensor(a ⊗ b, c, d...)

show(io::IO, a::TensorOp) =
    paren(io, Base.operator_precedence(:⊗)) do io
        join(io, a.subs, "⊗")
    end

isless(a::TensorOp, b::TensorOp) = isless(a.subs, b.subs)


############## Operator functions ###########

# PowOp

struct PowOp{R, N} <: GenericOp{R, N}
    arg::GenericOp{R, N}
    expo::Number
end

(a::GenericOp ^ b::Number) = PowOp(a, b)

show(io::IO, a::PowOp) =
    paren(io, Base.operator_precedence(:^)) do io
        print(io, a.arg, "^", a.expo)
    end

isless(a::PowOp, b::PowOp) =
    isless((a.arg, a.expo), (b.arg, b.expo))

# ExpOp

struct ExpOp{R, N} <: GenericOp{R, N}
    arg::GenericOp{R, N}
end

"""
    exp(::GenericOp)

exponential for generic tensors
"""
exp(a::GenericOp) = ExpOp(a)

show(io::IO, a::ExpOp) =
    paren(io, 1000, 0) do io
        show_func(io, "exp", a.arg)
    end

isless(a::ExpOp, b::ExpOp) = isless(a.arg, b.arg)

# SqrtOp

"""
    sqrt(::GenericOp)

square root for generic tensors
"""
struct SqrtOp{R, N} <: GenericOp{R, N}
    arg::GenericOp{R, N}
end

sqrt(a::GenericOp) = SqrtOp(a)

show(io::IO, a::SqrtOp) =
    paren(io, 1000, 0) do io
        show_func(io, "sqrt", a.arg)
    end

isless(a::SqrtOp, b::SqrtOp) = isless(a.arg, b.arg)

# DagOp

struct DagOp{R, N} <: GenericOp{R, N}
    arg::GenericOp{R, N}
end

"""
    dag(::GenericOp)

adjoint for generic operators
"""
dag(a::GenericOp) = DagOp(a)
dag(a::DagOp) = a.arg

show(io::IO, a::DagOp) =
    paren(io, 1000, 0) do io
        show_func(io, "dag", a.arg)
    end

isless(a::DagOp, b::DagOp) = isless(a.arg, b.arg)


################## Global Ordering ###############

ranking(a) = error("ranking not defined for ($a)")

# GenericOp
isless(a::GenericOp, b::GenericOp) = isless((ranking(a), a), (ranking(b), b))
ranking(::Operator) = 1
ranking(::Proj) = 2
ranking(::Gate) = 3
ranking(::ProdGenOp) = 4
ranking(::SumGenOp) = 5
ranking(::TensorOp) = 6
ranking(::PowOp) = 7
ranking(::ExpOp) = 8
ranking(::SqrtOp) = 9
ranking(::DagOp) = 10

# IndexedOp
isless(a::IndexedOp, b::IndexedOp) = isless((ranking(a), a), (ranking(b), b))
ranking(::Indexed) = 1
ranking(::LeftRight) = 2
ranking(::ProdIndOp) = 3
ranking(::SumIndOp) = 4

