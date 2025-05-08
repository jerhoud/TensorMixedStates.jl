export ExprOp, SimpleOp, ProdOp, SumOp, TensorOp, ExpOp, Operator, Indexed, IndexOp, ⊗, dag, DagOp, Proj
export Gate, Dissipator, Evolver, Left, Right, Pure, Mixed, PM


############# Types ################

"""
    type Pure
    Pure()

correspond to pure quantum representation
"""
struct Pure end


"""
    type Mixed
    Mixed()

correspond to mixed quantum representation
"""
struct Mixed end

PM = Union{Pure, Mixed}

"""
    abstract type ExprOp{T, N}

the abstract type of operators, `T` is `Pure` or `Mixed` and `N` is an `Int`
for generic operators (the number of sites on which the operator may apply)
and `IndexOp` for indexed operators

    SimpleOp = ExprOp{Pure, 1}
    ExprIndexed{T} = ExprOp{T, IndexOp}
"""
abstract type ExprOp{T, N} end
SimpleOp = ExprOp{Pure, 1}

isless(a::ExprOp, b::ExprOp) = isless((ranking(a), a), (ranking(b), b))

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

"""
    type Operator{T, N} <: ExprOp{T, N}

the type of base operators (like `X`, `Swap`, `C` ...), `T` is `Pure` or `Mixed` and `N` is an Int.

# Example

    Operator("X")                      a base operator whose value is predefined by the sites
    Operator("Z", [1 0 ; 0 -1])
    Operator{Pure, 2}("Swap", [ 1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1])
    Operator("CX", controlled(X))
    Operator("Sx", (Sp + Sm) / 2)
    Operator("C", [0 1 ; 0 0], true)   a fermionic operator
"""
struct Operator{T, N} <: ExprOp{T, N}
    name::String
    expr::Union{Nothing, Matrix, Function, ExprOp{T, N}}
    fermionic::Bool
    Operator{T, N}(name, expr = nothing, fermionic = false) where {T, N} =
    new(name, expr, fermionic)
end

Operator(args...) = Operator{Pure, 1}(args...)
Operator(name, expr::ExprOp{T, N}, fermionic=false) where {T, N} = Operator{T, N}(name, expr, fermionic)

show(io::IO, op::Operator) =
    print(io, op.name)

ranking(::Operator) = 1

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

ranking(::Proj) = 2

isless(a::Proj, b::Proj) = isless(a.state, b.state)

############ Indexed ################

struct IndexOp end

ExprIndexed{T} = ExprOp{T, IndexOp}

"""
    type Indexed{T, N} <: ExprIndexed{T}

represent a base indexed operator (like `X(1)` or `Swap(2, 4)`)
"""
struct Indexed{T, N} <: ExprIndexed{T}
    op::ExprOp{T, N}
    index::NTuple{N, Int}
end

(op::ExprOp{T, N})(index::Vararg{Int, N}) where {T, N} =
    Indexed{T, N}(op, index)

show(io::IO, ind::Indexed) =
    if ind.op isa Operator
        show_func(io, ind.op.name, collect(ind.index))
    else
        show_func(io, repr(ind.op; context=:precedence=>500), collect(ind.index))
    end

ranking(::Indexed) = 3

isless(a::Indexed, b::Indexed) = isless((a.index, fermionic(a.op), a.op), (b.index, fermionic(b.op), b.op))

apply_expr(f, a::Indexed) = Indexed(f(a.op), a.index)

############## Mixers ###############

# Gate

"""
    Gate(op)

a generic operator acting as a gate on states in mixed representation. Usefull for building noisy gates

# Examples

    G = 0.9 * Gate(Id) + 0.1 Gate(X)
"""
struct Gate{N} <: ExprOp{Mixed, N}
    arg::ExprOp{Pure, N}
end

show(io::IO, a::Gate) =
    paren(io, 1000, 0) do io
        show_func(io, "Gate", a.arg)
    end

ranking(::Gate) = 10

isless(a::Gate, b::Gate) = isless(a.arg, b.arg)

apply_expr(f, a::Gate) = Gate(f(a.arg))

# Dissipator

"""
    Dissipator(op)

a Lindbladian dissipator based on `op` to be used in evolver for time evolution
"""
struct Dissipator{N} <: ExprOp{Mixed, N}
    arg::ExprOp{Pure, N}
end

Dissipator(::ExprIndexed{Pure}) = error("Dissipator cannot be applied to indexed expressions")

show(io::IO, a::Dissipator) =
    paren(io, 1000, 0) do io
        show_func(io, "Dissipator", a.arg)
    end

ranking(::Dissipator) = 11

isless(a::Dissipator, b::Dissipator) = isless(a.arg, b.arg)

apply_expr(f, a::Dissipator) = Dissipator(f(a.arg))

# Evolver

struct Evolver <: ExprIndexed{Mixed}
    arg::ExprIndexed{Pure}
end

show(io::IO, a::Evolver) =
    paren(io, 1000, 0) do io
        show_func(io, "Evolver", a.arg)
    end

ranking(::Evolver) = 12

isless(a::Evolver, b::Evolver) = isless(a.arg, b.arg)

apply_expr(f, a::Evolver) = Evolver(f(a.arg))

# Left

struct Left{N} <: ExprOp{Mixed, N}
    arg::ExprOp{Pure, N}
end

show(io::IO, a::Left) =
    paren(io, 1000, 0) do io
        show_func(io, "Left", a.arg)
    end
    
ranking(::Left) = 13

isless(a::Left, b::Left) = isless(a.arg, b.arg)

apply_expr(f, a::Left) = Left(f(a.arg))

# Right

struct Right{N} <: ExprOp{Mixed, N}
    arg::ExprOp{Pure, N}
end

show(io::IO, a::Right) =
    paren(io, 1000, 0) do io
        show_func(io, "Right", a.arg)
    end

ranking(::Right) = 14

isless(a::Right, b::Right) = isless(a.arg, b.arg)
    
apply_expr(f, a::Right) = Right(f(a.arg))
    
############### Operator Products ###############

struct ProdOp{T, N} <: ExprOp{T, N}
    coef::Number
    subs::Vector{<:ExprOp{T, N}}
    ProdOp{T, N}(c::Number, s::Vector{<:ExprOp{T, N}}) where {T, N} =
        if isempty(s)
            error("bug: empty operator product")
        elseif c == 1 && length(s) == 1
            return s[1]
        else
            return new{T, N}(c, s)
        end
end

prodsubs(a::ProdOp) = a.subs
prodsubs(a::ExprOp) = [a]
prodcoef(a::ProdOp) = a.coef
prodcoef(a::ExprOp) = 1

(a::ExprOp{T, N} * b::ExprOp{T, N}) where {T, N} =
    ProdOp{T, N}(prodcoef(a) * prodcoef(b), [prodsubs(a) ; prodsubs(b)])
(a::ExprIndexed{Mixed} * b::ExprIndexed{Pure}) = a * Gate(b)
(a::ExprIndexed{Pure} * b::ExprIndexed{Mixed}) = Gate(a) * b
(a::Number * b::ExprOp{T, N}) where {T, N} = 
    ProdOp{T, N}(a * prodcoef(b), prodsubs(b))
(a::ExprOp * b::Number) = b * a
(a::ExprOp / b::Number) = inv(b) * a
-(a::ExprOp) = -1 * a

show(io::IO, a::ProdOp) =
    paren(io, Base.operator_precedence(:*)) do io
        print_coef(io, a.coef)
        join(io, a.subs, "*")
    end

ranking(::ProdOp) = 20

isless(a::ProdOp, b::ProdOp) = isless(a.subs, b.subs)

apply_expr(f, a::ProdOp) = a.coef * reduce(*, map(f, a.subs))

################ Operator Sums ##############

struct SumOp{T, N} <: ExprOp{T, N}
    subs::Vector{<:ExprOp{T, N}}
    SumOp{T, N}(s::Vector{<:ExprOp{T, N}}) where {T, N} =
        if isempty(s)
            error("bug: empty operator sum")
        elseif length(s) == 1
            return s[1]
        else
            return new{T, N}(s)
        end
end

sumsubs(a::SumOp) = a.subs
sumsubs(a::ExprOp) = [a]

(a::ExprOp{T, N} + b::ExprOp{T, N}) where {T, N} = SumOp{T, N}([ sumsubs(a) ; sumsubs(b)])
(a::ExprIndexed{Mixed} + b::ExprIndexed{Pure}) = a + Evolver(b)
(a::ExprIndexed{Pure} + b::ExprIndexed{Mixed}) = Evolver(a) + b
(a::ExprOp - b::ExprOp) = a + (-b)

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

ranking(::SumOp) = 21

isless(a::SumOp, b::SumOp) = isless(a.subs, b.subs)
        
apply_expr(f, a::SumOp) = reduce(+, map(f, a.subs))

    
############### Tensor products ############

struct TensorOp{T, N} <: ExprOp{T, N}
    subs::Vector{ExprOp{T}}
end

tensorsubs(a::TensorOp) = a.subs
tensorsubs(a::ExprOp) = [a]

"""
    op1 ⊗ op2

tensor product for generic operators, alternative syntax: tensor(op1, op2)

# Examples

    controlled(op) = Proj("Up") ⊗ Id + Proj("Dn") ⊗ op
    Rxy(t) = exp(im * t * (X ⊗ X + Y ⊗ Y) / 4)

"""
(a::ExprOp{T, N} ⊗ b::ExprOp{T, M}) where {T, N, M} =
    TensorOp{T, N + M}([tensorsubs(a) ; tensorsubs(b)])
tensor(a::ExprOp, b::ExprOp, c::ExprOp...) = ⊗(a, b, c...)

show(io::IO, a::TensorOp) =
    paren(io, Base.operator_precedence(:⊗)) do io
        join(io, a.subs, "⊗")
    end


ranking(::TensorOp) = 22

isless(a::TensorOp, b::TensorOp) = isless(a.subs, b.subs)

apply_expr(f, a::TensorOp) = TensorOp(map(f, a.subs))


############## Operator functions ###########

# PowOp

struct PowOp{T, N} <: ExprOp{T, N}
    arg::ExprOp{T, N}
    expo::Number
end

(a::ExprOp ^ b::Number) = PowOp(a, b)
(::ExprIndexed ^ ::Number) = error("cannot raise indexed expressions to a power, do it directly on operators")

show(io::IO, a::PowOp) =
    paren(io, Base.operator_precedence(:^)) do io
        print(io, a.arg, "^", a.expo)
    end

ranking(::PowOp) = 30

isless(a::PowOp, b::PowOp) = isless(a.arg, b.arg)

apply_expr(f, a::PowOp) = PowOp(f(a.arg), a.expo)

# ExpOp

struct ExpOp{T, N} <: ExprOp{T, N}
    arg::ExprOp{T, N}
end

"""
    exp(::ExprOp)

exponential for generic tensors
"""
exp(a::ExprOp) = ExpOp(a)
exp(::ExprIndexed) = error("cannot use exp on indexed expressions, use it directly on operators") 

show(io::IO, a::ExpOp) =
    paren(io, 1000, 0) do io
        show_func(io, "exp", a.arg)
    end

ranking(::ExpOp) = 31

isless(a::ExpOp, b::ExpOp) = isless(a.arg, b.arg)
    
apply_expr(f, a::ExpOp) = ExpOp(f(a.arg))

# SqrtOp

"""
    sqrt(::ExprOp)

square root for generic tensors
"""
struct SqrtOp{T, N} <: ExprOp{T, N}
    arg::ExprOp{T, N}
end

sqrt(a::ExprOp) = SqrtOp(a)
sqrt(::ExprIndexed) = error("cannot use sqrt on indexed expressions, use it directly on operators")

show(io::IO, a::SqrtOp) =
    paren(io, 1000, 0) do io
        show_func(io, "sqrt", a.arg)
    end


ranking(::SqrtOp) = 32

isless(a::SqrtOp, b::SqrtOp) = isless(a.arg, b.arg)

apply_expr(f, a::SqrtOp) = SqrtOp(f(a.arg))
       
# DagOp

struct DagOp{T, N} <: ExprOp{T, N}
    arg::ExprOp{T, N}
end

"""
    dag(::ExprOp)

adjoint for generic operators
"""
dag(a::ExprOp) = DagOp(a)
dag(a::DagOp) = a.arg

show(io::IO, a::DagOp) =
    paren(io, 1000, 0) do io
        show_func(io, "dag", a.arg)
    end

ranking(::DagOp) = 33

isless(a::DagOp, b::DagOp) = isless(a.arg, b.arg)

apply_expr(f, a::DagOp) = DagOp(f(a.arg))

eval_expr(f, a::Union{SumOp, ProdOp, TensorOp}) = any(f, a.subs)
eval_expr(f, a::Union{Gate, Left, Right, Evolver, Dissipator, PowOp, ExpOp, SqrtOp, DagOp}) = f(a.arg)
eval_expr(f, a::Indexed) = f(a.op)

