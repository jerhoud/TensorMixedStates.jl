export ExprOp, ProdOp, SumOp, TensorOp, ExpOp, Operator, ⊗, Gate, Dissipator

abstract type ExprOp{N} end

struct Operator{N} <: ExprOp{N}
    name::String
    expr::Union{Nothing, Matrix, Function, ExprOp{N}}
end
Operator{N}(name) where N = Operator{N}(name, nothing)
Operator(name) = Operator{1}(name)
Operator(name, expr) = Operator{1}(name, expr)


show(io::IO, op::Operator) =
    print(io, op.name)

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

struct ProdOp{N} <: ExprOp{N}
    coef::Number
    subs::Vector{<:ExprOp{N}}
end

function show(io::IO, a::ProdOp)
    print_coef(io, a.coef)
    join(io, a.subs, "*")
end

prodsubs(a::ProdOp) = a.subs
prodsubs(a::ExprOp) = [a]
coefsubs(a::ProdOp) = a.coef
coefsubs(a::ExprOp) = 1

(a::ExprOp{N} * b::ExprOp{N}) where N =
    ProdOp{N}(coefsubs(a) * coefsubs(b), [prodsubs(a) ; prodsubs(b)])
(a::Number * b::ExprOp{N}) where N =
    ProdOp{N}(a * coefsubs(b), prodsubs(b))
(a::ExprOp * b::Number) = b * a
(a::ExprOp / b::Number) = inv(b) * a
-(a::ExprOp) = -1 * a

struct SumOp{N} <: ExprOp{N}
    subs::Vector{<:ExprOp{N}}
end

function show(io::IO, a::SumOp)
    print(io, "(")
    n = length(a.subs)
    i = 1
    while (i <= n)
        s = sprint(print, a.subs[i])
        if i == 1
            print(io, s)
        else
            if s[1] ≠ '-'
                print(io, "+")
            elseif i == 4 && n > 6
                print(io, "-")
            end
            if i == 4 && n > 6
                print(io, "...")
                i = n - 1
                continue
            else
                print(io, s)
            end
        end
        i = i + 1
    end
    print(io, ")")
end

sumsubs(a::SumOp) = a.subs
sumsubs(a::ExprOp) = [a]

(a::ExprOp{N} + b::ExprOp{N}) where N =
    SumOp{N}([sumsubs(a) ; sumsubs(b)])
(a::ExprOp - b::ExprOp) = a + (-b)

struct TensorOp{N} <: ExprOp{N}
    subs::Vector{<:ExprOp}
end

function show(io::IO, a::TensorOp)
    print(io, "(")
    join(io, a.subs, "⊗")
    print(io, ")")
end

tensorsubs(a::TensorOp) = a.subs
tensorsubs(a::ExprOp) = [a]

(a::ExprOp{N} ⊗ b::ExprOp{M}) where {N, M} =
    TensorOp{N + M}([tensorsubs(a) ; tensorsubs(b)])

struct ExpOp{N} <: ExprOp{N}
    arg::ExprOp{N}
end

show(io::IO, a::ExpOp) = show_func(io, "exp", a.arg)

exp(a::ExprOp{N}) where N = ExpOp{N}(a)

struct Gate{N} <: ExprOp{N}
    arg::ExprOp{N}
end

show(io::IO, a::Gate) = show_func(io, "Gate", a.arg)

struct Dissipator{N} <: ExprOp{N}
    arg::ExprOp{N}
end

show(io::IO, a::Dissipator) = show_func(io, "Dissipator", a.arg)