export ExprOp, ProdOp, SumOp, TensorOp, ExpOp, Operator, Indexed, IndexOp, ⊗, Gate, Dissipator

abstract type ExprOp{N} end

struct ProdOp{N} <: ExprOp{N}
    coef::Number
    subs::Vector{<:ExprOp{N}}
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

sumsubs(a::SumOp) = a.subs
sumsubs(a::ExprOp) = [a]

(a::ExprOp{N} + b::ExprOp{N}) where N =
    SumOp{N}([sumsubs(a) ; sumsubs(b)])
(a::ExprOp - b::ExprOp) = a + (-b)

struct TensorOp{N} <: ExprOp{N}
    subs::Vector{<:ExprOp}
end

tensorsubs(a::TensorOp) = a.subs
tensorsubs(a::ExprOp) = [a]

(a::ExprOp{N} ⊗ b::ExprOp{M}) where {N, M} =
    TensorOp{N + M}([tensorsubs(a) ; tensorsubs(b)])

struct ExpOp{N} <: ExprOp{N}
    arg::ExprOp{N}
end

exp(a::ExprOp{N}) where N = ExpOp{N}(a)

struct Gate{N} <: ExprOp{N}
    arg::ExprOp{N}
end

struct Dissipator{N} <: ExprOp{N}
    arg::ExprOp{N}
end




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

show(io::IO, a::ProdOp) =
    paren(io, Base.operator_precedence(:*)) do io
        print_coef(io, a.coef)
        join(io, a.subs, "*")
    end

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


show(io::IO, a::TensorOp) =
    paren(io, Base.operator_precedence(:⊗)) do io
        join(io, a.subs, "⊗")
    end

show(io::IO, a::ExpOp) =
    paren(io, 1000) do io
        show_func(io, "exp", a.arg)
    end

#show(io::IO, a::Gate) = show_func(io, "Gate", a.arg)

#show(io::IO, a::Dissipator) = show_func(io, "Dissipator", a.arg)



struct Operator{N} <: ExprOp{N}
    name::String
    expr::Union{Nothing, Matrix, Function, ExprOp{N}}
    fermionic::Bool
end

Operator{N}(name, expr = nothing) where N = Operator{N}(name, expr, false)
Operator(name, expr = nothing, fermionic = false) = Operator{1}(name, expr, fermionic)

show(io::IO, op::Operator) =
    print(io, op.name)

struct IndexOp end

struct Indexed{N} <: ExprOp{IndexOp}
    op::ExprOp{N}
    index::NTuple{N, Int}
end

(op::ExprOp{N})(index::Vararg{Int, N}) where N =
    Indexed{N}(op, index)

show(io::IO, ind::Indexed) =
    if ind.op isa Operator
        show_func(io, ind.op.name, collect(ind.index))
    else
        show_func(io, repr(ind.op; context=:precedence=>500), collect(ind.index))
    end