import Base: +, *, -, /, show

struct Lit
    name::String
    index::Tuple
    param::NamedTuple
    fermionic::Bool
end

function show(io::IO, a::Lit)
    print(io, "$(a.name)(")
    join(io, a.index, ", ")
    print(io, ")")
end

struct ProdLit
    coef::Number
    ls::Vector{Lit}
    ProdLit(c, l) =
        if c == 0
            new(0, [])
        else
            new(c, l)
        end
    ProdLit() = new(0, [])
end

function show(io::IO, a::ProdLit)
    if a.coef ≠ 1
        if isa(a.coef, Complex) || a.coef < 0
            print(io, "($(a.coef))")
        else
            print(io, "$(a.coef)")
        end
    end
    join(io, a.ls)
end


(a::Number * b::ProdLit) =
    ProdLit(a * b.coef, b.ls)
(a::ProdLit * b::Number) = b * a
(a::ProdLit * b::ProdLit) =
    ProdLit(a.coef * b.coef, vcat(a.ls, b.ls))
(a::ProdLit / b::Number) = (1 / b) * a
-(a::ProdLit) = -1 * a

struct SumLit
    ps::Vector{ProdLit}
end

show(io::IO, a::SumLit) =
    if a.ps == []
        print(io, 0)
    else
        join(io, a.ps, "+")
    end

SumLit() = SumLit([])
SumLit(a::ProdLit) =
    if a.coef == 0
        SumLit()
    else 
        SumLit([a])
    end
SumLit(a::SumLit) = a

(a::SumLit + b::SumLit) = SumLit(vcat(a.ps, b.ps))
(a::SumLit + b::ProdLit) = a + SumLit(b)
(a::ProdLit + b::SumLit) = SumLit(a) + b
(a::ProdLit + b::ProdLit) = SumLit(a) + SumLit(b)
(a::ProdLit - b::ProdLit) = a + (-1) * b

(a::Number * b::SumLit) =
    if (a==0) SumLit()
    else SumLit(map(pl->ProdLit(a * pl.coef, pl.ls), b.ps))
    end
(a::SumLit * b::Number) = b * a
(a::SumLit * b::ProdLit) = a * SumLit(b)
(a::ProdLit * b::SumLit) = SumLit(a) * b
(a::SumLit / b::Number) = (1 / b) * a

-(a::SumLit) = -1 * a

(a::SumLit - b) = a + (-1) * b
(a - b::SumLit) = a + (-1) * b
(a::SumLit - b::SumLit) = a + (-1) * b

function (a::SumLit * b::SumLit)
    ps = ProdLit[]
    for ap in a.ps, bp in b.ps
        push!(ps, ProdLit(ap.coef * bp.coef, vcat(ap.ls, bp.ls)))
    end
    return SumLit(ps)
end

Lits = Union{ProdLit, SumLit}

macro opLit(name::String, fermionic::Bool = false)
    sname = Symbol(name)
    if isdefined(TensorMixedStates, sname)
        return quote end
    else
        return quote
            $(esc(sname))(index...; kwargs...) = ProdLit(1, [Lit($name, Tuple(index), NamedTuple(kwargs), $fermionic)])
            export $(esc(sname))
        end
    end
end

Lit_to_OpSum(a::SumLit) =
    sum(p.coef * prod(Op(l.name, l.index...; l.param...) for l in p.ls) for p in a.ps; init=OpSum())

Lit_to_OpSum(a::ProdLit) = Lit_to_OpSum(SumLit(a))

function Lit_to_ops(a::ProdLit, sites)
    r = [op(sites, l.name, l.index...; l.param...) for l in a.ls]
    r[1] *= a.coef
    return r
end
