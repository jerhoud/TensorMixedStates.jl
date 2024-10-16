import Base: +, *, -, /, show

struct Lit
    name::String
    index::Tuple
    param::NamedTuple
end

function show(io::IO, a::Lit)
    print(io, "$(a.name)(")
    join(io, a.index, ", ")
    print(io, ")")
end

struct ProdLit
    coef::Number
    ls::Vector{Lit}
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
    if (a==0) 0 else ProdLit(a * b.coef, b.ls) end
(a::ProdLit * b::Number) = b * a
(a::ProdLit * b::ProdLit) = ProdLit(a.coef * b.coef, vcat(a.ls, b.ls))
(a::ProdLit / b::Number) = (1 / b) * a
-(a::ProdLit) = -1 * a

struct SumLit
    ps::Vector{ProdLit}
end

show(io::IO, a::SumLit) = join(io, a.ps, "+")

SumLit(a::ProdLit) = SumLit([a])
SumLit(a::SumLit) = a

(a::SumLit + b::SumLit) = SumLit(vcat(a.ps, b.ps))
(a::SumLit + b::ProdLit) = a + SumLit(b)
(a::ProdLit + b::SumLit) = SumLit(a) + b
(a::ProdLit + b::ProdLit) = SumLit([a, b])
(a::ProdLit - b::ProdLit) = a + (-1) * b

(a::Number * b::SumLit) =
    if (a==0) 0
    else SumLit(map(pl->ProdLit(a * pl.coef, pl.ls), b.ps))
    end
(a::SumLit * b::Number) = b * a
(a::SumLit * b::ProdLit) = a * SumLit(b)
(a::ProdLit * b::SumLit) = SumLit(a) * b
(a::SumLit / b::Number) = (1 / b) * a

-(a::SumLit) = -1 * a

(a::SumLit - b) = a + (-1) * b
(a - b::SumLit) = a + (-1) * b

function (a::SumLit * b::SumLit)
    ps = ProdLit[]
    for ap in a.ps, bp in b.ps
        push!(ps, ProdLit(ap.coef * bp.coef, vcat(ap.ls, bp.ls)))
    end
    return SumLit(ps)
end

Lits = Union{ProdLit, SumLit}

macro opLit(name::String)
    sname = Symbol(name)
    if isdefined(TensorMixedStates, sname)
        return quote end
    else
        return quote
            $(esc(sname))(index...; kwargs...) = ProdLit(1, [Lit($name, Tuple(index), NamedTuple(kwargs))])
            export $(esc(sname))
        end
    end
end

function Lit_to_OpSum(a::SumLit, ext::String="")
    return sum(p.coef * prod(Op(ext * l.name, l.index...; l.param...) for l in p.ls) for p in a.ps; init=OpSum())
end

Lit_to_OpSum(a::ProdLit, ext::String="") = Lit_to_OpSum(SumLit(a), ext)

function Lit_to_ops(a::ProdLit, sites, ext::String="")
    r = [op(sites, ext * l.name, l.index...; l.param...) for l in a.ls]
    r[1] *= a.coef
    return r
end
