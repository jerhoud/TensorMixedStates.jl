import Base: +, *, -, /, show, isless, ==

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

isless(a::Lit, b::Lit) =
    isless((a.name, a.index), (b.name, b.index))

(a::Lit == b::Lit) = 
    (a.name, a.index) == (b.name, b.index)

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

isless(a::ProdLit, b::ProdLit) =
    isless((a.ls, a.coef), (b.ls, b.coef))

(a::ProdLit == b::ProdLit) =
    (a.ls, a.coef) == (b.ls, b.coef)

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

function signature(p)
    n = length(p)
    t = fill(false, n)
    r = 1
    for i in 1:n
        if t[i] continue end
        j = p[i]
        while j ≠ i
            r = -r
            t[j] = true
            j = p[j]
        end
    end
    return r
end

function reorder(a::ProdLit)
    fs = filter(t->t.fermionic, a.ls)
    p = sortperm(fs; by=l->l.index)
    s = signature(p)
    ProdLit(s * a.coef, sort(a.ls; by=l->l.index))
end

function reorder(a::SumLit)
    ps = sort(map(reorder, a.ps))
    pps = ProdLit[]
    c = 0
    l = Lit[]
    for p in ps
        if p.ls == l
            c += p.coef
        else
            if c ≠ 0
                push!(pps, ProdLit(c, l))
            end
            c = p.coef
            l = p.ls
        end
    end
    if c ≠ 0
        push!(pps, ProdLit(c, l))
    end
    return SumLit(pps)
end

litF(idx) = Lit("F", (idx,), (;), false)

function insertFfactors(a::ProdLit)
    nls = []
    fermion_idx = 0
    for l in reverse(a.ls)
        idx = first(l.index)
        if fermion_idx ≠ 0
            append!(nls, (litF(i) for i in reverse(idx + 1:fermion_idx - 1)))
            if l.fermionic
                push!(nls, litF(idx))
                fermion_idx = 0
            else
                fermion_idx = idx
            end
        elseif l.fermionic
            fermion_idx = idx
        end
        push!(nls, l)
    end
    if fermion_idx ≠ 0
        append!(nls, (litF(i) for i in reverse(1:fermion_idx - 1)))
    end
    return ProdLit(a.coef, reverse(nls))
end