import Base: +, *, -, /, show, isless

export Lit, ProdLit, SumLit, Lits, @opLit, multipleLit, dissipLit

"""
    struct Lit

    Represent a quantum operator
"""
struct Lit
    name::String
    opname::String
    index::Tuple
    param::NamedTuple
    fermionic::Bool
    dissipator::Bool
end

function show(io::IO, a::Lit)
    print(io, "$(a.name)(")
    join(io, a.index, ", ")
    print(io, ")")
end

isless(a::Lit, b::Lit) =
    isless((a.name, a.index, a.param), (b.name, b.index, b.param))

"""
    struct ProdLit
    ProdLit()
    ProdLit(::Number, ::Vector{Lit})

Represent a product of quantum operators with a scalar factor

# Fields
- `coef::Number`: scalar factor
- `ls::Vector{Lit}`: quantum operator factors
"""
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

(a::Number * b::ProdLit) =
    ProdLit(a * b.coef, b.ls)
(a::ProdLit * b::Number) = b * a
(a::ProdLit * b::ProdLit) =
    ProdLit(a.coef * b.coef, vcat(a.ls, b.ls))
(a::ProdLit / b::Number) = (1 / b) * a
-(a::ProdLit) = -1 * a

"""
    struct SumLit
    SumLit()
    SumLit(::Vector{ProdLit})

Represent a sum of products of quantum oprators

# Fields
- `ps::Vector{ProdLit}`: terms of the sum
"""
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

"""
    @opLit(names, fermionic::Bool=false, dissipator::Bool=false)

Define quantum operators with the given names (names should be iterable)

These operators are callable on `Int...` or on `Index...` to produce Lit or ITensor
"""
macro opLit(names, fermionic::Bool=false, dissipator::Bool=false)
    ns = eval(names)
    e = Expr(:block)
    for n in ns
        if n isa String
            name = opname = n
        else
            name = n[1]
            opname = n[2]
        end
        sname = Symbol(name)
        push!(e.args,
        quote
            $(esc(sname))() = $name
            $(esc(sname))(index::Int...; kwargs...) =
                ProdLit(1, [Lit($name, $opname, Tuple(index), NamedTuple(kwargs), $fermionic, $dissipator)])
            $(esc(sname))(index::Index...; kwargs...) =
                op($opname, index...; kwargs...)
            export $(esc(sname))
        end)
    end
    return e
end

"""
    multipleLit(::ProdLit)
    multipleLit(::SumLit)

Return true if any contained Lit has multiple indices
"""
multipleLit(a::Lit) = length(a.index) ≠ 1
multipleLit(a::ProdLit) = any(multipleLit, a.ls)
multipleLit(a::SumLit) = any(multipleLit, a.ps)

"""
    dissipLit(a::ProdLit)
    dissipLit(a::SumLit)

Return true if any contained Lit is a dissipator
"""
dissipLit(a::Lit) = a.dissipator
dissipLit(a::ProdLit) = any(dissipLit, a.ls)
dissipLit(a::SumLit) = any(dissipLit, a.ps)
