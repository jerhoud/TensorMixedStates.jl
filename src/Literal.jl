import Base: +, *, -, /, show, isless

export Operator, Lit, ProdLit, SumLit, Lits, @make_operators
export multipleLit, dissipLit, reorder, insertFfactors

struct Operator
    name::String
    opname::String
    dim::Int
    fermionic::Bool
    dissipator::Bool 
end

show(io::IO, op::Operator) =
    print(io, op.name)

"""
    struct Lit

    Represent a quantum operator
"""
struct Lit
    op::Operator
    index::Tuple
    param::NamedTuple
end

function show_func(io::IO, name, args, kwargs)
    print(io, name, "(")
    join(io, (repr(a) for a in args), ",")
    if !isempty(kwargs)
        print(io, ";")
        join(io, ("$sym=$(repr(val))" for (sym, val) in pairs(kwargs)), ",")
    end
    print(io, ")")
end

show(io::IO, a::Lit) =
    show_func(io, a.op.name, a.index, a.param)

isless(a::Lit, b::Lit) =
    isless((a.op.name, a.index, a.param), (b.op.name, b.index, b.param))


(oper::Operator)(index::Int...; kwargs...) =
    if length(index) == oper.dim
        ProdLit(1, [Lit(oper, Tuple(index), NamedTuple(kwargs))])
    else
        error("Operator $(oper.name) has dim = $(oper.dim) and was called with $(index)")
    end

(oper::Operator)(index::Index...; kwargs...) =   
    if length(index) == oper.dim
        op(oper.opname, index...; kwargs...)
    else
        error("Operator $(oper.name) has dim = $(oper.dim) and was called with $(index)")
    end


(l::Lit)(sites::Vector{Index}) =
    l.op((sites[i] for i in l.index)...; l.param...)
(l::Lit)(args...) = l.op(args...; l.param...)

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

function show(io::IO, a::ProdLit)
    print_coef(io, a.coef)
    join(io, a.ls)
end

isless(a::ProdLit, b::ProdLit) =
    isless((a.ls, a.coef), (b.ls, b.coef))

(a::Number * b::ProdLit) =
    ProdLit(a * b.coef, b.ls)
(a::ProdLit * b::Number) = b * a
(a::ProdLit * b::ProdLit) =
    ProdLit(a.coef * b.coef, [a.ls; b.ls])
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

function show(io::IO, a::SumLit)
    if a.ps == []
        print(io, 0)
        return
    end
    n = length(a.ps)
    i = 1
    while (i <= n)
        s = sprint(print, a.ps[i])
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
end

SumLit() = SumLit([])
SumLit(a::ProdLit) =
    if a.coef == 0
        SumLit()
    else 
        SumLit([a])
    end
SumLit(a::SumLit) = a

(a::SumLit + b::SumLit) = SumLit([a.ps; b.ps])
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
        push!(ps, ProdLit(ap.coef * bp.coef, [ap.ls; bp.ls]))
    end
    return SumLit(ps)
end

Lits = Union{ProdLit, SumLit}

"""
    @make_operators(names, dim::Int, fermionic::Bool, dissipator::Bool)

Define quantum operators with the given names (names should be iterable)

These operators are callable on `Int...` or on `Index...` to produce Lit or ITensor
"""
macro make_operators(names, dim::Int, fermionic::Bool, dissipator::Bool)
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
            $(esc(sname)) = Operator($name, $opname, $dim, $fermionic, $dissipator)
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
dissipLit(a::Lit) = a.op.dissipator
dissipLit(a::ProdLit) = any(dissipLit, a.ls)
dissipLit(a::SumLit) = any(dissipLit, a.ps)



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

"""
    reorder(ops)

Put the operators in good order and simplifies for use in mpo creation.
Take in account fermionic nature of the ops.
Work only if `multipleLit(ops)==false``.
Not to be used directly.

# Examples
```julia-repl
julia> reorder(X(3)Y(1) + X(1) + Y(5)Z(2)X(1) - X(1))
X(1)Z(2)Y(5)+Y(1)X(3)

julia> reorder(C(3)C(1)+Cdag(4)Cdag(1)+2C(1)C(3))
C(1)C(3)+(-1)Cdag(1)Cdag(4)
```
"""
function reorder(a::ProdLit)
    fs = filter(t->t.op.fermionic, a.ls)
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

const opF = Operator("F", "F", 1, false, false)
litF(idx) = Lit(opF, (idx,), (;))

function insertFfactorsCore(a::ProdLit)
    nls = Lit[]
    fermion_idx = 0
    for l in reverse(a.ls)
        idx = first(l.index)
        if fermion_idx ≠ 0
            append!(nls, (litF(i) for i in reverse(idx + 1:fermion_idx - 1)))
            if l.op.fermionic
                push!(nls, litF(idx))
                fermion_idx = 0
            else
                fermion_idx = idx
            end
        elseif l.op.fermionic && !l.op.dissipator
            fermion_idx = idx
        end
        push!(nls, l)
    end
    if fermion_idx ≠ 0
        append!(nls, (litF(i) for i in reverse(1:fermion_idx - 1)))
    end
    return ProdLit(a.coef, reverse(nls))
end

insertFfactors(a::ProdLit) =
    insertFfactorsCore(reorder(a))

"""
    insertFfactors(ops)

Reorder ops and insert the needed F factors according to fermion string rule.
Not to be used directly.

# Examples
```julia-repl
julia> insertFfactors(C(3)C(1))
(-1)C(1)F(1)F(2)C(3)

julia> insertFfactors(C(1)C(3)C(5)C(7))
C(1)F(1)F(2)C(3)C(5)F(5)F(6)C(7)
```
"""
insertFfactors(a::SumLit) =
    SumLit(map(insertFfactorsCore, (reorder(a)).ps))

insertFfactors(a) = map(insertFfactors, a)