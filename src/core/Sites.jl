# Site types: what every site provides, its dimension, the parity of its basis states, and the
# states and operators declared for its type, which become ordinary methods.

export AbstractSite, dim, fermion_parities, fermion_species, string_state
export Level, F, Swap, controlled, @def_states, @def_operators, @create_site_module

"""
    abstract type AbstractSite

the supertype of all site types. A site type defines `dim`, declares its states and operators
with `@def_states` and `@def_operators`, and, if fermionic, `fermion_parities`. It keeps the
constructor taking the values of its fields in their order, which reading a state file and
`weaken` use.

# Examples

    struct MySite <: AbstractSite
        conserve::Conservation
    end
    MySite(; conserve = ()) = conserving(MySite; conserve)
    TensorMixedStates.dim(::MySite) = 2

# Operators

Declared once for every site type, each of these holds on a site type that declares what it is
made of, `Sp` and `Sz`, `A`, or `Cup` and `Cdn`, unless it declares the name itself:

- `Sm`             : `dag(Sp)`
- `Sx, Sy`         : `(Sp + Sm) / 2` and `(Sp - Sm) / 2im`
- `S2`             : `(Sp * Sm + Sm * Sp) / 2 + Sz^2`
- `Q, P`           : `(A + dag(A)) / √2` and `(A - dag(A)) / (im * √2)`
- `Nup, Ndn, Ntot` : `dag(Cup) * Cup`, `dag(Cdn) * Cdn` and `Nup + Ndn`
"""
abstract type AbstractSite end

"""
    dim(site)

the dimension of the Hilbert space of a site, which every site type defines.
"""
dim(site::AbstractSite) = error("dim not implemented on site $site")

"""
    fermion_parities(site)

the fermionic parity, 0 or 1, of each basis state of `site`, all 0 on a site that is not
fermionic. A fermionic site type defines it, `F` being derived from it.

# Examples

    TensorMixedStates.fermion_parities(::MySite) = (0, 1)
"""
fermion_parities(site::AbstractSite) = fill(0, dim(site))

"""
    fermion_species(site)

the annihilation operators of the fermion species of a site, as `(C,)` or `(Cup, Cdn)`, for
`slater_state` and `fermi_sea`; `nothing` when the site type declares none.
"""
fermion_species(::AbstractSite) = nothing

"""
    Level

the operator numbering the basis states of a site, ``\\mathrm{diag}(0, 1, \\dots, d - 1)``,
defined on every site type.

# Examples

    State{Pure}(System(3, Qudit(4)), Level .=> [0, 3, 1])
    Proj(Level => 2)
"""
const Level = Operator{1}("Level", s -> diagm(Float64.(0:dim(s) - 1)))

"""
    F

the fermionic parity of a site, ``(-1)^p`` on a basis state of parity `p`, see
`fermion_parities`: the identity on a site that is not fermionic.
"""
const F = Operator{1}("F", s -> diagm(Float64[ (-1)^p for p in fermion_parities(s) ]))

"""
    swap_matrix(a, b)

the matrix of `Swap` on the sites `a` and `b`, refused when they are not of one type and one
dimension
"""
function swap_matrix(a::AbstractSite, b::AbstractSite)
    d = dim(a)
    if typeof(a) ≠ typeof(b) || dim(b) ≠ d
        error("Swap exchanges the states of two sites of one type and one dimension, not those " *
              "of $a and $b")
    end
    p = fermion_parities(a)
    m = zeros(d^2, d^2)
    for i in 1:d, j in 1:d
        m[(j - 1) * d + i, (i - 1) * d + j] = (-1)^(p[i] * p[j])
    end
    return m
end

"""
    Swap

the exchange of the states of two sites of one type, ``|a, b\\rangle \\mapsto (-1)^{p_a p_b}
|b, a\\rangle``, the sign coming between two odd states, see `fermion_parities`.

# Examples

    Swap(1, 2)
    Proj(Swap => -1)       # the antisymmetric part of two sites
"""
const Swap = Operator{2}("Swap", swap_matrix)

"""
    controlled_name(op)

the default name of `controlled(op)`: `C` followed by the name of an `Operator`, as `CZ`, and
`controlled(op)` for an expression
"""
controlled_name(a::Operator) = "C" * a.name
controlled_name(a::Op) = "controlled($a)"

"""
    controlled(op; name)

the gate leaving the following sites as they are when the first one is in its first basis
state, `Level => 0`, and applying `op` to them otherwise: on a qubit, `op` when the first is in
`"Dn"`, `"1"`. Its name is `"C"` followed by the name of a named operator, `CZ` for
`controlled(Z)`, and `"controlled(op)"` otherwise.

# Examples

    CZ = controlled(Z)
    Toffoli = controlled(controlled(X))
    CZ(1, 2)
"""
controlled(op::GenericOp{IsOperator, N}; name::String = controlled_name(op)) where N =
    Operator{N+1}(name, Proj(Level => 0) ⊗ IdentityOp(op) + (Id - Proj(Level => 0)) ⊗ op)

"""
    operator_definition(site, ::Val{name})
    state_definition(site, ::Val{name})

the definition `@def_operators` or `@def_states` gave the operator or the state `name` on the
type of `site`, and `nothing` when there is none
"""
operator_definition(::AbstractSite, ::Val) = nothing
state_definition(::AbstractSite, ::Val) = nothing

"""
    string_state(site, name)

the local state a name gives by a rule of the site rather than by a declaration, or `nothing`
when the name is not one of its forms. By default `"0"` gives the first basis state, `"1"` the
second, and so on. A site type overloads it to read names of its own, as `Spin` reads `"1/2"`,
or to read none.

# Examples

    TensorMixedStates.string_state(::MySite, ::String) = nothing
"""
function string_state(site::AbstractSite, name::String)
    n = tryparse(Int, name)
    if isnothing(n) || !(0 ≤ n < dim(site))
        return nothing
    end
    v = fill(0., dim(site))
    v[n + 1] = 1.
    return v
end

"""
    common_states

the states every site has, as functions of the site, looked up after those the site declares
or names by its rule: `"FullyMixed"`, the density matrix proportional to the identity
"""
const common_states = Dict{String, Function}(
    "FullyMixed" => s -> Matrix{Float64}(I, dim(s), dim(s)) / dim(s),
)

"""
    rebuild_site(t, values, why)

the site of type `t` built from the values of its fields, refused with a message naming `why`
when `t` has no constructor taking them
"""
function rebuild_site(t::Type{<:AbstractSite}, values, why::String)
    if !hasmethod(t, Tuple{map(typeof, values)...})
        error("$t must keep its constructor $t($(join(fieldnames(t), ", "))) for $why")
    end
    return t(values...)
end

"""
    declared_type(t, what)

the site type a declaration is for, refused when given as a site, as version 1 did
"""
function declared_type(t, what::String)
    if t isa Expr && t.head == :call
        error("@$what takes the site type, not a site: write @$what $(t.args[1]) [...]")
    end
    return t
end

"""
    @def_states(T, [name => def, ...])

declare states for the site type `T`. A name is a string, or a vector of strings naming the
same state. A definition is a vector (a pure state in the basis of the site), a matrix (a
density matrix), a function of the site returning either, or the name of another state.

# Examples

    @def_states MySite [
        ["Emp", "0"] => [1, 0],
        "Occ" => [0, 1],
    ]
"""
macro def_states(t, defs)
    t = declared_type(t, "def_states")
    if !(defs isa Expr) || defs.head ≠ :vect
        error("@def_states takes a vector of declarations, name => definition")
    end
    e = Expr(:block)
    seen = String[]
    for d in defs.args
        if !(d isa Expr && d.head == :call && d.args[1] == :(=>))
            error("@def_states takes declarations name => definition, not $d")
        end
        # literal names, a declaration becoming methods when the macro expands
        names = d.args[2]
        if names isa String
            names = [names]
        elseif names isa Expr && names.head == :vect && all(x -> x isa String, names.args)
            names = String.(names.args)
        else
            error("@def_states names a state by a literal string or a vector of them, not $names")
        end
        for name in names
            if name in seen
                error("@def_states declares the state $name twice")
            end
            push!(seen, name)
            push!(e.args, :($(GlobalRef(@__MODULE__, :state_definition))(
                ::$(esc(t)), ::Val{$(QuoteNode(Symbol(name)))}) = $(esc(d.args[3]))))
        end
    end
    return e
end

"""
    check_shared_operator(existing, name)

refuse to declare the operator `name` when the name already stands for something else than
the operator of that name each site type defines
"""
check_shared_operator(existing, name::String) =
    if existing isa Operator{1} && existing.name == name && !isnothing(existing.expr)
        error("cannot declare operator $name: the name already stands for an operator with " *
              "a definition of its own")
    elseif !(existing isa Operator{1} && existing.name == name)
        error("cannot declare operator $name: the name already stands for $(repr(existing))")
    end

"""
    @def_operators(T, [name = def, ...])

declare operators for the site type `T`. A definition is a matrix, an expression of operators
of one site, or a function of the site returning either. Each name becomes a `const` of the
calling module the first time it is seen, so the macro is used at the top level; a name
already standing for the operator of that name is shared between site types.

# Examples

    @def_operators MySite [
        C = [0 1 ; 0 0],
        N = dag(C) * C,
    ]
"""
macro def_operators(t, defs)
    t = declared_type(t, "def_operators")
    if !(defs isa Expr) || defs.head ≠ :vect
        error("@def_operators takes a vector of declarations, name = definition")
    end
    e = Expr(:block)
    seen = Symbol[]
    for d in defs.args
        if !(d isa Expr && d.head == :(=) && d.args[1] isa Symbol)
            error("@def_operators takes declarations name = definition, not $d")
        end
        sym, def = d.args
        if sym == :F
            error("F is not declared: a fermionic site type gives the parity of its basis " *
                  "states with fermion_parities, F being derived from it")
        elseif sym in seen
            error("@def_operators declares the operator $sym twice")
        end
        push!(seen, sym)
        name = string(sym)
        if isdefined(__module__, sym)
            # not bound again: the name is shared with the site types already declaring it
            push!(e.args, :($(GlobalRef(@__MODULE__, :check_shared_operator))($(esc(sym)), $name)))
        else
            push!(e.args, :(const $(esc(sym)) = $(GlobalRef(@__MODULE__, :Operator)){1}($name)))
        end
        push!(e.args, :($(GlobalRef(@__MODULE__, :operator_definition))(
            ::$(esc(t)), ::Val{$(QuoteNode(sym))}) = $(esc(def))))
    end
    return e
end

"""
    @create_site_module(name, symbols)

define a submodule `name` importing the given symbols from the module the macro is called in
and exporting them, so that `using .name` brings a site type and its operators into scope.

# Examples

    @create_site_module(Spins, [Spin, Sp, Sm, Sx, Sy, Sz, S2, N])
"""
macro create_site_module(name, symbols)
    if !(symbols isa Expr) || symbols.head ≠ :vect
        error("@create_site_module takes a vector of symbols")
    end
    imports = [ Expr(:., :., :., s) for s in symbols.args ]
    block = Expr(:block, Expr(:import, imports...), Expr(:export, symbols.args...))
    doc = "    using .$name\n\nmodule giving access to the `$(symbols.args[1])` site type " *
        "and the names it exports: " * join(map(s -> "`$s`", symbols.args[2:end]), ", ")
    mod = Expr(:module, true, name, block)
    return esc(Expr(:macrocall, GlobalRef(Core, Symbol("@doc")), __source__, doc, mod))
end

############### Relations declared once for every site type ###############

# what the relations below are made of, each site type defining it in its own way
const Sp = Operator{1}("Sp")
const Sz = Operator{1}("Sz")
const A = Operator{1}("A")
const Cup = Operator{1}("Cup")
const Cdn = Operator{1}("Cdn")

# a relation holds on a site type declaring what it is made of, and a declaration of the same
# name by the site type wins over it
@def_operators AbstractSite [
    Sm = dag(Sp),
    Sx = (Sp + Sm) / 2,
    Sy = (Sp - Sm) / (2im),
    S2 = (Sp * Sm + Sm * Sp) / 2 + Sz^2,
    Q = (A + dag(A)) / √2,
    P = (A - dag(A)) / (im * √2),
    Nup = dag(Cup) * Cup,
    Ndn = dag(Cdn) * Cdn,
    Ntot = Nup + Ndn,
]

############### The help of the operators of the site types ###############

"""
    site_types()

the site types declaring operators, read off the methods `@def_operators` added, in the order of
their names, and last the abstract ones declaring some, as `AbstractSite` its relations
"""
function site_types()
    ts = DataType[]
    for m in methods(operator_definition)
        sig = Base.unwrap_unionall(m.sig)
        if length(sig.parameters) == 3
            t, v = sig.parameters[2], sig.parameters[3]
            if t isa DataType && t <: AbstractSite && !(t in ts) &&
               (!isabstracttype(t) || (v isa DataType && !isempty(v.parameters)))
                push!(ts, t)
            end
        end
    end
    return sort!(ts; by = t -> (isabstracttype(t), string(nameof(t))))
end

"""
    declaring_types(name)

the names of the site types declaring the operator `name`, for a message
"""
function declaring_types(name::String)
    v = Val{Symbol(name)}
    ts = String[]
    for m in methods(operator_definition)
        sig = Base.unwrap_unionall(m.sig)
        if length(sig.parameters) == 3 && sig.parameters[3] === v &&
           sig.parameters[2] isa DataType && !isabstracttype(sig.parameters[2])
            push!(ts, string(nameof(sig.parameters[2])))
        end
    end
    return sort!(unique!(ts))
end

"""
    own_rule(site)

whether the site type of `site` names states by a rule of its own, a method of `string_state`
"""
own_rule(site::AbstractSite) =
    Base.unwrap_unionall(which(string_state, Tuple{typeof(site), String}).sig).parameters[2] ≠
    AbstractSite

"""
    known_states(site)

the names of the states of `site`, for a message: those it declares and those of every site,
with `"0"` to `"d-1"` when it keeps the default rule of `string_state`, its help named when it
has a rule of its own
"""
function known_states(site::AbstractSite)
    names = String[]
    for m in methods(state_definition)
        sig = Base.unwrap_unionall(m.sig)
        if length(sig.parameters) == 3
            t, v = sig.parameters[2], sig.parameters[3]
            if t isa Type && site isa t && v isa DataType && v <: Val && !isempty(v.parameters)
                push!(names, string(v.parameters[1]))
            end
        end
    end
    if own_rule(site)
        return join([ sort!(names); sort!(collect(keys(common_states))) ], ", ") *
               " (see ?$(nameof(typeof(site))))"
    end
    levels = string.(0:dim(site) - 1)
    named = sort!(filter(n -> n ∉ levels, names))
    return join([ named; sort!(collect(keys(common_states))) ], ", ") *
           " and \"0\" to \"$(dim(site) - 1)\""
end

"""
    markdown_blocks(md)

the blocks of the help `md`, those of the helps it gathers taken one by one
"""
markdown_blocks(md::Markdown.MD) =
    reduce(vcat, [ markdown_blocks(x) for x in md.content ]; init = Any[])
markdown_blocks(x) = Any[x]

"""
    type_help(t)

the help of the type `t`, its docstrings read off the documentation of its module, where
`Base.Docs.doc` needs the REPL
"""
function type_help(t::DataType)
    m = parentmodule(t)
    docs = get(Base.Docs.meta(m), Base.Docs.Binding(m, nameof(t)), nothing)
    if isnothing(docs)
        return Markdown.MD(Any[])
    end
    return Markdown.MD(Any[ Markdown.parse(join(string.(docs.docs[s].text))) for s in docs.order ])
end

"""
    operator_lines(t)

the items of the section `Operators` of the help of the site type `t`, as the site types of
this package write them: the names given as code before the colon of each, and what follows it
"""
function operator_lines(t::DataType)
    lines = Tuple{Vector{String}, Vector{Any}}[]
    inside = false
    for b in markdown_blocks(type_help(t))
        if b isa Markdown.Header
            inside = b.text == ["Operators"]
        elseif inside && b isa Markdown.List
            for item in b.items
                if isempty(item) || !(first(item) isa Markdown.Paragraph)
                    continue
                end
                names, rest, colon = String[], Any[], false
                for x in first(item).content
                    if colon
                        push!(rest, x)
                    elseif x isa Markdown.Code
                        append!(names, strip.(split(x.code, ",")))
                    elseif x isa AbstractString && occursin(":", x)
                        colon = true
                        after = lstrip(split(x, ":"; limit = 2)[2])
                        if !isempty(after)
                            push!(rest, after)
                        end
                    end
                end
                push!(lines, (names, rest))
            end
        end
    end
    return lines
end

"""
    getdoc(op::Operator{1})

the help of an operator of one site that has none of its own, as `X`, `?X` at the prompt: what
it is on each site type naming it in the section `Operators` of its help, then on every site
type as `AbstractSite` names it, and the modules of this package exporting it; `nothing`, which
leaves the help to Julia, for another operator
"""
function Base.Docs.getdoc(a::Operator{1})
    name = Symbol(a.name)
    ts = site_types()
    homes = unique([@__MODULE__; [ parentmodule(t) for t in ts ]])
    bound = [ m for m in homes if isdefined(m, name) && getfield(m, name) === a ]
    if isempty(bound) || any(m -> haskey(Base.Docs.meta(m), Base.Docs.Binding(m, name)), bound)
        return nothing
    end
    lines = Any[]
    for t in ts, (names, rest) in operator_lines(t)
        if a.name in names
            group = length(names) > 1 ? Any[" (", Markdown.Code(join(names, ", ")), ")"] : Any[]
            on = isabstracttype(t) ? Any["on every site type declaring what it is made of"] :
                                     Any["on ", Markdown.Code(string(nameof(t)))]
            push!(lines, Any[Markdown.Paragraph(Any[on..., group..., ": ", rest...])])
        end
    end
    if isempty(lines)
        return nothing
    end
    head = Any[Markdown.Code(a.name), ", an operator of one site"]
    mods = [ string(n) for n in names(@__MODULE__)
             if getfield(@__MODULE__, n) isa Module && n ≠ nameof(@__MODULE__) &&
                name in names(getfield(@__MODULE__, n)) ]
    for (i, m) in enumerate(mods)
        push!(head, i == 1 ? ", exported by " : i == length(mods) ? " and " : ", ")
        push!(head, Markdown.Code(m))
    end
    push!(head, ":")
    return Markdown.MD(Any[Markdown.Paragraph(head), Markdown.List(lines),
                           Markdown.Paragraph(Any["The help of a site type, as ",
                                                  Markdown.Code("?Qubit"),
                                                  ", gives its states and its other operators."])])
end
