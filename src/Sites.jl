# Site types: what every site provides, its dimension, the parity of its basis states, and the
# states and operators declared for its type, which become ordinary methods.

export AbstractSite, dim, fermion_parities, fermion_species, string_state
export Level, @def_states, @def_operators, @create_site_module

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
    declaration_names(name, what)

the names a declaration gives, from a literal string or a vector of them, refused otherwise:
a declaration becomes methods when the macro expands
"""
function declaration_names(name, what::String)
    if name isa String
        return [name]
    elseif name isa Expr && name.head == :vect && all(x -> x isa String, name.args)
        return String.(name.args)
    end
    error("@$what names a state by a literal string or a vector of them, not $name")
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
        for name in declaration_names(d.args[2], "def_states")
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
