export AbstractSite, mix, dim, Index, string_state, identity_operator, state
export @def_operators, @def_states, @create_site_module

"""
    abstract type AbstractSite

An abstract type which is the super type of all site types
"""
abstract type AbstractSite end

"""
    dim(::AbstractSite)

return the dimension of the given site
"""
dim(site::AbstractSite) = error("dim not implemented on site $site")

"""
    Index(::AbstractSite)

return an ITensor.Index for the given site for pure representations
"""
Index(site::AbstractSite) = Index(dim(site); tags="$(string(typeof(site))), Site")

"""
    string_state(::AbstractSite, ::String)

Do not call directly. It returns a local state corresponding to the string,
this is tried first before trying specifically defined states.

The default implementation returns the first state for "0", the second for "1" and so on.

This should be overloaded if necessary when defining new site types. It should return an error when not needed.
"""
string_state(site::AbstractSite, st::String) =
    state(site, parse(Int, st))

"""
    mix(::Index)

return an ITensor.Index for a mixed representation corresponding to the pure representation Index given
"""
mix(i::Index) =
    addtags(combinedind(combiner(i, i'; tags = tags(i))), "Mixed")

"""
    operator_library::Dict

a global variable containing the site dependent definitions
of implicit operators as defined by `@def_operators`
"""
const operator_library::Dict{Tuple{DataType, String}, Union{Matrix, Function, GenericOp}} = Dict()

"""
    state_library::Dict

a global variable containing the definitions of local states as defined by `@def_states`
"""
const state_library::Dict{Tuple{DataType, String}, Union{String, Vector, Matrix, Function}} = Dict()

"""
    F_info(site)

return the matrix value of `F`` for the `site`` as stored in `operator_library`

`get!` here is a cache rather than a lookup: a site with no `F` of its own is not
fermionic, and the `Id` it falls back to is written into the library so that the next call
finds it.
"""
function F_info(site::AbstractSite)
    name = typeof(site)
    t = (name, "F")
    return get!(operator_library, t, Id)
end

"""
    operator_info(site, op)

return the definition of `op` for the given `site` as stored in `operator_library`,
it may be an `Op` a matrix or a site function
"""
function operator_info(site::AbstractSite, op::String)
    name = typeof(site)
    t = (name, op)
    r = get(operator_library, t, nothing)
    if isnothing(r)
        error("operator $op is not defined for site $name")
    else
        return r
    end
end

"""
    state_info(site, statename)

return the state definition for `site` as stored in `state_library`,
it may be a `Vector` (for pure state), a `Matrix` for mixed states or a site function
"""
function state_info(site::AbstractSite, st::String)
    name = typeof(site)
    t = (name, st)
    r = get(state_library, t, nothing)
    if isnothing(r)
        error("state $st is not defined for site $name")
    else
        return r
    end
end

"""
    identity_operator(::AbstractSite)

return a matrix representing the identity operator for the given site
"""
identity_operator(dim::Int) = Matrix{Float64}(I, dim, dim)
identity_operator(site::AbstractSite) = identity_operator(dim(site))

"""
    add_operator(site, op, r, type = plain_op)

register the definition `r` of the operator named `op` for the given `site`, and return the
`Operator{1}` standing for that name

Do not call directly, use `@def_operators`, which is what keeps the operator name, its
`OpType` and the definitions made for the other site types consistent.
"""
function add_operator(site::AbstractSite, op::String, r::Union{Matrix, Function, GenericOp{Pure, 1}}, type::OpType=plain_op)
    name = typeof(site)
    t = (name, op)
    if haskey(operator_library, t)
        error("operator $op is already defined for site $name")
    else
        operator_library[t] = r
        return Operator{1}(op, nothing, type)
    end
end

"""
    check_shared_operator(existing, name, type, site)

check that a name already in scope can stand for the operator about to be registered
"""
function check_shared_operator(existing, name::String, type::OpType, site::AbstractSite)
    if !(existing isa Operator{1})
        error("cannot declare operator $name for site $(typeof(site)): the name $name is " *
              "already used in this module for something else, of type $(typeof(existing)). " *
              "Give your operator another name, or make sure $name is not brought into scope")
    elseif existing.name ≠ name
        error("cannot declare operator $name for site $(typeof(site)): in this module the " *
              "name $name already stands for the operator $(existing.name)")
    elseif existing.type ≠ type
        error("operator $name is declared as $(existing.type) by a site type already in " *
              "scope, and as $type for site $(typeof(site)). A name stands for one " *
              "operator, shared by every site type that declares it, so the two " *
              "declarations must agree on the OpType")
    end
    return existing
end

"""
    @def_operators(site, symbols)

define the given operators for the given site, see also `OpType`

Each operator name becomes a `const` of the module the macro is called from, but only the
first time that name is seen: a name already in scope is registered for the new site and
checked against what it already stands for, not bound again. Declaring an operator whose name is already used for something else, or declared
with another `OpType`, is an error rather than a silent redefinition.

# Examples

    @def_operators(Fermion(),
    [
        fermionic_op => 
        [
            C = [0. 1. ; 0. 0.],
        ],
        selfadjoint_op =>
        [
            N = dag(C) * C,
        ],
        plain_op =>
        [
            A = C,
        ],
        involution_op =>
        [
            F = Float64[1 0 ; 0 -1]
        ]
    ])
"""
macro def_operators(site, symbols)
    e = Expr(:block)
    if !(symbols isa Expr) || symbols.head ≠ :vect
        error("syntax error in @def_operators second argument should be a vector")
    end
    for types in symbols.args
        if !(types isa Expr) || types.head ≠ :call || types.args[1] ≠ :(=>)
            error("syntax error in @def_operators second argument should contain pairs : plain_op => [...]")
        end

        type = types.args[2]
        for expr in types.args[3].args
            if !(expr isa Expr) || expr.head ≠ :(=)
                error("syntax error in @def_operators item expressions must be assignments (sym = val)")
            end
            sym = first(expr.args)
            nsym = string(sym)
            val = last(expr.args)
            if nsym == "F"
                # `F` is the Jordan-Wigner operator of `Operators.jl`, shared by every
                # fermionic site and not an `Operator{1}`: the site is registered and the
                # name is left alone
                push!(e.args,
                quote
                    add_operator($(esc(site)), $nsym, $(esc(val)), $(esc(type)))
                end)
            elseif isdefined(__module__, sym)
                # the name is already in scope, so it is registered for this site and
                # checked, but not bound again. Binding it again would rebind it for every
                # site already using it, and up to Julia 1.11 rebinding a name brought in by
                # `using` is a hard error of the language. The decision is taken here, at
                # expansion time, so that no binding is emitted at all in that case
                push!(e.args,
                    quote
                        check_shared_operator($(esc(sym)), $nsym, $(esc(type)), $(esc(site)))
                        add_operator($(esc(site)), $nsym, $(esc(val)), $(esc(type)))
                    end)
            else
                push!(e.args,
                    quote
                        const $(esc(sym)) = add_operator($(esc(site)), $nsym, $(esc(val)), $(esc(type)))
                    end)
            end
        end
    end
    return e
end

"""
    add_state(site, st, r)
    add_state(site, sts, r)

register the definition `r` of the state named `st` for the given `site`, or the same
definition for every name of the vector `sts`

Do not call directly, use `@def_states`.
"""
function add_state(site::AbstractSite, st::String, r::Union{String, Vector, Matrix, Function})
    name = typeof(site)
    t = (name, st)
    if haskey(state_library, t)
        error("state $st is already defined for site $name")
    else
        state_library[t] = r
    end
end

add_state(site::AbstractSite, sts::Vector{String}, r::Union{String, Vector, Matrix, Function}) =
    foreach(sts) do st
        add_state(site, st, r)
    end

"""
    @def_states(site, symbols)

define the given states for the given site

# Example

    @def_states(Fermion(),
    [
        ["Emp", "0"] => [1., 0.],
        "Occ" => [0., 1.],
    ])

"""
macro def_states(site, symbols)
    e = Expr(:block)
    if !(symbols isa Expr) || symbols.head ≠ :vect
        error("syntax error in @def_states second argument should be a vector")
    end
    for expr in symbols.args
        if !(expr isa Expr) || expr.head ≠ :call || expr.args[1] ≠ :(=>)
            error("syntax error in @def_states item expressions must be pairs (\"sym\" => val or [\"sym1\", \"sym2\", ...] => val)")
        end
        sym = expr.args[2]
        val = expr.args[3]
        push!(e.args,
            quote
                add_state($(esc(site)), $(esc(sym)), $(esc(val)))
            end)
    end
    return e
end

"""
    @create_site_module(name, symbols)

define a submodule named `name` which imports and re-exports the given symbols from the parent module

# Example

    @create_site_module(Spins, [Spin, Sp, Sm, Sx, Sy, Sz, S2])
"""
macro create_site_module(name, symbols)
    if !(symbols isa Expr) || symbols.head ≠ :vect
        error("syntax error in @create_site_module second argument should be a vector of symbols")
    end
    imports = [ Expr(:., :., :., s) for s in symbols.args ]
    block = Expr(:block,
        Expr(:import, imports...),
        Expr(:export, symbols.args...))
    doc = "    using .$name\n\nmodule giving access to the `$(symbols.args[1])` site type " *
        "and its operators: " * join(map(s -> "`$s`", symbols.args[2:end]), ", ")
    mod = Expr(:module, true, name, block)
    return esc(Expr(:macrocall, GlobalRef(Core, Symbol("@doc")), __source__, doc, mod))
end

state(::AbstractSite, a::Union{Vector, Matrix}) = a
state(site::AbstractSite, a::Function) = a(site)
function state(site::AbstractSite, a::Int)
    if a < 0 || a >= dim(site)
        error("invalid state number")
    end
    v = fill(0., dim(site))
    v[a + 1] = 1.0
    return v
end

"""
    state(::AbstractSite, ::String)

return the local state (as a vector or matrix) corresponding to the site and name given
the special name "FullyMixed" gives the infinite temperature state

# Examples

```jldoctest
julia> using TensorMixedStates, .Qubits, .Fermions

julia> state(Qubit(), "+")
2-element Vector{Float64}:
 0.7071067811865475
 0.7071067811865475

julia> state(Fermion(), "FullyMixed")
2×2 Matrix{Float64}:
 0.5  0.0
 0.0  0.5
```
"""
state(site::AbstractSite, st::String) =
    if st == "FullyMixed"
        identity_operator(site) / dim(site)
    else
        try
            string_state(site, st)
        catch
            state(site, state_info(site, st))
        end
    end
