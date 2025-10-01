export AbstractSite, mix, dim, Index, string_state, identity_operator, @def_operators, @def_states, state

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
function string_state(site::AbstractSite, st::String)
    i = 1 + parse(Int, st)
    v = zeros(Float64, dim(site))
    v[i] = 1.0
    return v
end

"""
    mix(::Index)

return an ITensor.Index for a mixed representation corresponding to the pure representation Index given
"""
mix(i::Index) =
    addtags(combinedind(combiner(i, i'; tags = tags(i))), "Mixed")

operator_library::Dict{Tuple{DataType, String}, Union{Matrix, Function, GenericOp}} = Dict()
state_library::Dict{Tuple{DataType, String}, Union{String, Vector, Matrix, Function}} = Dict()

function F_info(site::AbstractSite)
    name = typeof(site)
    t = (name, "F")
    return get!(operator_library, t, Id)
end

function operator_info(site::AbstractSite, op::String)
    name = typeof(site)
    t = (name, op)
    r = getkey(operator_library, t, nothing)
    if isnothing(r)
        error("operator $op is not defined for site $name")
    else
        return r
    end
end

function state_info(site::AbstractSite, st::String)
    name = typeof(site)
    t = (name, st)
    r = getkey(state_library, t, nothing)
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
identity_operator(dim::Int) = [ (i==j) ? 1. : 0. for i in 1:dim, j in 1:dim ]
identity_operator(site::AbstractSite) = identity_operator(dim(site))

function add_operator(site::AbstractSite, op::String, r::Union{Matrix, Function}, type::OpType=plain_op, N::Int=1)
    name = typeof(site)
    t = (name, op)
    if haskey(operator_library, t)
        error("operator $op is already defined for site $name")
    else
        operator_library[t] = r
        return Operator{N}(op, nothing, type)
    end
end

function add_operator(site::AbstractSite, op::String, r::GenericOp{Pure, N}, type::OpType=plain_op, ::Int=1) where N
    name = typeof(site)
    t = (name, op)
    if haskey(operator_library, t)
        error("operator $op is already defined for site $name")
    else
        operator_library[t] = r
        return Operator{N}(op, nothing, type)
    end
end

"""
    @def_operators(site, symbols, type=plain)

define the given operator for the given site

# Examples

    @def_operators(Fermion(),
    [
        C = [0. 1. ; 0. 0.],
    ],
    fermionic_op)

    @def_operators(Fermion(),
    [
        F = [1. 0. ; 0. -1.],
        A = C,
        N = dag(C)*C
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
                push!(e.args,
                quote
                    add_operator($(esc(site)), $nsym, $(esc(val)), $(esc(type)))
                end)
            else
                push!(e.args,
                    quote
                        $(esc(sym)) = add_operator($(esc(site)), $nsym, $(esc(val)), $(esc(type)))
                    end)
            end
        end
    end
    return e
end

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

state(::AbstractSite, a::Union{Vector, Matrix}) = a
state(site::AbstractSite, a::Function) = a(site)

"""
    state(::AbstractSite, ::String)

return the local state (as a vector or matrix) corresponding to the site and name given

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
