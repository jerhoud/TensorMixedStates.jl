export AbstractSite, mix, dim, Index, identity_operator, @def_operator, @def_operators, @def_states, state, Id, F

abstract type AbstractSite end

dim(site::AbstractSite) = error("dim not defined on site $site")
Index(site::AbstractSite) = Index(dim(site); tags="$(string(typeof(site))), Site")
function generic_state(::AbstractSite, st::String)
    i = 1 + parse(Int, st)
    v = zeros(Float64, dim(site))
    v[i] = 1.0
    return v
end

mix(i::Index) =
    addtags(combinedind(combiner(i, i'; tags = tags(i))), "Mixed")

operator_library::Dict{Tuple{DataType, String}, Union{Matrix, Function, ExprOp}} = Dict()
state_library::Dict{Tuple{DataType, String}, Union{String, Vector, Matrix, Function}} = Dict()

function operator_info(site::AbstractSite, op::String)
    name = typeof(site)
    t = (name, op)
    if haskey(operator_library, t)
        return operator_library[t]
    else
        error("operator $op is not defined for site $name")
    end
end

function state_info(site::AbstractSite, st::String)
    name = typeof(site)
    t = (name, st)
    if haskey(state_library, t)
        return state_library[t]
    else
        error("state $st is not defined for site $name")
    end
end

identity_operator(dim::Int) = [ (i==j) ? 1. : 0. for i in 1:dim, j in 1:dim ]
identity_operator(site::AbstractSite) = identity_operator(dim(site))

function add_operator(site::AbstractSite, op::String, r::Union{Matrix, Function}, fermionic::Bool=false, N::Int=1)
    name = typeof(site)
    t = (name, op)
    if haskey(operator_library, t)
        error("operator $op is already defined for site $name")
    else
        operator_library[t] = r
        return Operator{Pure, N}(op, nothing, fermionic)
    end
end

function add_operator(site::AbstractSite, op::String, r::ExprOp{T, N}, fermionic::Bool=false, ::Int=1) where {T, N}
    name = typeof(site)
    t = (name, op)
    if haskey(operator_library, t)
        error("operator $op is already defined for site $name")
    else
        operator_library[t] = r
        return Operator{T, N}(op, nothing, fermionic)
    end
end

macro def_operators(site, symbols, fermionic=false)
    e = Expr(:block)
    if !(symbols isa Expr) || symbols.head ≠ :vect
        error("syntax error in @def_operators second argument should be a vector")
    end
    for expr in symbols.args
        if !(expr isa Expr) || expr.head ≠ :(=)
            error("syntax error in @def_operators item expressions must be assignments (sym = val)")
        end
        sym = first(expr.args)
        nsym = string(sym)
        val = last(expr.args)
        push!(e.args,
            quote
                $(esc(sym)) = add_operator($(esc(site)), $nsym, $(esc(val)), $(esc(fermionic)))
            end)
    end
    return e
end

macro def_operator(site, expr, N=1, fermionic=false)
    if !(expr isa Expr) || expr.head ≠ :(=)
        error("syntax error in @def_operator expression must be an assignment (sym = val)")
    end
    sym = first(expr.args)
    nsym = string(sym)
    val = last(expr.args)
    quote
        $(esc(sym)) = add_operator($(esc(site)), $nsym, $(esc(val)), $(esc(fermionic)), $(esc(N)))
    end
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

state(site::AbstractSite, st::String) =
    if st == "FullyMixed"
        identity_operator(site) / dim(site)
    else
        try
            generic_state(site, st)
        catch
            state(site, state_info(site, st))
        end
    end

Id = Operator("Id", identity_operator)
F = Operator("F")