export AbstractSite, mix, site_dim, site_index, mix, @def_operator, @def_operators

abstract type AbstractSite end

site_dim(site::AbstractSite) = error("site_dim not defined on site $site")
site_index(site::AbstractSite) = Index(site_dim(site); tags="$(string(typeof(site))), Site")

mix(i::Index) =
    addtags(combinedind(combiner(i, i'; tags = tags(i))), "Mixed")

operator_library::Dict{Tuple{DataType, String}, Tuple{Union{Nothing, Matrix, Function, ExprOp}, String}} = Dict()

function operator_info(site::AbstractSite, op::String)
    name = typeof(site)
    t = (name, op)
    if haskey(operator_library, t)
        return operator_library[t]
    else
        error("operator $op is not defined for site $name")
    end
end


operator_doc(site::AbstractSite, op::String) =
    last(operator_info(site, op))

identity_operator(site::AbstractSite) = [ (i==j) ? 1. : 0. for i in 1:site_dim(site), j in 1:site_dim(site) ]

function add_operator(site::AbstractSite, op::String, r::Union{Nothing, Matrix, Function},
    fermionic::Bool=false, doc::String="", N::Int=1)

    name = typeof(site)
    t = (name, op)
    if haskey(operator_library, t)
        error("operator $op is already defined for site $name")
    else
        operator_library[t] = (r, doc)
        return Operator{N}(op, nothing, fermionic)
    end
end

function add_operator(site::AbstractSite, op::String, r::ExprOp{N}, fermionic::Bool=false, doc::String="", _::Int=1) where N
    name = typeof(site)
    t = (name, op)
    if haskey(operator_library, t)
        error("operator $op is already defined for site $name")
    else
        operator_library[t] = (r, doc)
        return Operator{N}(op, nothing, fermionic)
    end
end

macro def_operators(site, symbols, fermionic=false)
    e = Expr(:block)
    if !(symbols isa Expr) || symbols.head ≠ :vect
        error("syntax error in @def_operators second argument should be a vector")
    end
    for p in symbols.args
        if !(p isa Expr) || p.head ≠ :tuple || !(2 <= length(p.args) <= 3)
            error("syntax error in @def_operators items should be 2-tuple")
        end
        expr = first(p.args)
        doc = last(p.args)
        n = length(p.args) == 3 ? p.args[2] : 1
        if !(expr isa Expr) || expr.head ≠ :(=)
            error("syntax error in @def_operators item expression must be an assignment")
        end
        sym = first(expr.args)
        nsym = string(sym)
        val = last(expr.args)
        push!(e.args,
            quote
                $(esc(sym)) = add_operator($(esc(site)), $nsym, $(esc(val)), $(esc(fermionic)), $(esc(doc)), $n)
            end)
    end
    return e
end

macro def_operator(site, expr, fermionic=false, doc="", N=1)
    if !(expr isa Expr) || expr.head ≠ :(=)
        error("syntax error in @def_operator expression must be an assignment")
    end
    sym = first(expr.args)
    nsym = string(sym)
    val = last(expr.args)
    quote
        $(esc(sym)) = add_operator($(esc(site)), $nsym, $(esc(val)), $(esc(fermionic)), $(esc(doc)), $(esc(N)))
    end
end
