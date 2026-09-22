export AbstractSite, mix, dim, Index, string_state, identity_operator, state
export @def_operators, @def_states, @create_site_module, conserve_string

"""
    abstract type AbstractSite

An abstract type which is the super type of all site types.

A site type defines `dim`, possibly `string_state`, and its states and operators through
`@def_states` and `@def_operators`. It may also carry a field `conserve::String` recording
what it conserves, which its constructor fills. That field is optional: declare it only if
your site can have conserved quantities, a site without it conserving nothing.
"""
abstract type AbstractSite end

"""
    dim(::AbstractSite)

return the dimension of the given site
"""
dim(site::AbstractSite) = error("dim not implemented on site $site")

# `nameof` rather than `string(typeof(site))`: the latter prints the module prefix when
# the site module is not in scope, and ITensors silently cuts a tag at 16 characters, so
# every site type ended up tagged "TensorMixedState" depending on what the user imported
"""
    Index(::AbstractSite)

return an ITensor.Index for the given site for pure representations
"""
Index(site::AbstractSite) = Index(dim(site); tags="$(nameof(typeof(site))), Site")

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
        error("cannot declare operator $name for site $(typeof(site)): the name already " *
              "stands for a $(typeof(existing)). Choose another one")
    elseif existing.name ≠ name
        error("cannot declare operator $name for site $(typeof(site)): the name already " *
              "stands for the operator $(existing.name)")
    elseif existing.type ≠ type
        error("operator $name is $(existing.type) for a site already in scope and $type " *
              "for site $(typeof(site)): a shared name must agree on the OpType")
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


################ Conserved quantities ################

"""
    charge_tol

how far the eigenvalues of a conserved quantity may sit from the charges they stand for.
Every quantity the built in sites carry lands exactly on an integer, and the roots of unity
of `Zd` miss the unit circle by at most five `eps`, so this is a rounding tolerance and
nothing wider. It is not a setting: a quantity that misses it by more is genuinely
approximate, and the answer is to define it exactly rather than to let it through.
"""
const charge_tol = 1e-14

"""
    site_charges(op, site)

the modulus and the charge each basis state of `site` carries for the conserved quantity
`op`, as `(modulus, charges)`. A modulus of `1` is an ordinary additive charge over the
integers, a modulus of `m` a charge of the cyclic group of order `m`.

A conserved quantity has to be diagonal in the basis the site is written in, and its
eigenvalues have to be readable as charges, which leaves exactly two cases. Integer
eigenvalues are the charge itself. Eigenvalues on the unit circle are roots of unity, the
charge is the exponent and the modulus is read from the denominators — which is what makes
`Zd` work on a `Qudit` with no modulus written anywhere.

The two cases overlap on ±1, which is as much a pair of integers as a pair of square roots
of unity, and they are different conservations: two sites carrying -1 make -2 over the
integers and 0 modulo 2. The integer reading wins, and `parity` is how the other one is
asked for.
"""
function site_charges(op, site::AbstractSite; tol::Float64 = charge_tol)
    m = matrix(op, site)
    d = diag(m)
    off = norm(m - Diagonal(d))
    if off > tol
        error("$op is not diagonal on site $(typeof(site)), off by $(short(off)), so it " *
              "cannot be a conserved quantity: a charge is carried by each basis state")
    end
    to_int = maximum(max(abs(imag(x)), abs(real(x) - round(real(x)))) for x in d)
    if to_int ≤ tol
        return (1, Int.(round.(real.(d))))
    end
    to_circle = maximum(abs(abs(x) - 1) for x in d)
    if to_circle ≤ tol
        θ = angle.(d) ./ (2π)
        modulus = reduce(lcm, denominator.(rationalize.(Int, θ; tol = 1e-8)))
        q = Int.(round.(θ .* modulus))
        to_root = maximum(abs.([exp(2im * π * k / modulus) for k in q] .- d))
        if to_root > tol
            error("the eigenvalues of $op on site $(typeof(site)) are on the unit circle " *
                  "but miss the roots of unity by $(short(to_root)), so they are not charges")
        end
        return (modulus, mod.(q, modulus))
    end
    error("the eigenvalues of $op on site $(typeof(site)) miss the integers by " *
          "$(short(to_int)) and the unit circle by $(short(to_circle)), so they are not " *
          "charges. Half integer ones are written doubled, 2Sz rather than Sz")
end

"""
    short(x)

a number as an error message shows it, two significant digits being all one reads of a
deviation
"""
short(x::Real) = round(x; sigdigits = 2)

# the modulus of a ModOp is carried rather than read back, ±1 being unreadable, and the
# charges are those of its argument taken modulo it
function site_charges(a::ModOp, site::AbstractSite; tol::Float64 = charge_tol)
    m, q = site_charges(a.arg, site; tol)
    if m ≠ 1
        error("cannot take $(a.arg) modulo $(a.modulus) on site $(typeof(site)): it " *
              "already carries a charge modulo $m")
    end
    return (a.modulus, mod.(q, a.modulus))
end

"""
    conserve_string(site, spec)

the form in which a site records what it conserves: for each quantity, its name, its
modulus when that is not 1, and the charge of every basis state.

`spec` is what the user wrote, one operator or a tuple of them, and it is read here rather
than kept, because an operator cannot be written to a state file: the name a conserved
quantity prints under is an expression, `2Sz` or `parity(N)`, and not a key of the operator
library, so it could not be looked up again. What the operator is needed for is the charges,
and those are what travel.

This is what a site type of your own calls to fill its `conserve` field, the site being
built bare first since the charges depend on its type and not on that field.

# Examples

    MySite(; conserve = ()) = MySite(conserve_string(MySite(""), conserve))

    conserve_string(Fermion(""), N)              # "N:0,1"
    conserve_string(Fermion(""), parity(N))      # "parity(N)%2:0,1"
    conserve_string(Electron(""), (Ntot, 2Sz))   # "Ntot:0,1,1,2;2Sz:0,1,-1,0"
"""
function conserve_string(site::AbstractSite, spec)
    ops = spec isa Tuple ? collect(spec) : [spec]
    if isempty(ops)
        return ""
    end
    parts = map(ops) do op
        modulus, q = site_charges(op, site)
        head = modulus == 1 ? obs_name(op) : "$(obs_name(op))%$modulus"
        return head * ":" * join(q, ",")
    end
    return join(parts, ";")
end

"""
    decode_conserve(s)

the conserved quantities a site records, as a vector of `(name, modulus, charges)`.
See `conserve_string`.
"""
function decode_conserve(s::AbstractString)
    if isempty(s)
        return Tuple{String, Int, Vector{Int}}[]
    end
    map(split(s, ';')) do part
        i = findfirst(==(':'), part)
        if isnothing(i)
            error("a site records \"$part\" as a conserved quantity, which has no charges")
        end
        head, tail = part[1:i-1], part[i+1:end]
        j = findlast(==('%'), head)
        name, modulus = isnothing(j) ? (String(head), 1) :
                        (String(head[1:j-1]), parse(Int, head[j+1:end]))
        return (name, modulus, parse.(Int, split(tail, ',')))
    end
end

"""
    conserved(site)

what a site conserves, in the form `conserve_string` produces, and the empty string when it
conserves nothing.

The `conserve` field is optional, a site type that can have no conserved quantity simply not
declaring it, so this is what everything else reads rather than the field itself.
"""
conserved(site::AbstractSite) =
    hasfield(typeof(site), :conserve) ? site.conserve : ""

"""
    conserve_names(s)

the names of the conserved quantities a site records, as they are written back when the site
is printed. A string it cannot read is given back as it is, `show` having to print something
whatever a site put in its field.
"""
function conserve_names(s::AbstractString)
    ns = String[]
    for part in split(s, ';')
        i = findfirst(==(':'), part)
        if isnothing(i)
            return String(s)
        end
        head = part[1:i-1]
        j = findlast(==('%'), head)
        push!(ns, isnothing(j) ? String(head) : String(head[1:j-1]))
    end
    return length(ns) == 1 ? ns[1] : "(" * join(ns, ", ") * ")"
end

"""
    show(io, ::AbstractSite)

print a site as the call that builds it, leaving out the trailing fields that carry nothing,
that is an empty string or `nothing`.

A site conserving nothing therefore goes on printing as it always did, `Qubit()` rather than
`Qubit("")`. Only the trailing ones are left out: dropping a field in the middle would print
a call whose arguments no longer line up with the fields, `Site(nothing, 3)` coming out as
`Site(3)` and reading as something else. The conserved quantities print under their names,
`Fermion(conserve = N)`, the charges they record being an implementation detail.
"""
function show(io::IO, site::AbstractSite)
    t = typeof(site)
    c = conserved(site)
    args = String[]
    nothings = Bool[]
    for f in fieldnames(t)
        v = getfield(site, f)
        if f === :conserve && !isempty(c)
            push!(args, "conserve = " * conserve_names(c))
            push!(nothings, false)
        else
            push!(args, repr(v))
            push!(nothings, isnothing(v) || (v isa AbstractString && isempty(v)))
        end
    end
    n = length(args)
    while n > 0 && nothings[n]
        n -= 1
    end
    print(io, nameof(t), "(", join(args[1:n], ", "), ")")
end
