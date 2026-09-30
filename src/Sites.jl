# Site types: what every site provides (its dimension, its ITensor index, its local states and
# the library of its operators), and the conserved quantities a site records, strong or weak,
# with the relabellings of the charges they need.

export AbstractSite, dim, Index, string_state, identity_operator, state, weaken,
       symmetries
export @def_states, @create_site_module, strong

"""
    abstract type AbstractSite

the supertype of all site types.

A site type defines `dim`, possibly `string_state`, and its states and operators with
`@def_states` and `@def_operators`. It may carry a field `conserve::String`, which its
constructor fills through `conserve_string` with what the site conserves. The field is
optional: a site type that can conserve nothing leaves it out.

# Examples

    struct MySite <: AbstractSite
        conserve::String
    end

    MySite(; conserve = ()) = MySite(conserve_string(MySite(""), conserve))
    TensorMixedStates.dim(::MySite) = 2
"""
abstract type AbstractSite end

"""
    dim(site)

the dimension of the Hilbert space of a site, which every site type must define.

# Examples

    dim(Spin(1))    # 3
"""
dim(site::AbstractSite) = error("dim not implemented on site $site")

"""
    conserved(site)

what a site conserves, as `conserve_string` records it, or `""` when it conserves nothing.
Everything reads this rather than the optional field `conserve`: a site type without that
field, or with one that is not a string, conserves nothing.
"""
conserved(site::AbstractSite) =
    if hasfield(typeof(site), :conserve) && getfield(site, :conserve) isa AbstractString
        site.conserve
    else
        ""
    end

"""
    decode_conserve(s)

the quantities recorded in the string `s` of `conserve_string`, as a vector of
`(name, modulus, charges, strong)`: `charges` holds the charge of each basis state, `modulus`
is 1 for a charge taken modulo nothing, and `strong` tells whether it is conserved strongly.
"""
function decode_conserve(s::AbstractString)
    if isempty(s)
        return Tuple{String, Int, Vector{Int}, Bool}[]
    end
    map(split(s, ';')) do part
        i = findfirst(==(':'), part)
        if isnothing(i)
            error("a site records \"$part\" as a conserved quantity, which has no charges")
        end
        head, tail = part[1:i-1], part[i+1:end]
        st = endswith(head, '!')
        if st
            head = head[1:end-1]
        end
        j = findlast(==('%'), head)
        name, modulus = isnothing(j) ? (String(head), 1) :
                        (String(head[1:j-1]), parse(Int, head[j+1:end]))
        return (name, modulus, parse.(Int, split(tail, ',')), st)
    end
end

"""
    sorted_conserve(s)

the string `s` of `conserve_string` with its quantities sorted by name, as `Conserved` holds
them. A state file written before they were sorted holds them in the order of their
declaration.
"""
function sorted_conserve(s::AbstractString)
    if isempty(s)
        return String(s)
    end
    return join(sort(split(s, ';'); by = p -> first(only(decode_conserve(p)))), ";")
end

"""
    strong_names(site)

the names of the quantities the site conserves strongly. See `strong`.
"""
strong_names(site::AbstractSite) =
    [ q[1] for q in decode_conserve(conserved(site)) if q[4] ]

"""
    basis_charges(site)

the charge of each basis state of `site`, as a `QN` holding every quantity it conserves, and
`QN()` for each when it conserves nothing
"""
function basis_charges(site::AbstractSite)
    n = dim(site)
    qs = decode_conserve(conserved(site))
    for (name, _, charges, _) in qs
        if length(charges) ≠ n
            error("site $(typeof(site)) records $(length(charges)) charges for $name but " *
                  "has $n basis states")
        end
    end
    return [ QN([ (name, charges[k], modulus) for (name, modulus, charges, _) in qs ]...)
             for k in 1:n ]
end

"""
    site_index(site, charged)

the ITensor index of a site in the pure representation, `charged` telling whether its system
carries quantum numbers.

A site conserving nothing in a charged system takes a trivial index, one sector of charge zero
holding the whole space, since an MPS cannot mix indices with and without charges; every
operator keeps its matrix, in that one block. The sectors come one per basis state rather than
merged by charge: merging would reorder the basis whenever equal charges are not contiguous,
as those of `parity(N)` on a boson, 0, 1, 0, 1.
"""
function site_index(site::AbstractSite, charged::Bool)
    n = dim(site)
    # `nameof` rather than `string(typeof(site))`: the latter prints the module prefix when
    # the site module is not in scope, and ITensors silently cuts a tag at 16 characters, so
    # every site type ended up tagged "TensorMixedState" depending on what the user imported
    tg = "$(nameof(typeof(site))), Site"
    if isempty(conserved(site))
        return charged ? Index(QN() => n; tags = tg) : Index(n; tags = tg)
    end
    return Index([ q => 1 for q in basis_charges(site) ]...; tags = tg)
end

"""
    qn_components(q::QN)

the components of the charge `q`, as a vector of `(name, value, modulus)`, a modulus of 1
meaning an integer charge. A `QN` always has four slots, the unused ones having an empty
name: they are left out. `QN(components...)` gives the charge back, which is how `weak_qn`,
`star` and `adjoint_qn` rebuild a charge after changing its components.
"""
function qn_components(q::QN)
    cs = Tuple{String, Int, Int}[]
    for v in q.data
        n = String(ITensors.name(v))
        if !isempty(n)
            push!(cs, (n, ITensors.val(v), ITensors.modulus(v)))
        end
    end
    return cs
end

"""
    map_charges(f, i)

the index `i` with the charge `q` of each block replaced by `f(q)`, the dimensions of the
blocks, the tags, the prime level and the direction kept.
"""
map_charges(f, i::Index) =
    Index([ f(q) => d for (q, d) in space(i) ]...; tags = tags(i), plev = plev(i), dir = dir(i))

"""
    weak_qn(q, collapse, drop)
    weak_index(i, collapse, drop)

the charge, or the index, with the quantities of `drop` left out and the strong quantities of
`collapse` made weak: their components `X` and `X*`, for the ket and the daggered bra, are
summed into `X`, which gives the difference of the two charges, what a weak symmetry records.

Both operations are homomorphisms of the charge group: they map a flux to a flux, so every
tensor on a relabelled index stays consistent, with no data moved. The blocks are not merged,
several may share a charge, which keeps the relabelling free.
"""
function weak_qn(q::QN, collapse, drop)
    vals = Tuple{String, Int, Int}[]
    for (nm, v, m) in qn_components(q)
        base = endswith(nm, "*") ? nm[1:end-1] : nm
        if base in drop
            continue
        end
        if !(base in collapse)
            push!(vals, (nm, v, m))
            continue
        end
        k = findfirst(x -> x[1] == base, vals)
        if isnothing(k)
            push!(vals, (base, v, m))
        else
            vals[k] = (base, vals[k][2] + v, vals[k][3])
        end
    end
    return QN(vals...)
end

weak_index(i::Index, collapse, drop) =
    if (isempty(collapse) && isempty(drop)) || !hasqns(i)
        i
    else
        map_charges(q -> weak_qn(q, collapse, drop), i)
    end

"""
    star(q::QN, names)
    star(i::Index, names)

the charge, or the index, with each component `X` named in `names` renamed `X*`.

Under a strong symmetry the bra carries its charges under the starred names, so that pairing
it with the ket keeps both charges rather than their difference. The links of a state carry
the same charges and are renamed with the sites. With no name, the index is unchanged.
"""
star(q::QN, names) =
    QN([ (n in names ? n * "*" : n, v, m) for (n, v, m) in qn_components(q) ]...)

function star(i::Index, names)
    if isempty(names) || !hasqns(i)
        return i
    end
    return map_charges(q -> star(q, names), i)
end

"""
    adjoint_qn(q::QN, names)
    adjoint_index(i::Index, names)

the charge, or the index, relabelled for the adjoint: the element ``|x\\rangle\\langle y|``
of a mixed index takes the charge ``|y\\rangle\\langle x|`` had, `names` being the quantities
conserved strongly.

A strong quantity carries `X`, the charge of `x`, and `X*`, minus that of `y`, so `(X, X*)`
becomes `(-X*, -X)`; a weak one holds their difference and is only negated. This is an
automorphism of the charge group: relabelling every index of a state, links included, keeps
each tensor consistent with no data moved, and what remains of the adjoint is a permutation of
zero flux, see `adj_map`.
"""
adjoint_qn(q::QN, names) =
    QN([ (endswith(n, "*") ? n[1:end-1] : n in names ? n * "*" : n, -v, m)
         for (n, v, m) in qn_components(q) ]...)

adjoint_index(i::Index, names) =
    if !hasqns(i)
        i
    else
        map_charges(q -> adjoint_qn(q, names), i)
    end

"""
    bra_index(i, site)

the index carrying the bra of `i`: `i` with the quantities `site` conserves strongly starred,
which is `i` itself when there are none. See `strong` and `star`.
"""
bra_index(i::Index, site::AbstractSite) = star(i, strong_names(site))

"""
    Index(site)

the ITensor index of a site in the pure representation, carrying the quantum numbers of what
the site conserves, and none when it conserves nothing.

In a system where another site conserves something, a site conserving nothing takes instead
an index of a single sector of charge zero, which depends on the system and not on the site.
"""
Index(site::AbstractSite) = site_index(site, !isempty(conserved(site)))

"""
    mixed_index(i, site)

a new index of the mixed representation, combining the ket index `i` with the bra of the same
site.

The bra is daggered, so that the charge of ``|m\\rangle\\langle n|`` is the difference of
those of ``m`` and ``n``: this is the pairing `mix(::State)` produces, contracting `t` with
`dag(t')`. Under a strong symmetry the bra is also starred, which keeps the two charges
apart, see `strong`.

The index is fresh, so it contracts with nothing already built: a caller wants the one its
system drew, `SysIndex{Mixed}(system, i)`.
"""
mixed_index(i::Index, site::AbstractSite) =
    addtags(combinedind(combiner(i, dag(bra_index(i, site)'); tags = tags(i))), "Mixed")

"""
    operator_library

the definitions of the operators declared by `@def_operators`, keyed by site type and
operator name: a matrix, a function of the site returning one, or an operator expression.
"""
const operator_library::Dict{Tuple{DataType, String}, Union{Matrix, Function, GenericOp}} = Dict()

"""
    state_library

the definitions of the states declared by `@def_states`, keyed by site type and state name:
a vector, a density matrix, a function of the site returning either, or the name of another
state.
"""
const state_library::Dict{Tuple{DataType, String}, Union{String, Vector, Matrix, Function}} = Dict()

"""
    F_info(site)

return the matrix value of `F` for the `site` as stored in `operator_library`, the
identity for a site with no `F` of its own, which is not fermionic.
"""
function F_info(site::AbstractSite)
    name = typeof(site)
    t = (name, "F")
    return get(operator_library, t, Id)
end

"""
    operator_info(site, op)

the definition of the operator named `op` for the type of `site`, as `operator_library` holds
it, and an error when there is none.
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

the definition of the state `statename` for the type of `site`, as `state_library` holds it,
and `nothing` when the site declares no state of that name.
"""
state_info(site::AbstractSite, st::String) = get(state_library, (typeof(site), st), nothing)

"""
    identity_operator(site)
    identity_operator(dim::Int)

the identity matrix of a site, or of dimension `dim`, as a `Matrix{Float64}`.
"""
identity_operator(dim::Int) = Matrix{Float64}(I, dim, dim)
identity_operator(site::AbstractSite) = identity_operator(dim(site))

"""
    add_operator(site, op, r, type = plain_op)

register `r` as the definition of the operator named `op` for the type of `site`, refusing a
second one, and return the `Operator{1}` standing for that name.

Only `@def_operators` calls it, which keeps the name, its `OpType` and the definitions made
for the other site types consistent.
"""
function add_operator(site::AbstractSite, op::String, r::Union{Matrix, Function, SimpleOp}, type::OpType=plain_op)
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

check that `existing`, the value already bound to `name`, can stand for the operator about to
be registered for `site`, and return it. It must be an `Operator{1}` of that name and `type`,
with no definition of its own, since such an operator never reads the library of the sites.
"""
function check_shared_operator(existing, name::String, type::OpType, site::AbstractSite)
    if !(existing isa Operator{1})
        error("cannot declare operator $name for site $(typeof(site)): the name already " *
              "stands for a $(typeof(existing)). Choose another one")
    elseif existing.name ≠ name
        error("cannot declare operator $name for site $(typeof(site)): the name already " *
              "stands for the operator $(existing.name)")
    elseif !isnothing(existing.expr)
        # an operator with a definition of its own never reads the library of the sites, so
        # the declaration would be recorded and never used
        error("cannot declare operator $name for site $(typeof(site)): the name already " *
              "stands for an operator with a definition of its own. Choose another one")
    elseif existing.type ≠ type
        error("operator $name is $(existing.type) for a site already in scope and $type " *
              "for site $(typeof(site)): a shared name must agree on the OpType")
    end
    return existing
end

"""
    add_state(site, st, r)
    add_state(site, sts, r)

register `r` as the definition of the state named `st`, or of every name of `sts`, for the
type of `site`, refusing a second one. Only `@def_states` calls it.
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
    @def_states(site, [name => def, ...])

declare states for the type of `site`. A name is a string, or a vector of strings naming the
same state. A definition is a vector (a pure state in the basis of the site), a matrix (a
density matrix), a function of the site returning either, or the name of another state.
Declaring a name twice for a site type is an error.

# Examples

For a fermionic site type of your own, `MySite`, declared as the package declares `Fermion`:

    @def_states(MySite(),
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

define a submodule `name` importing the given symbols from the module the macro is called in
and exporting them, so that `using .name` brings a site type and its operators into scope. The
submodule gets a docstring naming them, the first symbol as the site type.

# Examples

    @create_site_module(Spins, [Spin, Sp, Sm, Sx, Sy, Sz, S2, N])
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
        "and the names it exports: " * join(map(s -> "`$s`", symbols.args[2:end]), ", ")
    mod = Expr(:module, true, name, block)
    return esc(Expr(:macrocall, GlobalRef(Core, Symbol("@doc")), __source__, doc, mod))
end

function state(site::AbstractSite, a::Union{Vector, Matrix})
    d = dim(site)
    if a isa Vector && length(a) ≠ d
        error("a state of $(length(a)) components cannot be one of $site, whose dimension is $d")
    elseif a isa Matrix && size(a) ≠ (d, d)
        error("a $(size(a, 1))×$(size(a, 2)) density matrix cannot be one of $site, whose " *
              "dimension is $d")
    end
    return a
end

state(site::AbstractSite, a::Function) = state(site, a(site))
function state(site::AbstractSite, a::Int)
    if a < 0 || a >= dim(site)
        error("invalid state number")
    end
    v = fill(0., dim(site))
    v[a + 1] = 1.0
    return v
end

"""
    string_state(site, name)

the state a name gives by a rule of the site rather than by a declaration, tried by `state`
after the states the site declares and before those every site has. It is not called
directly.

By default `"0"` gives the first basis state, `"1"` the second, and so on. A site type
overloads it to read names of its own, as `Spin` reads `"1/2"` or `"X1/2"`, or to read none
by raising an error, which `state` takes to mean that the name is not one of its forms.

# Examples

    TensorMixedStates.string_state(::MySite, ::String) = error("no generic state for MySite")
"""
string_state(site::AbstractSite, st::String) =
    state(site, parse(Int, st))

"""
    common_states

the states every site has, as functions of the site: `"FullyMixed"`, the density matrix
proportional to the identity, that is the infinite temperature state. `state` looks them up
last, so that a site may declare its own under the same name.
"""
const common_states = Dict{String, Function}(
    "FullyMixed" => s -> identity_operator(s) / dim(s),
)

"""
    state(::AbstractSite, ::String)

the local state a name gives on a site: a vector, or a density matrix for a mixed one such as
`"FullyMixed"`.

The name is looked for among the states the site declares, see `@def_states`, then among its
generic forms, see `string_state`, then among the states every site has, `"FullyMixed"`.

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
function state(site::AbstractSite, st::String)
    declared = state_info(site, st)
    if !isnothing(declared)
        return state(site, declared)
    end
    generic = try
        string_state(site, st)
    catch e
        # an error is how `string_state` says the name is not one of its forms. An interrupt
        # is not one, and has to go on to where a simulation stops on it
        if e isa InterruptException
            rethrow()
        end
        nothing
    end
    if !isnothing(generic)
        return generic
    end
    common = get(common_states, st, nothing)
    if !isnothing(common)
        return state(site, common)
    end
    error("state $st is not defined for site $(typeof(site))")
end


################ Conserved quantities ################

"""
    charge_tol

how far the eigenvalues of a conserved quantity may lie from the charges they stand for. It is
a rounding tolerance and not a setting: a quantity missing its charges by more is
approximate, and has to be defined exactly rather than let through.
"""
const charge_tol = 1e-14

"""
    rounding_tol

the size, relative to the norm of a matrix, below which an element, a singular value or a
whole term is taken as rounding and treated as zero. The flux of a matrix (`charge_flux`),
its tensor on charged indices (`charged_itensor`) and the splitting of an operator into one
site factors (`Operator{N}(name, def, type, sites...)`) all go by it and have to agree.

It is a rounding tolerance and not a setting: operators are not compressed here, the states
they act on are truncated by the algorithms.
"""
const rounding_tol = 1e-13

"""
    short(x)

a number rounded to two significant digits, as an error message shows a deviation.
"""
short(x::Real) = round(x; sigdigits = 2)

"""
    show_charges(d)

a list of `(name, value, modulus)` charges as an error message shows it, `Ntot=-1,2Sz=-1`.
"""
show_charges(d) = join(["$name=$val" for (name, val, _) in d], ",")

"""
    charge_flux(m, what, site)

the flux of the matrix `m` of `what` on `site`, and an error naming `what` when it has none.
`flux` answers with it and the tensor of an operator is checked with it, so that a matrix not
fitting the charges of its site is refused by name rather than by the `Fluxes not all equal` of
ITensors. Elements below `tol` relative to the norm are rounding and ignored, as
`charged_itensor` does.
"""
function charge_flux(m::Matrix, what, site::AbstractSite; tol::Float64 = rounding_tol)
    qs = decode_conserve(conserved(site))
    if isempty(qs)
        return QN()
    end
    small = tol * norm(m)
    found = nothing
    for i in axes(m, 1), j in axes(m, 2)
        if abs(m[i, j]) ≤ small
            continue
        end
        d = [ (name, modulus == 1 ? ch[i] - ch[j] : mod(ch[i] - ch[j], modulus), modulus)
              for (name, modulus, ch, _) in qs ]
        if isnothing(found)
            found = d
        elseif d ≠ found
            error("$what on site $(typeof(site)) connects charges differing by " *
                  "$(show_charges(found)) and by $(show_charges(d)), so it has no definite flux")
        end
    end
    # an operator with no element at all carries no charge
    return isnothing(found) ? QN() : QN(found...)
end

"""
    charged_itensor(a, inds)

the ITensor of the array `a` on the indices `inds`, dropping what rounding leaves outside the
blocks of charged indices, see `rounding_tol`. Plain indices keep everything.

A matrix computed through an eigendecomposition, as an exponential or a non integer power,
holds elements of the size of rounding between charges the exact one keeps apart. ITensors,
which drops nothing by default, would make a block of each and refuse the tensor for its
fluxes. A tensor that genuinely has no definite charge is still refused.
"""
charged_itensor(a::AbstractArray, inds) =
    if any(hasqns, inds)
        ITensor(a, inds...; tol = rounding_tol * norm(a))
    else
        ITensor(a, inds...)
    end

"""
    index_charges(i)

the charge each value of the index `i` brings to the flux of a tensor: the charge of its
block, negated on an incoming index.
"""
index_charges(i::Index) =
    reduce(vcat, [ fill(dir(i) == ITensors.In ? -q : q, n) for (q, n) in space(i) ]; init = QN[])

"""
    common_charge(a, charges, tol)

the charge every element of the array `a` above `tol` carries, `charges[k][v]` being the
charge the value `v` of axis `k` brings: `QN()` when no element is above `tol`, and `nothing`
when two of them carry different charges
"""
function common_charge(a::AbstractArray, charges, tol)
    found = nothing
    for c in CartesianIndices(a)
        if abs(a[c]) ≤ tol
            continue
        end
        q = sum(charges[k][c[k]] for k in eachindex(charges))
        if isnothing(found)
            found = q
        elseif q ≠ found
            return nothing
        end
    end
    return isnothing(found) ? QN() : found
end

"""
    has_definite_flux(a, inds)

whether the array `a` on the indices `inds` carries a definite charge: every element above
rounding, by the rule of `charged_itensor`, connects states whose charges differ by the same
amount. An array on plain indices always does. It asks beforehand what ITensors answers with
`Fluxes not all equal`, so that a refusal can name what it refuses without catching an error
of ITensors.
"""
function has_definite_flux(a::AbstractArray, inds)
    if !any(hasqns, inds)
        return true
    end
    return !isnothing(common_charge(a, map(index_charges, inds), rounding_tol * norm(a)))
end

"""
    struct Strong

a conserved quantity declared strong, as `strong` returns it, wrapping its operator.
"""
struct Strong
    arg::SimpleOp
end

show(io::IO, a::Strong) = print(io, "strong(", a.arg, ")")

"""
    strong(op)

declare a conserved quantity as a strong symmetry rather than the weak one `conserve`
assumes by default.

A weak symmetry asks only that the density matrix commute with the charge: the mixed index
holds the difference of the charges of ket and bra, every jump operator of definite charge
preserves it, particle loss and gain included, and a state may mix several sectors.

A strong symmetry asks that every jump operator commute with the charge. Ket and bra are then
conserved separately, which gives finer blocks, but a state lives in a single sector, as a
pure one does, and a jump of non zero charge is refused.

Use it when every dissipator commutes with the quantity, as dephasing does.

# Examples

    Fermion(conserve = strong(N))          # dephasing, `L = N`
    Electron(conserve = (strong(Ntot), 2Sz))
"""
strong(a::SimpleOp) = Strong(a)
strong(a::Strong) = a
strong(a) = error("a conserved quantity is one operator acting on one site, and $a is not")

"""
    struct Conserved

what a site or a system conserves: a list of names, each marked strong or weak, without the
charges, which belong to the site and never change. It prints as the value `conserve` would
be given, so that what `symmetries` reports can be given back to `weaken`.

# Examples

    symmetries(system)                     # (2Sz, strong(Ntot))
    weaken(state, symmetries(system))      # the identity, by construction
"""
struct Conserved
    names::Vector{Tuple{String, Bool}}
    # sorted by name, as ITensors sorts the components of a QN: the order in which a
    # declaration or a target names its quantities says nothing about what is conserved
    Conserved(names) = new(sort(collect(Tuple{String, Bool}, names); by = first))
end

# defined together, as `Op` does: a `Set` or a `Dict` picks its bucket by `hash` and only
# then compares, so two equal values that hash apart would sit in different buckets. The
# field being a vector, the fallback would compare identities and call two equal lists
# different
==(a::Conserved, b::Conserved) = a.names == b.names
hash(a::Conserved, h::UInt) = hash(a.names, hash(Conserved, h))

show(io::IO, c::Conserved) =
    if isempty(c.names)
        print(io, "()")
    else
        one(n, st) = st ? "strong($n)" : n
        print(io, length(c.names) == 1 ? one(c.names[1]...) :
                  "(" * join([ one(n, st) for (n, st) in c.names ], ", ") * ")")
    end

"""
    spec_names(spec)

the quantities a target names, as the `(name, strong)` pairs `Conserved` holds.

It reads both the operators `conserve` takes and what `symmetries` returns. The names are all
`weaken` needs, since it recomputes no charge; a declaration does, which is why `conserve`
asks for the operators themselves.
"""
spec_names(c::Conserved) = c.names
spec_names(::Tuple{}) = Tuple{String, Bool}[]
spec_names(spec::Tuple) = reduce(vcat, map(spec_names, spec))
spec_names(a::Strong) = [ (obs_name(a.arg), true) ]
spec_names(a::SimpleOp) = [ (obs_name(a), false) ]
spec_names(a) = error("$a does not name a conserved quantity")

"""
    symmetries(::AbstractSite)
    symmetries(::System)

what a site or a system conserves, and whether strongly or weakly, as a `Conserved`. It prints
as the value `conserve` is given to declare it, and can be given back to `weaken` as a target.

A system conserves what its sites conserve: a site that does not conserve a quantity does not
keep the others from conserving it.

# Examples

    symmetries(System(4, Electron(conserve = (strong(Ntot), 2Sz))))   # (2Sz, strong(Ntot))
    symmetries(System(4, Fermion()))                                  # ()
    symmetries(System([Qubit(), Fermion(conserve = N)]))              # N
"""
symmetries(site::AbstractSite) =
    Conserved([ (q[1], q[4]) for q in decode_conserve(conserved(site)) ])

"""
    one_step_down(c::Conserved)

the target of `weaken` when none is given: every strong quantity made weak, or, when none is
strong, every quantity dropped. Repeating it walks strong, then weak, then nothing, and stops
there.
"""
one_step_down(c::Conserved) =
    any(last, c.names) ? Conserved([ (n, false) for (n, _) in c.names ]) :
                         Conserved(Tuple{String, Bool}[])

"""
    check_target(source, target, what)

refuse a target that is not a weakening of `source`: a quantity may be dropped or made weak,
not added or made strong, since the finer blocks of a strong symmetry cannot be recovered once
merged.
"""
function check_target(source::Conserved, target::Conserved, what)
    for (name, strong) in target.names
        k = findfirst(q -> q[1] == name, source.names)
        if isnothing(k)
            error("$what does not conserve $name, and weakening cannot start conserving " *
                  "what was not conserved")
        end
        if strong && !source.names[k][2]
            error("$name is conserved weakly, and weakening cannot make it strong: the finer " *
                  "blocks of a strong symmetry are not recoverable from the coarser ones")
        end
    end
    return nothing
end

"""
    retarget(site, target)

the `conserve` string of `site` once it conserves only what `target` names, a name the site
does not conserve being ignored. The charges are those the site already holds: weakening
recomputes none, which is why it needs no operator where a declaration does.
"""
function retarget(site::AbstractSite, target::Conserved)
    qs = decode_conserve(conserved(site))
    parts = String[]
    for (name, strong) in target.names
        k = findfirst(q -> q[1] == name, qs)
        if isnothing(k)
            continue
        end
        (_, modulus, charges, _) = qs[k]
        head = modulus == 1 ? name : "$name%$modulus"
        push!(parts, (strong ? head * "!" : head) * ":" * join(charges, ","))
    end
    return join(parts, ";")
end

"""
    transitions(source, target)

the names to collapse onto their weak form and the names to drop, going from `source` to
`target`. See `weak_qn`.
"""
function transitions(source::Conserved, target::Conserved)
    collapse = String[]
    drop = String[]
    for (name, strong) in source.names
        k = findfirst(q -> q[1] == name, target.names)
        if isnothing(k)
            push!(drop, name)
        elseif strong && !target.names[k][2]
            push!(collapse, name)
        end
    end
    return collapse, drop
end

function weaken(site::AbstractSite, target::Conserved)
    # a site conserving nothing has nothing to weaken, whatever its conserve field holds
    if isempty(conserved(site))
        return site
    end
    t = typeof(site)
    return t(( f === :conserve ? retarget(site, target) : getfield(site, f)
               for f in fieldnames(t) )...)
end

"""
    check_charges(sites)

refuse sites whose conserved quantities cannot live together on one system, which would
otherwise go wrong without a word:
- a name conserved strongly on one site and weakly on another, or with two different moduli,
  would stand for two different charges, which the flux of a state would add;
- the starred name of a strong quantity may already name another quantity, merging two charges;
- ITensors allows four components to a charge, and a strong quantity takes two.
"""
function check_charges(sites::Vector{<:AbstractSite})
    kind = Dict{String, Bool}()
    modulus = Dict{String, Int}()
    for site in sites, (name, m, _, st) in decode_conserve(conserved(site))
        if get(kind, name, st) ≠ st
            error("$name is conserved strongly on one site and weakly on another, so its " *
                  "name would stand for two different charges on the same system")
        end
        if get(modulus, name, m) ≠ m
            error("$name is conserved modulo $(modulus[name]) on one site and modulo $m on " *
                  "another: give the two quantities different names with named")
        end
        kind[name] = st
        modulus[name] = m
    end
    for (name, st) in kind
        if st && haskey(kind, name * "*")
            error("$name is conserved strongly, so its bra goes under $(name)*, which is " *
                  "already the name of another conserved quantity")
        end
    end
    n = sum(st -> st ? 2 : 1, values(kind); init = 0)
    if n > 4
        error("these sites conserve $(length(kind)) quantities, which take $n of the four " *
              "components ITensors allows, a strong one costing two")
    end
    return nothing
end

"""
    conserve_names(s)

the conserved quantities a site records, as `Conserved` prints them, which is how the site
prints them. A string `decode_conserve` cannot read is given back unchanged, since `show` must
print whatever a site holds.
"""
function conserve_names(s::AbstractString)
    try
        return sprint(show, Conserved([ (q[1], q[4]) for q in decode_conserve(s) ]))
    catch e
        # what `decode_conserve` raises on a string it cannot read, and nothing else
        if !(e isa Union{ErrorException, ArgumentError})
            rethrow()
        end
        return String(s)
    end
end

"""
    charged_state(a, inds, what, site)

the tensor of the local state `what`, of array `a`, on the indices `inds` of `site`, refused
by a message naming the state and the site when it has no definite charge. A state of a
charged site must lie in one sector: `"Up"` and `"Dn"` do, `"+"`, their sum, does not. The
question is asked here, see `has_definite_flux`, rather than left to the `Fluxes not all
equal` of ITensors, raised where neither the state nor the site is in sight.
"""
function charged_state(a::AbstractArray, inds, what, site::AbstractSite)
    if !has_definite_flux(a, inds)
        error("the state $(repr(what)) of site $(typeof(site)) spreads over several charges " *
              "of $(conserve_names(conserved(site))), so it has none of its own and cannot " *
              "be used where that is conserved")
    end
    return charged_itensor(a, inds)
end

"""
    show(io, ::AbstractSite)

print a site as the call that builds it, `Qubit()` or `Fermion(conserve = N)`, the conserved
quantities under their names rather than as the charges recorded.

The trailing fields holding `""` or `nothing` are left out, and only those: dropping one in
the middle would misalign the arguments with the fields, `Site(nothing, 3)` printing as
`Site(3)`.
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
