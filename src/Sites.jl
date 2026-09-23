export AbstractSite, dim, Index, string_state, identity_operator, state, flux, weaken,
       symmetries
export @def_operators, @def_states, @create_site_module, conserve_string, strong

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
    site_index(site, charged)

the ITensor index of a site for pure representations.

`charged` says whether the system the site belongs to carries quantum numbers, which a site
conserving nothing cannot know on its own: in such a system it takes a trivial index, a
single sector of charge zero holding the whole space. An MPS cannot mix indices that carry
charges with indices that do not, and every operator remains available on a trivial index,
every matrix element sitting in the one block.

The sectors come one per basis state rather than merged by charge, because merging would
reorder the basis whenever equal charges are not contiguous, as `parity(N)` on a boson gives
0, 1, 0, 1.
"""
function site_index(site::AbstractSite, charged::Bool)
    n = dim(site)
    tg = "$(nameof(typeof(site))), Site"
    qs = decode_conserve(conserved(site))
    if isempty(qs)
        return charged ? Index(QN() => n; tags = tg) : Index(n; tags = tg)
    end
    for (name, _, charges, _) in qs
        if length(charges) ≠ n
            error("site $(typeof(site)) records $(length(charges)) charges for $name but " *
                  "has $n basis states")
        end
    end
    return Index([ QN([(name, charges[k], modulus)
                       for (name, modulus, charges, _) in qs]...) => 1
                   for k in 1:n ]...; tags = tg)
end

"""
    weak_qn(q, names)
    weak_index(i, names)

the charge, or the index, with each strong quantity in `names` collapsed onto its weak form.

Keeping the ket and the bra apart records `X` and `X*`; asking for the same quantity weakly
records their difference, which is what the sum of the two components is, the bra having been
daggered. This map is a homomorphism of the charge group, so it carries a flux to a flux and
a relation between blocks to the same relation: an index may be relabelled with it and every
tensor built on it stays consistent, with no data moved.

The blocks are not merged. Several may end up under one charge, which an index allows, and
that is what lets the relabelling cost nothing.
"""
function weak_qn(q::QN, collapse, drop)
    vals = Tuple{String, Int, Int}[]
    for v in q.data
        nm = String(ITensors.name(v))
        if isempty(nm)
            continue
        end
        base = endswith(nm, "*") ? nm[1:end-1] : nm
        if base in drop
            continue
        end
        if !(base in collapse)
            push!(vals, (nm, ITensors.val(v), ITensors.modulus(v)))
            continue
        end
        k = findfirst(x -> x[1] == base, vals)
        if isnothing(k)
            push!(vals, (base, ITensors.val(v), ITensors.modulus(v)))
        else
            vals[k] = (base, vals[k][2] + ITensors.val(v), vals[k][3])
        end
    end
    return isempty(vals) ? QN() : QN(vals...)
end

weak_index(i::Index, collapse, drop) =
    if (isempty(collapse) && isempty(drop)) || !hasqns(i)
        i
    else
        Index([ weak_qn(q, collapse, drop) => d for (q, d) in space(i) ]...;
              tags = tags(i), plev = plev(i), dir = dir(i))
    end

"""
    bra_index(i, site)

the index the bra of `i` is carried by, which is `i` itself unless the site declares a
strong symmetry. See `strong` and `star`.
"""
bra_index(i::Index, site::AbstractSite) = star(i, strong_names(site))

"""
    Index(::AbstractSite)

return an ITensor.Index for the given site for pure representations.

It carries the quantum numbers the site declares, and none when it declares none. A site
conserving nothing inside a system where another one does takes a trivial index instead, see
`site_index`, which is a property of the system rather than of the site.
"""
Index(site::AbstractSite) = site_index(site, !isempty(conserved(site)))

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
    mixed_index(i, site)

the index of the mixed representation pairing the ket `i` with the bra of the same site.

The bra is daggered, so that the charge of ``|m\\rangle\\langle n|`` is the difference of
those of ``m`` and ``n`` rather than their sum. This is the pairing `mix(::State)` produces,
its tensors being contracted as `t * dag(t')`, and on an index without charges the dag is a
no operation. A site conserving something strongly stars its bra instead, which keeps the two
charges apart; see `strong`.

This is internal: the index it draws is a fresh one, of the right space but of an identity of
its own, so it contracts with nothing. What a caller wants is the index the system drew,
`SysIndex{Mixed}(system, i)`.
"""
mixed_index(i::Index, site::AbstractSite) =
    addtags(combinedind(combiner(i, dag(bra_index(i, site)'); tags = tags(i))), "Mixed")

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
function site_charges(op::SimpleOp, site::AbstractSite; tol::Float64 = charge_tol)
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

site_charges(op::GenericOp{Pure}, ::AbstractSite) =
    error("a conserved quantity acts on one site, and $op acts on several")

# the modulus of a ModOp is carried rather than read back, ±1 being unreadable, and the
# charges are those of its argument taken modulo it
function site_charges(a::ModOp{1}, site::AbstractSite; tol::Float64 = charge_tol)
    m, q = site_charges(a.arg, site; tol)
    if m ≠ 1
        error("cannot take $(a.arg) modulo $(a.modulus) on site $(typeof(site)): it " *
              "already carries a charge modulo $m")
    end
    return (a.modulus, mod.(q, a.modulus))
end

"""
    show_charges(d)

a list of charge differences as an error message shows it, `2Sz=2` or `Ntot=-1,2Sz=-1`
"""
show_charges(d) = join(["$name=$val" for (name, val, _) in d], ",")

"""
    flux(op, site)

the charge an operator carries on a site, as a `QN`: the difference between the charges of
the states it connects. It is the zero charge for an operator commuting with everything the
site conserves, and `QN()` for a site conserving nothing.

An operator connecting states whose charges differ in more than one way has no flux at all
and cannot be used where that quantity is conserved. `X` raises and lowers `N` at once and is
refused, while under `parity(N)` it carries `1`, the two differences becoming the same one
modulo 2.

The difference is taken modulo the charge, without which `Xd` would be refused although it
generates the very symmetry `Zd` records: its wrap around connects the last state to the
first, a difference of `1 - d` rather than of 1.

# Examples

    flux(N, Fermion(conserve = N))        # QN("N",0)
    flux(Sp, Qubit(conserve = 2Sz))       # QN("2Sz",2), in units of the declared charge
    flux(Xd, Qudit(3, conserve = Zd))     # QN("Zd",1,3)
"""
flux(op::SimpleOp, site::AbstractSite; tol::Float64 = charge_tol) =
    charge_flux(matrix(op, site), op, site; tol)

"""
    charge_flux(m, what, site)

the flux of a matrix already computed, `what` being what to name if it has none. This is what
`flux` answers and what the tensor of an operator is checked with, so that a matrix which does
not fit the charges of its site is refused by a message naming the operator rather than by the
`Fluxes not all equal` of ITensors, raised from somewhere neither the operator nor the site is
in sight.
"""
function charge_flux(m::Matrix, what, site::AbstractSite; tol::Float64 = charge_tol)
    qs = decode_conserve(conserved(site))
    if isempty(qs)
        return QN()
    end
    found = nothing
    for i in axes(m, 1), j in axes(m, 2)
        if abs(m[i, j]) ≤ tol
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

flux(op::GenericOp{Pure}, site::AbstractSite) =
    error("flux is only defined for one site operators, and $op acts on several")

struct Strong
    arg::SimpleOp
end

show(io::IO, a::Strong) = print(io, "strong(", a.arg, ")")

"""
    strong(op)

declare a conserved quantity as a strong symmetry rather than the weak one `conserve`
assumes by default.

A weak symmetry only asks that the density matrix commute with the charge, which is what
the mixed index records when it holds the difference of the two charges of
``|m\\rangle\\langle n|``. Every jump operator of definite charge preserves it, particle
loss and gain included, and a state may mix several sectors.

A strong symmetry asks more: that every jump operator commute with the charge. The ket and
the bra are then conserved separately, the mixed index keeps them apart instead of holding
their difference, and the blocks are finer. In exchange a state lives in a single sector, as
a pure one does, and a jump of non zero charge is refused: it is not a strong symmetry.

Use it when every dissipator commutes with the quantity, as dephasing does, and leave it out
otherwise.

# Examples

    Fermion(conserve = strong(N))          # dephasing, `L = N`
    Electron(conserve = (strong(Ntot), 2Sz))
"""
strong(a::SimpleOp) = Strong(a)
strong(a::Strong) = a
strong(a) = error("a conserved quantity is one operator acting on one site, and $a is not")

"""
    struct Conserved

what a site or a system conserves, as a list of names each marked strong or weak.

It prints as the expression that would declare it, so that what a system reports can be read
back and given to `weaken`. The charges themselves are left out: they belong to the site and
never change, only the way the ket is paired with the bra does.

# Examples

    symmetries(system)                     # (strong(Ntot), 2Sz)
    weaken(state, symmetries(system))      # the identity, by construction
"""
struct Conserved
    names::Vector{Tuple{String, Bool}}
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

the quantities a target names, as `Conserved` holds them.

Both vocabularies are accepted: the operators one writes by hand, as `conserve` takes them,
and what `symmetries` reports. `weaken` needs no more than the names, since it recomputes no
charge; only a declaration does, which is why `conserve` asks for the operators themselves.
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

what is conserved, and how. It prints as the value `conserve` would be given to declare it,
and it can be given back to `weaken` as a target.

# Examples

    symmetries(System(4, Electron(conserve = (strong(Ntot), 2Sz))))   # (strong(Ntot), 2Sz)
    symmetries(System(4, Fermion()))                                  # ()
"""
symmetries(site::AbstractSite) =
    Conserved([ (q[1], q[4]) for q in decode_conserve(conserved(site)) ])

"""
    one_step_down(c::Conserved)

the target `weaken` aims at when none is given: every strong quantity asked for weakly, or,
when none is strong, every quantity dropped.

The level of a system is the strongest of its quantities, and this takes it down one notch.
Repeating it walks strong, then weak, then nothing, and stops there.
"""
one_step_down(c::Conserved) =
    any(last, c.names) ? Conserved([ (n, false) for (n, _) in c.names ]) :
                         Conserved(Tuple{String, Bool}[])

"""
    check_target(source, target, what)

refuse a target that is not a weakening of `source`.

A quantity may be dropped or asked for less strongly; it may not be invented, nor made
stronger, the finer blocks of a strong symmetry not being recoverable from the coarser ones
once they have been merged.
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

what a site records once it conserves what `target` names, and only that.

The charges are the ones the site already holds: weakening never recomputes them, which is
why it needs no operator where a declaration does.
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

weaken(site::AbstractSite, target::Conserved) =
    let t = typeof(site)
        t(( f === :conserve ? retarget(site, target) : getfield(site, f)
            for f in fieldnames(t) )...)
    end

"""
    check_charges(sites)

refuse a list of sites whose conserved quantities cannot live together on one system.

Three things would otherwise go wrong without a word. A name conserved strongly on one site
and weakly on another would stand for the charge of the ket on the first and for a difference
on the second, and the flux of a state would add the two. The star a strong quantity gives
its bra may be the name of another quantity, which would merge two charges into one. And the
links of a state carry every component of every site, which ITensors limits to four, a strong
quantity costing two of them.
"""
function check_charges(sites::Vector{<:AbstractSite})
    kind = Dict{String, Bool}()
    for site in sites, (name, _, _, st) in decode_conserve(conserved(site))
        if get(kind, name, st) ≠ st
            error("$name is conserved strongly on one site and weakly on another, so its " *
                  "name would stand for two different charges on the same system")
        end
        kind[name] = st
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
    strong_names(site)

the names of the quantities the site conserves strongly, empty when it conserves none that
way. See `strong`.
"""
strong_names(site::AbstractSite) =
    [ q[1] for q in decode_conserve(conserved(site)) if q[4] ]

"""
    star(q::QN, names)
    star(i::Index, names)

the charge, or the index, with every component named in `names` renamed to carry a star.

This is what separates the bra from the ket: a strong symmetry conserves the two sides
apart, so the bra holds its charges under other names and combining the pair keeps them
rather than subtracting them. It applies to a whole index and not only to a site one,
because the links of a state carry the same charges and must be renamed with it, or the
two halves of the same tensor would count in two different ways. Renaming nothing gives the
index back as it is, which is every case without a strong symmetry.
"""
function star(q::QN, names)
    vs = Tuple{String, Int, Int}[]
    for v in q.data
        n = String(ITensors.name(v))
        if !isempty(n)
            push!(vs, (n in names ? n * "*" : n, ITensors.val(v), ITensors.modulus(v)))
        end
    end
    return isempty(vs) ? QN() : QN(vs...)
end

function star(i::Index, names)
    if isempty(names) || !hasqns(i)
        return i
    end
    return Index([ star(q, names) => d for (q, d) in space(i) ]...;
                 tags = tags(i), plev = plev(i), dir = dir(i))
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
    parts = map(ops) do spec
        op = spec isa Strong ? spec.arg : spec
        modulus, q = site_charges(op, site)
        name = obs_name(op)
        if endswith(name, '!')
            error("cannot conserve $name: a name ending in ! cannot be told from the mark " *
                  "a site puts on a strong symmetry")
        end
        head = modulus == 1 ? name : "$name%$modulus"
        # a strong symmetry is marked on the quantity and not on the site, so that one site
        # may hold both kinds, and at the end of the head so that the name and the modulus
        # are read exactly as before
        return (spec isa Strong ? head * "!" : head) * ":" * join(q, ",")
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
    conserved(site)

what a site conserves, in the form `conserve_string` produces, and the empty string when it
conserves nothing.

The `conserve` field is optional, a site type that can have no conserved quantity simply not
declaring it, so this is what everything else reads rather than the field itself. A field of
that name holding something other than a string is not one, and the site conserves nothing.
"""
conserved(site::AbstractSite) =
    if hasfield(typeof(site), :conserve) && getfield(site, :conserve) isa AbstractString
        site.conserve
    else
        ""
    end

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
        st = endswith(head, '!')
        if st
            head = head[1:end-1]
        end
        j = findlast(==('%'), head)
        name = isnothing(j) ? String(head) : String(head[1:j-1])
        push!(ns, st ? "strong($name)" : name)
    end
    return length(ns) == 1 ? ns[1] : "(" * join(ns, ", ") * ")"
end

"""
    charged_state(f, i::Index, what, site)

the tensor `f` builds for the local state `what` on the site index `i`, refused by a message
naming the state and its site when the charges of the site cannot carry it.

A state of a charged site belongs to one sector: `"Up"` and `"Dn"` do, `"+"` does not, being
their sum, and no amount of bookkeeping gives a superposition of two charges a charge of its
own. ITensors says `Fluxes not all equal` from a place where neither the state nor the site
is in sight, so it is said here instead. The tensor comes from a function rather than being
passed in because a state vector, a density matrix and the target of a `SetState` are laid
on their indices in three different ways, and one message covers them all.
"""
function charged_state(f, i::Index, what, site::AbstractSite)
    if !hasqns(i)
        return f()
    end
    try
        return f()
    catch e
        if !(e isa ErrorException)
            rethrow()
        end
        error("the state $(repr(what)) of site $(typeof(site)) spreads over several charges " *
              "of $(conserve_names(conserved(site))), so it has none of its own and cannot " *
              "be used where that is conserved")
    end
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
