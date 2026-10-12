# What a site conserves: the charges read off the conserved operators, kept typed in the site,
# strong or weak, the relabellings of the charges they need, and the check of a site at its
# first use. The later layers lay matrices on charged indices with `conserved`, `strong_names`,
# `star`, `basis_charges`, `has_definite_flux`, `charged_itensor`, `common_charge`,
# `charge_flux`, `FluxError`, `weak_qn` and `adjoint_index`.

export Conservation, conserving, strong, symmetries, weaken, flux

"""
    charge_tol

how far the eigenvalues of a conserved quantity may lie from the charges they stand for,
relative to the largest of them when it is above one
"""
const charge_tol = 1e-14

"""
    short(x)

a number rounded to two significant digits, as an error message shows a deviation
"""
short(x::Real) = round(x; sigdigits = 2)

"""
    site_charges(op, site)

the modulus and the charge of each basis state of `site` for the conserved quantity `op`, as
`(modulus, charges)`. `op` is a charge, diagonal on the site with integer eigenvalues, or a charge
modulo `m`, written `mod(Q, m)` or `parity(Q)`: the expression carries the modulus, which the
spectrum is never asked for.
"""
function site_charges(op::SimpleOp, site::AbstractSite; tol::Float64 = charge_tol)
    # a renamed operator, or one declared for the site, reads its charges off its definition,
    # which may carry a modulus, as parity(N) does
    if op isa Operator
        def = isnothing(op.expr) ? operator_definition(site, Val(Symbol(op.name))) : op.expr
        if def isa Function
            def = def(site)
        end
        if def isa GenericOp
            return site_charges(def, site; tol)
        end
    end
    m = matrix(op, site)
    d = diag(m)
    # relative to the largest charge, whose rounding grows with it
    tol *= max(1, maximum(abs, d))
    off = norm(m - Diagonal(d))
    if off > tol
        error("$op is not diagonal on site $(typeof(site)), off by $(short(off)), so it " *
              "cannot be a conserved quantity: a charge is carried by each basis state")
    end
    to_int = maximum(max(abs(imag(x)), abs(real(x) - round(real(x)))) for x in d)
    if to_int > tol
        error("the eigenvalues of $op on site $(typeof(site)) miss the integers by " *
              "$(short(to_int)), so they are not charges. Half integer ones are written " *
              "doubled, 2Sz rather than Sz, and a charge modulo m is written mod(Q, m)")
    end
    return (1, Int.(round.(real.(d))))
end

site_charges(op::GenericOp{IsOperator}, ::AbstractSite) =
    error("a conserved quantity acts on one site, and $op acts on several")

# a power or the adjoint of a charge modulo m is one modulo m, whatever its eigenvalues: those
# of Z^2 on Qudit(4), all ±1, would be read as integers
function site_charges(a::Union{IntPowOp{IsOperator, 1}, GenPowOp{IsOperator, 1}, DagOp{1}}, site::AbstractSite;
                      tol::Float64 = charge_tol)
    p = a isa DagOp ? -1 : a.expo
    if isreal(p) && isinteger(real(p))
        mq = try
            site_charges(a.arg, site; tol)
        catch e
            if !(e isa ErrorException)
                rethrow()
            end
            nothing
        end
        if !isnothing(mq) && first(mq) ≠ 1
            m, q = mq
            return (m, mod.(Int(real(p)) .* q, m))
        end
    end
    return invoke(site_charges, Tuple{SimpleOp, AbstractSite}, a, site; tol)
end

# the modulus of a ModOp is carried rather than read back, ±1 being ambiguous
function site_charges(a::ModOp{1}, site::AbstractSite; tol::Float64 = charge_tol)
    m, q = site_charges(a.arg, site; tol)
    if m ≠ 1
        error("cannot take $(a.arg) modulo $(a.modulus) on site $(typeof(site)): it " *
              "already carries a charge modulo $m")
    end
    return (a.modulus, mod.(q, a.modulus))
end

"""
    Quantity(name, modulus, charges, strong)

a conserved quantity of a site: its name, its modulus, 1 for a charge over the integers, the
charge of each basis state, and whether it is conserved strongly
"""
struct Quantity
    name::String
    modulus::Int
    charges::Vector{Int}
    strong::Bool
end

(a::Quantity == b::Quantity) =
    (a.name, a.modulus, a.charges, a.strong) == (b.name, b.modulus, b.charges, b.strong)
hash(a::Quantity, h::UInt) = hash((a.name, a.modulus, a.charges, a.strong), hash(Quantity, h))

"""
    Conservation(quantities = Quantity[])
    Conservation(site, spec)

what a site conserves, its quantities sorted by name, as ITensors sorts the components of a
charge: built from the operators `spec`, one or a tuple of them, each possibly marked `strong`,
whose charges are read on `site`. A site type fills its field with `conserving`.
"""
struct Conservation
    quantities::Vector{Quantity}
    Conservation(qs::Vector{Quantity} = Quantity[]) = new(sort(qs; by = q -> q.name))
end

(a::Conservation == b::Conservation) = a.quantities == b.quantities
hash(a::Conservation, h::UInt) = hash(a.quantities, hash(Conservation, h))

"""
    Strong

a conserved quantity declared strong, as `strong` returns it
"""
struct Strong
    arg::SimpleOp
end

"""
    strong(op)

declare a conserved quantity as a strong symmetry rather than the weak one `conserve` assumes.
Under a strong symmetry every jump operator commutes with the charge: ket and bra are
conserved separately, which gives finer blocks, and a state lives in a single sector.

# Examples

    Fermion(conserve = strong(N))
    Electron(conserve = (strong(Ntot), 2Sz))
"""
strong(a::SimpleOp) = Strong(a)
strong(a::Strong) = a
strong(a) = error("a conserved quantity is one operator acting on one site, and $a is not")

function Conservation(site::AbstractSite, spec)
    ops = spec isa Tuple ? collect(spec) : [spec]
    names = [ obs_name(o isa Strong ? o.arg : o) for o in ops ]
    for n in names
        if count(==(n), names) > 1
            error("$n is given twice in conserve")
        end
    end
    qs = map(ops) do o
        op = o isa Strong ? o.arg : o
        if !(op isa SimpleOp)
            error("a conserved quantity is one operator acting on one site, and $op is not")
        end
        modulus, charges = site_charges(op, site)
        name = obs_name(op)
        if endswith(name, '*')
            error("cannot conserve $name: a name ending in * is that of the bra of a strong quantity")
        end
        limit = ITensors.SmallStrings.smallLength - (o isa Strong ? 1 : 0)
        if length(name) > limit
            error("cannot conserve $name: ITensors takes names of at most $limit characters " *
                  "here, give it a shorter one with named")
        end
        Quantity(name, modulus, charges, o isa Strong)
    end
    return Conservation(Vector{Quantity}(qs))
end

"""
    conserving(T, params...; conserve)

the site of type `T` and parameters `params` conserving `conserve`, one operator or a tuple of
them, each possibly marked `strong`: the charges are read on the site built conserving
nothing. A site type writes its keyword constructor with it.

# Examples

    Qudit(dim::Int; conserve = ()) = conserving(Qudit, dim; conserve)
"""
function conserving(t::Type{<:AbstractSite}, params...; conserve)
    bare = t(params..., Conservation())
    return t(params..., Conservation(bare, conserve))
end

"""
    conserved(site)

what a site conserves, empty for a site type without a `Conservation` field
"""
conserved(site::AbstractSite) =
    if hasfield(typeof(site), :conserve) && getfield(site, :conserve) isa Conservation
        getfield(site, :conserve)
    else
        Conservation()
    end

"""
    strong_names(site)

the names of the quantities the site conserves strongly
"""
strong_names(site::AbstractSite) = [ q.name for q in conserved(site).quantities if q.strong ]

"""
    basis_charges(site)

the charge of each basis state of `site`, a `QN` holding every quantity it conserves
"""
function basis_charges(site::AbstractSite)
    n = dim(site)
    qs = conserved(site).quantities
    for q in qs
        if length(q.charges) ≠ n
            error("site $(typeof(site)) records $(length(q.charges)) charges for $(q.name) " *
                  "but has $n basis states")
        end
    end
    return [ QN([ (q.name, q.charges[k], q.modulus) for q in qs ]...) for k in 1:n ]
end

"""
    Conserved

what a site or a system conserves, a list of names each marked strong or weak, without the
charges. It prints as `conserve` would be given, so that what `symmetries` reports can be given
back to `weaken`.
"""
struct Conserved
    names::Vector{Tuple{String, Bool}}
    Conserved(names) = new(sort(collect(Tuple{String, Bool}, names); by = first))
end

(a::Conserved == b::Conserved) = a.names == b.names
hash(a::Conserved, h::UInt) = hash(a.names, hash(Conserved, h))

"""
    spec_names(spec)

the quantities a target names, as the `(name, strong)` pairs `Conserved` holds
"""
spec_names(c::Conserved) = c.names
spec_names(::Tuple{}) = Tuple{String, Bool}[]
function spec_names(spec::Tuple)
    names = reduce(vcat, map(spec_names, spec))
    for (n, _) in names
        if count(x -> first(x) == n, names) > 1
            error("$n is given twice")
        end
    end
    return names
end
spec_names(a::Strong) = [ (obs_name(a.arg), true) ]
spec_names(a::SimpleOp) = [ (obs_name(a), false) ]
spec_names(a) = error("$a does not name a conserved quantity")

"""
    symmetries(::AbstractSite)
    symmetries(::System)

what a site or a system conserves, and whether strongly or weakly, as a `Conserved`, which
prints as `conserve` is given and can be given back to `weaken`.

# Examples

    symmetries(System(4, Electron(conserve = (strong(Ntot), 2Sz))))   # (2Sz, strong(Ntot))
"""
symmetries(site::AbstractSite) = Conserved([ (q.name, q.strong) for q in conserved(site).quantities ])

"""
    one_step_down(c::Conserved)

the target of `weaken` when none is given: every strong quantity made weak or, when none is
strong, every quantity dropped
"""
one_step_down(c::Conserved) =
    any(last, c.names) ? Conserved([ (n, false) for (n, _) in c.names ]) :
                         Conserved(Tuple{String, Bool}[])

"""
    check_target(source, target, what)

refuse a target that is not a weakening of `source`: a quantity may be dropped or made weak,
not added or made strong
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
    transitions(source, target)

the names to collapse onto their weak form and the names to drop, going from `source` to
`target`
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

"""
    retarget(site, target)

what `site` conserves once restricted to what `target` names, as strongly as it asks, the
charges kept
"""
retarget(site::AbstractSite, target::Conserved) =
    Conservation(Quantity[ Quantity(q.name, q.modulus, q.charges, strong)
                           for (name, strong) in target.names
                           for q in conserved(site).quantities if q.name == name ])

function weaken(site::AbstractSite, target::Conserved)
    if isempty(conserved(site).quantities)
        return site
    end
    t = typeof(site)
    return rebuild_site(t, [ f === :conserve ? retarget(site, target) : getfield(site, f)
                             for f in fieldnames(t) ], "its conserved quantities to be weakened")
end

"""
    check_charges(sites)

refuse sites whose conserved quantities cannot live together on one system: one name strong
on a site and weak on another, or of two moduli; a starred name of a strong quantity naming
another; more than the four components ITensors allows, a strong quantity taking two
"""
function check_charges(sites::Vector{<:AbstractSite})
    kind = Dict{String, Bool}()
    modulus = Dict{String, Int}()
    for site in sites, q in conserved(site).quantities
        if get(kind, q.name, q.strong) ≠ q.strong
            error("$(q.name) is conserved strongly on one site and weakly on another, so its " *
                  "name would stand for two different charges on the same system")
        end
        if get(modulus, q.name, q.modulus) ≠ q.modulus
            error("$(q.name) is conserved modulo $(modulus[q.name]) on one site and modulo " *
                  "$(q.modulus) on another: give the two quantities different names with named")
        end
        kind[q.name] = q.strong
        modulus[q.name] = q.modulus
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

############### Charges of indices and tensors ###############

"""
    qn_components(q::QN)

the components of the charge `q`, as `(name, value, modulus)`, the unused slots left out
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

the index `i` with the charge `q` of each block replaced by `f(q)`, the dimensions, tags, prime
level and direction kept
"""
map_charges(f, i::Index) =
    Index([ f(q) => d for (q, d) in space(i) ]...; tags = tags(i), plev = plev(i), dir = dir(i))

"""
    weak_qn(q, collapse, drop)
    weak_index(i, collapse, drop)

the charge, or the index, with the quantities of `drop` left out and the strong quantities of
`collapse` made weak, their components `X` and `X*` summed into `X`. Both are homomorphisms of
the charge group, so every tensor on a relabelled index stays consistent with no data moved.
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

the charge, or the index, with each component `X` named in `names` renamed `X*`: under a
strong symmetry the bra carries its charges under the starred names
"""
star(q::QN, names) = QN([ (n in names ? n * "*" : n, v, m) for (n, v, m) in qn_components(q) ]...)

function star(i::Index, names)
    if isempty(names) || !hasqns(i)
        return i
    end
    return map_charges(q -> star(q, names), i)
end

"""
    adjoint_qn(q::QN, names)
    adjoint_index(i::Index, names)

the charge, or the index, relabelled for the adjoint, `names` being the quantities conserved
strongly: `(X, X*)` becomes `(-X*, -X)`, a weak quantity is negated
"""
adjoint_qn(q::QN, names) =
    QN([ (endswith(n, "*") ? n[1:end-1] : n in names ? n * "*" : n, -v, m)
         for (n, v, m) in qn_components(q) ]...)

adjoint_index(i::Index, names) = hasqns(i) ? map_charges(q -> adjoint_qn(q, names), i) : i

"""
    bra_index(i, site)

the index carrying the bra of `i`: `i` with the quantities `site` conserves strongly starred
"""
bra_index(i::Index, site::AbstractSite) = star(i, strong_names(site))

"""
    index_charges(i)

the charge each value of the index `i` brings to the flux of a tensor, negated on an incoming
index
"""
index_charges(i::Index) =
    reduce(vcat, [ fill(dir(i) == ITensors.In ? -q : q, n) for (q, n) in space(i) ]; init = QN[])

"""
    common_charge(a, charges, tol)

the charge every element of the array `a` above `tol` carries, `charges[k][v]` that of value
`v` of axis `k`: `QN()` for no element, `nothing` when two carry different charges
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

whether the array `a` on the indices `inds` carries a definite charge, up to rounding
"""
has_definite_flux(a::AbstractArray, inds) =
    !any(hasqns, inds) || !isnothing(common_charge(a, map(index_charges, inds), rounding_tol * norm(a)))

"""
    charged_itensor(a, inds)

the ITensor of the array `a` on the indices `inds`, dropping what rounding leaves outside the
blocks of charged indices
"""
charged_itensor(a::AbstractArray, inds) =
    if any(hasqns, inds)
        ITensor(a, inds...; tol = rounding_tol * norm(a))
    else
        ITensor(a, inds...)
    end

"""
    show_charges(d)

a list of `(name, value, modulus)` charges as an error message shows it, `Ntot=-1,2Sz=-1`
"""
show_charges(d) = join(["$name=$val" for (name, val, _) in d], ",")

"""
    FluxError(msg)

the refusal of an operator that changes the conserved charges by no definite amount, which
`expect` answers by measuring it on the state conserving nothing
"""
struct FluxError <: Exception
    msg::String
end

showerror(io::IO, e::FluxError) = print(io, e.msg)

"""
    charge_flux(m, what, site)

the flux of the matrix `m` of `what` on `site`, refused when it has none, see `FluxError`
"""
function charge_flux(m::Matrix, what, site::AbstractSite; tol::Float64 = rounding_tol)
    qs = conserved(site).quantities
    if isempty(qs)
        return QN()
    end
    small = tol * norm(m)
    found = nothing
    for i in axes(m, 1), j in axes(m, 2)
        if abs(m[i, j]) ≤ small
            continue
        end
        d = [ (q.name, q.modulus == 1 ? q.charges[i] - q.charges[j] : mod(q.charges[i] - q.charges[j], q.modulus), q.modulus)
              for q in qs ]
        if isnothing(found)
            found = d
        elseif d ≠ found
            throw(FluxError("$what on site $(typeof(site)) connects charges differing by " *
                            "$(show_charges(found)) and by $(show_charges(d)), so it has no " *
                            "definite flux"))
        end
    end
    return isnothing(found) ? QN() : QN(found...)
end

"""
    flux(op, site)

the charge the one site operator `op` carries on `site`, as a `QN`: the difference between the
charges of the states it connects. An operator connecting states whose charges differ in more
than one way is refused.

# Examples

    flux(N, Fermion(conserve = N))        # QN("N",0)
    flux(X, Qudit(3, conserve = Z))       # QN("Z",1,3)
"""
flux(op::SimpleOp, site::AbstractSite; tol::Float64 = rounding_tol) =
    charge_flux(matrix(op, site), op, site; tol)

flux(op::GenericOp{IsOperator}, ::AbstractSite) =
    error("flux is only defined for one site operators, and $op acts on several")

############### Sites as values ###############

# a site is compared by its fields, a Conservation among them
(a::AbstractSite == b::AbstractSite) =
    typeof(a) == typeof(b) && all(f -> getfield(a, f) == getfield(b, f), fieldnames(typeof(a)))

hash(a::AbstractSite, h::UInt) =
    foldl((h, f) -> hash(no_signed_zero(getfield(a, f)), h), fieldnames(typeof(a)); init = hash(typeof(a), h))

show(io::IO, c::Conserved) =
    if isempty(c.names)
        print(io, "()")
    else
        one(n, st) = st ? "strong($n)" : n
        print(io, length(c.names) == 1 ? one(c.names[1]...) :
                  "(" * join([ one(n, st) for (n, st) in c.names ], ", ") * ")")
    end

show(io::IO, a::Strong) = print(io, "strong(", a.arg, ")")

show(io::IO, c::Conservation) = show(io, Conserved([ (q.name, q.strong) for q in c.quantities ]))

"""
    show(io, ::AbstractSite)

print a site as the call that builds it, `Qubit()` or `Fermion(conserve = N)`, the trailing
fields holding nothing left out
"""
function show(io::IO, site::AbstractSite)
    t = typeof(site)
    args = String[]
    empty = Bool[]
    for f in fieldnames(t)
        v = getfield(site, f)
        if v isa Conservation
            push!(args, "conserve = " * sprint(show, v))
            push!(empty, isempty(v.quantities))
        else
            push!(args, repr(v))
            push!(empty, isnothing(v) || (v isa AbstractString && isempty(v)))
        end
    end
    n = length(args)
    while n > 0 && empty[n]
        n -= 1
    end
    print(io, nameof(t), "(", join(args[1:n], ", "), ")")
end

############### The check of a site at its first use ###############

"""
    declared_names(f, site)

the names `@def_states` or `@def_operators`, as `f`, declared for the type of `site`, read off
the methods they added, those declared on an abstract type, made of the others, left out
"""
function declared_names(f, site::AbstractSite)
    names = String[]
    for m in methods(f)
        sig = Base.unwrap_unionall(m.sig)
        if length(sig.parameters) == 3
            t, v = sig.parameters[2], sig.parameters[3]
            if v isa DataType && v <: Val && !isempty(v.parameters) && typeof(site) <: t &&
                    !isabstracttype(t)
                push!(names, string(v.parameters[1]))
            end
        end
    end
    return names
end

"""
    check_declarations(site)

refuse a site whose declarations do not fit it: parities of the wrong number or not 0 or 1,
states or operators of the wrong size, fermion species not satisfying the anticommutation
relations
"""
function check_declarations(site::AbstractSite)
    d = dim(site)
    p = collect(fermion_parities(site))
    if length(p) ≠ d || !all(in((0, 1)), p)
        error("$site has $d basis states, and fermion_parities gives $p")
    end
    for name in declared_names(state_definition, site)
        local_state(site, name)
    end
    for name in declared_names(operator_definition, site)
        matrix(Operator{1}(name), site)
    end
    species = fermion_species(site)
    if !isnothing(species)
        cs = [ matrix(c, site) for c in species ]
        f = matrix(F, site)
        id = Matrix{Float64}(I, d, d)
        for (a, x) in enumerate(cs), (b, y) in enumerate(cs)
            car = norm(x * y' + y' * x - (a == b ? id : zero(id))) + norm(x * y + y * x)
            if car > rounding_tol * d || !nearly(f * x * f, -x)
                error("the fermion species $(species) of $site are not anticommuting " *
                      "annihilation operators odd under F")
            end
        end
    end
    return nothing
end

const checked_sites = Set{AbstractSite}()
const checked_lock = ReentrantLock()

"""
    check_site(site)

`site`, its declarations checked the first time this value of the site is used
"""
function check_site(site::AbstractSite)
    lock(checked_lock) do
        if !(site in checked_sites)
            check_declarations(site)
            push!(checked_sites, site)
        end
    end
    return site
end
