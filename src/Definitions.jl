# Defining operators and conserved quantities: operators of several sites split into one site
# factors, the declaration of operators for a site type with @def_operators, renaming with
# named, and the charges a site records for what it conserves.

export flux, conserve_string, @def_operators


############### Operators of several sites split into one site factors ###############

"""
    superscripts

the digits 0 to 9 as superscripts, which give the site of a factor in its name, see
`split_matrix`
"""
const superscripts = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")

"""
    subscripts

the digits 0 to 9 as subscripts, which number the factors of a site in their names, see
`split_matrix`
"""
const subscripts = collect("₀₁₂₃₄₅₆₇₈₉")

"""
    script(digits, n)

the number `n` written with `digits`, `superscripts` or `subscripts`
"""
script(digits, n::Int) = join(digits[c - '0' + 1] for c in string(n))

"""
    check_even(name, m, sites, tol)

refuse a matrix that does not commute with `F` on each of its sites. Its one site factors are
placed with no Jordan-Wigner string between them, and the strings of other factors cross them
as if they were even: right for a density, a spin or a pair, all even on each site, wrong for
an operator moving a fermion from one site to another.
"""
function check_even(name, m, sites, tol)
    d = [ dim(s) for s in sites ]
    for (j, s) in enumerate(sites)
        f = matrix(F, s)
        if f == I
            continue
        end
        fj = kron(identity_operator(prod(d[1:j-1])), f, identity_operator(prod(d[j+1:end])))
        if norm(fj * m * fj - m) > tol
            error("$name does not commute with F on its site $j, $s: it moves a fermion " *
                  "there, which its one site factors, placed with no Jordan-Wigner string, " *
                  "cannot do. Give it without its sites, or develop it")
        end
    end
    return nothing
end

"""
    pair_charges(site)

the charge of each element ``|a\\rangle\\langle b|`` of a site, in the order of a vectorised
matrix: the flux of a one site operator made of that element alone
"""
function pair_charges(site::AbstractSite)
    c = basis_charges(site)
    return vec([ c[a] - c[b] for a in eachindex(c), b in eachindex(c) ])
end

"""
    total_charge(name, x, charges, sites, tol)

the charge the operator `name` carries, read off its elements above `tol`, and refused when
they disagree: such an operator could act on these sites neither whole nor split
"""
function total_charge(name, x, charges, sites, tol)
    q = common_charge(x, charges, tol)
    if isnothing(q)
        no_definite_charge(name, sites)
    end
    return q
end

"""
    along(f, y, j)

the linear map `f` applied along axis `j` of the array `y`
"""
function along(f::AbstractMatrix, y::AbstractArray, j::Int)
    p = [j; setdiff(1:ndims(y), j)]
    z = permutedims(y, p)
    w = reshape(f * reshape(z, size(z, 1), :), size(f, 1), size(z)[2:end]...)
    return permutedims(w, invperm(p))
end

"""
    part_on(x, d, on)

the part of an operator acting on the sites where `on` is true and as the identity on the
others, with an axis for each of the former. On those, the component along the identity is
projected out; the others are traced out and divided by their dimension, so that what the part
stands for there is the identity itself.
"""
function part_on(x::AbstractArray, d, on)
    y = x
    # from the last site, so that dropping an axis leaves those still to come in place
    for j in length(d):-1:1
        e = vec(identity_operator(d[j]))
        if on[j]
            y = along(I - e * transpose(e) / d[j], y, j)
        else
            y = dropdims(along(transpose(e) / d[j], y, j); dims = j)
        end
    end
    return y
end

"""
    svd_terms(y, charges, q, tol, make)

the terms of an operator acting on each of its sites, as pairs of a coefficient and of one
factor per site: `y` has an axis per site, `charges` gives the charge of each value of each
axis, `q` is the charge of the operator and `make(k, v)` builds the factor of its `k`-th site
from a vectorised matrix.

The first site is split from the others by a singular value decomposition, made charge by
charge so that every factor has a definite one, and each right singular vector is split in
turn in the same way.
"""
function svd_terms(y::AbstractArray, charges, q, tol, make)
    if ndims(y) == 1
        # what rounding left outside the charge of the operator is dropped, so that the
        # factor carries exactly that charge
        v = [ charges[1][p] == q ? y[p] : zero(eltype(y)) for p in eachindex(y) ]
        c = norm(v)
        return c ≤ tol ? [] : [ (c, [ make(1, v / c) ]) ]
    end
    rest = size(y)[2:end]
    rq = vec([ sum(charges[l + 1][i[l]] for l in eachindex(rest)) for i in CartesianIndices(rest) ])
    ym = reshape(y, size(y, 1), :)
    terms = []
    for c1 in unique(charges[1])
        rows = findall(==(c1), charges[1])
        cols = findall(==(q - c1), rq)
        if isempty(cols)
            continue
        end
        f = svd(ym[rows, cols])
        for k in eachindex(f.S)
            if f.S[k] ≤ tol
                break
            end
            v = zeros(eltype(f.Vt), length(rq))
            v[cols] = f.Vt[k, :]
            sub = svd_terms(reshape(v, rest), charges[2:end], q - c1, tol / f.S[k],
                            (l, w) -> make(l + 1, w))
            if isempty(sub)
                continue
            end
            u = zeros(eltype(f.U), size(y, 1))
            u[rows] = f.U[:, k]
            lead = make(1, u)
            for (c, fs) in sub
                push!(terms, (f.S[k] * c, [ lead; fs ]))
            end
        end
    end
    return terms
end

"""
    split_matrix(name, m, sites)

the definition `Operator{N}(name, def, type, sites...)` gives its operator, `m` being the
matrix of `def` on `sites`: a sum of tensor products of one site operators, each of a definite
charge. Each part acting on some of the sites and as the identity on the others is split on
its own, so that no term spans more sites than it acts on.
"""
function split_matrix(name::String, m::AbstractMatrix, sites::Vector)
    n = length(sites)
    d = [ dim(s) for s in sites ]
    # a real operator keeps real factors, which gives a real MPO and halves the cost of
    # every contraction with it
    if eltype(m) <: Complex && all(x -> iszero(imag(x)), m)
        m = real(m)
    end
    m = float(m)
    tol = rounding_tol * norm(m)
    check_even(name, m, sites, tol)
    # one axis per site, holding the vectorised matrix of a one site operator. The axes of a
    # matrix reshaped put the last site first and every output before every input, so each
    # site has its two brought together, the output varying faster
    x = reshape(permutedims(reshape(m, (reverse(d)..., reverse(d)...)),
                            [ k for j in 1:n for k in (n - j + 1, 2n - j + 1) ]),
                Tuple(d .^ 2))
    charges = [ pair_charges(s) for s in sites ]
    q = total_charge(name, x, charges, sites, tol)
    counts = zeros(Int, n)
    factor(j, v, k) = Operator{1}(name * script(superscripts, j) * script(subscripts, k),
                                  reshape(v, d[j], d[j]), plain_op)
    terms = GenericOp{Pure, n}[]
    for on in Iterators.product(fill((false, true), n)...)
        y = part_on(x, d, on)
        js = findall(collect(on))
        if isempty(js)
            if abs(y[]) > tol
                push!(terms, y[] * IdentityOp{Pure, Generic, n}())
            end
            continue
        end
        # a part acting on a single site is the one of index 0 there, the others are
        # numbered site by site in the order they come
        make = length(js) == 1 ? (l, v) -> factor(js[l], v, 0) :
                                 (l, v) -> factor(js[l], v, counts[js[l]] += 1)
        for (c, fs) in svd_terms(y, charges[js], q, tol, make)
            ops = GenericOp{Pure, 1}[ Id for _ in 1:n ]
            ops[js] = fs
            push!(terms, c * TensorOp{n}(ops))
        end
    end
    return SumOp(terms)
end

"""
    Operator{N}(name, def, type, sites...)

an operator of `N` sites whose definition `simplify` cannot develop (a matrix, a function of
its sites, or an expression such as `exp(X ⊗ X)`) split once and for all into a sum of tensor
products of one site operators. That sum becomes its definition, which `simplify` substitutes
as it does for `Swap`, and this is what lets it into a hamiltonian, a lindbladian or `expect`.
Created without its sites, such an operator can only be applied as a gate.

On a single site nothing is split: the definition is replaced by its matrix on that site,
computed once rather than each time a tensor is built.

The sites come in the order of the indices, a single one standing for `N` identical ones, and
a matrix is written in their basis with the last site varying fastest, as `matrix` gives it.
The split depends on the sites through:

- their dimensions, which the size of a matrix does not give when the sites differ;
- what they conserve: every factor carries a definite charge, so that the operator acts on
  these sites, and on the same sites once weakened;
- whether they are fermionic. A matrix is taken as it is, with no Jordan-Wigner string, which
  is only right for an operator commuting with `F` on each of its sites. Any other is refused,
  and is to be written as an expression of `C` and `dag(C)`, into which `simplify` inserts the
  strings.

The factors are named after the operator: `P2¹₂` is its second factor on its first site, and
the subscript 0 marks the part acting on that site alone. What acts as the identity on a site
is taken out first, so that no term spans more sites than it acts on, and the rest is split by
singular value decompositions.

# Examples

    P2 = Operator{2}("P2", m, selfadjoint_op, Spin(1))
    K = Operator{2}("K", mk, plain_op, Spin(1, conserve = 2Sz), Qubit(conserve = 2Sz))
    R = Operator{2}("R", exp(-0.3im * (X ⊗ X)), plain_op, Qubit())
"""
function Operator{N}(name::String, def::Union{Matrix, Function, GenericOp{Pure, N}},
                     type::OpType, site::AbstractSite, sites::AbstractSite...) where N
    ss = expand_sites(name, N, (site, sites...))
    m = checked_type(name, type, check_size(name, matrix(def, ss...), ss), ss)
    if N > 1
        return Operator{N}(name, split_matrix(name, m, ss), type)
    end
    # one site has nothing to be split: its matrix is only computed once and for all, and
    # its type, fermionic included, says what simplify does with it
    return Operator{1}(name, m, type)
end


############### Declaring and naming operators ###############

"""
    check_declared(site, declared)

check the operators `@def_operators` has just declared for `site`, given as `(name, type)`
pairs, against their types on that site, and `F` against being an involution. Other sites of
the same type are checked where the operators are used, see `checked_type`. A refusal takes
the declarations back, so that they can be made again once corrected.
"""
function check_declared(site::AbstractSite, declared)
    try
        for (name, type) in declared
            matrix(name == "F" ? F : Operator{1}(name, nothing, type), site)
        end
    catch
        for (name, _) in declared
            remove_definition(operator_definition, site, name)
        end
        rethrow()
    end
end

"""
    @def_operators(site, symbols)

define operators for `site`, grouped by `OpType` as in the example below.

Each name becomes a `const` of the calling module the first time it is seen. A name already in
scope is not bound again: it is registered for the new site and checked against what it
already stands for, so that reusing a name that stands for something else, or declaring it
with another `OpType`, is an error rather than a silent redefinition.

A definition is a matrix, an expression of operators of one site, or a function of the site
giving either. Every operator declared here acts on one site: one of several sites is defined
with `named` or `Operator{N}`. `F` declares the Jordan-Wigner operator of a fermionic site,
shared by every site, and binds no name of its own.

An operator neither fermionic nor self adjoint is `plain_op`. On a fermionic site an operator
is placed with no Jordan-Wigner string, so every type but `fermionic_op` has to commute with
`F`, and a `fermionic_op`, which moves a fermion, has to anticommute with it. The types, and
`F` being an involution, are checked on the site given, and again on each site an operator is
placed on.

# Examples

For a fermionic site type of your own, `MySite`, declared as the package declares `Fermion`:

    @def_operators(MySite(),
    [
        fermionic_op => 
        [
            C = [0. 1. ; 0. 0.],
        ],
        selfadjoint_op =>
        [
            N = dag(C) * C,
        ],
        involution_op =>
        [
            F = Float64[1 0 ; 0 -1]
        ]
    ])
"""
macro def_operators(site, symbols)
    e = Expr(:block)
    declared = []
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
            push!(declared, :(($nsym, $(esc(type)))))
            if nsym == "F"
                # `F` is the Jordan-Wigner operator of `Operators.jl`, shared by every
                # fermionic site and not an `Operator{1}`: the site is registered and the
                # name is left alone
                push!(e.args,
                quote
                    add_operator($__module__, $(esc(site)), $nsym, $(esc(val)), $(esc(type)))
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
                        add_operator($__module__, $(esc(site)), $nsym, $(esc(val)), $(esc(type)))
                    end)
            else
                push!(e.args,
                    quote
                        const $(esc(sym)) = add_operator($__module__, $(esc(site)), $nsym,
                                                         $(esc(val)), $(esc(type)))
                    end)
            end
        end
    end
    # once they are all declared, since a function of the site may use one declared after it
    push!(e.args, :(check_declared($(esc(site)), [$(declared...)])))
    return e
end

"""
    matrix_type(m[, sites, what])

the strongest `OpType` the matrix `m` has on `sites`, see `violation`: `involution_op`,
`selfadjoint_op`, `plain_op`, or `fermionic_op` when it anticommutes with the `F` of its
single site. A matrix of no definite parity on its site is refused, named `what`; one whose
size does not fit the sites is `plain_op`, to be refused where the operator is built.
"""
function matrix_type(m::AbstractMatrix, sites = AbstractSite[], what = "the matrix")
    n = isempty(sites) ? size(m, 1) : prod(dim, sites)
    # a size the sites do not fit is refused where the operator is built, naming it
    if size(m) ≠ (n, n)
        return plain_op
    end
    f = site_F(sites)
    for t in (involution_op, selfadjoint_op, plain_op, fermionic_op)
        if isnothing(violation(t, m, f))
            return t
        end
    end
    error("$what has no definite fermionic parity on $(only(sites)): write it with C and dag(C)")
end

function named(m::Matrix, name::String; type::OpType = matrix_type(m))
    # what simplify reads off the type, as an involution squaring to the identity, is checked
    # at once, where it was only checked once a matrix of the operator was computed, which
    # simplify could have spared. Parity needs the F of a site, where it is checked
    v = violation(type, m, nothing)
    if !isnothing(v)
        error("$name is declared $type but $v")
    end
    return Operator{1}(name, m, type)
end

named(f::Function, name::String; type::OpType = plain_op) =
    Operator{1}(name, f, type)

function named(def::Union{Matrix, Function, GenericOp{Pure}}, name::String,
               site::AbstractSite, sites::AbstractSite...; type::Union{Nothing, OpType} = nothing)
    n = named_sites(def, site, sites)
    ss = expand_sites(name, n, (site, sites...))
    t = isnothing(type) ? matrix_type(matrix(def, ss...), ss, name) : type
    return Operator{n}(name, def, t, ss...)
end


############### Conserved quantities ###############

"""
    site_charges(op, site)

the modulus and the charge of each basis state of `site` for the conserved quantity `op`, as
`(modulus, charges)`. A modulus of `1` is an additive charge over the integers, a modulus of
`m` a charge of the cyclic group of order `m`.

`op` has to be diagonal on the site, with eigenvalues that are either integers, taken as the
charges, or roots of unity, the charge being the exponent and the modulus read off the
denominators: this is what makes `Zd` work on a `Qudit` with no modulus written anywhere.

The two readings overlap on ±1 and are different conservations: two sites carrying -1 make -2
over the integers and 0 modulo 2. The integer reading wins, and `parity` asks for the other.
"""
function site_charges(op::SimpleOp, site::AbstractSite; tol::Float64 = charge_tol)
    # a renamed operator reads its charges off what it renames, which is where a modulus is
    # carried: named(parity(N), "P") read its ±1 as integers and conserved their sum instead
    # of a parity
    # and so does one declared for the site, whose definition is in the library rather than
    # in the operator: parity(N) given to @def_operators was read as a U(1) charge
    if op isa Operator
        def = isnothing(op.expr) ? operator_info(site, op.name) : op.expr
        if def isa Function
            def = def(site)
        end
        if def isa GenericOp
            return site_charges(def, site; tol)
        end
    end
    m = matrix(op, site)
    d = diag(m)
    # relative to the largest charge: the 80 bosons of Boson(80) missed theirs by 1.4e-14
    tol *= max(1, maximum(abs, d))
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
    flux(op, site)

the charge the one site operator `op` carries on `site`, as a `QN`: the difference between the
charges of the states it connects, taken modulo the modulus of the charge when it has one. It
is the zero charge for an operator commuting with everything the site conserves, and `QN()` on
a site conserving nothing.

An operator connecting states whose charges differ in more than one way has no flux and is
refused. `X` raises and lowers `N` at once, while under `parity(N)` both differences are 1
modulo 2. The modulus is also what gives `Xd` a flux under `Zd`: its wrap around, from the last
state to the first, is a difference of `1 - d`, which is 1 modulo `d`.

# Examples

    flux(N, Fermion(conserve = N))        # QN("N",0)
    flux(Sp, Qubit(conserve = 2Sz))       # QN("2Sz",2), in units of the declared charge
    flux(Xd, Qudit(3, conserve = Zd))     # QN("Zd",1,3)
"""
flux(op::SimpleOp, site::AbstractSite; tol::Float64 = rounding_tol) =
    charge_flux(matrix(op, site), op, site; tol)

flux(op::GenericOp{Pure}, site::AbstractSite) =
    error("flux is only defined for one site operators, and $op acts on several")

"""
    conserve_string(site, spec)

the string in which a site records what it conserves: for each quantity, its name, its
modulus when that is not 1, and the charge of every basis state. `spec` is what the user
wrote, one operator or a tuple of them, each possibly marked `strong`.

The operators are read here rather than kept, since they could not be written to a state file
and read back: the name a conserved quantity prints under, as `2Sz` or `parity(N)`, is an
expression and not the name of an operator declared for the site. Only their charges are
needed, and those are what the string holds.

A site type of your own calls it to fill its `conserve` field, the site being built bare
first, since the charges depend on its type and on its other fields, not on that one.

# Examples

    MySite(; conserve = ()) = MySite(conserve_string(MySite(""), conserve))

    conserve_string(Fermion(""), N)              # "N:0,1"
    conserve_string(Fermion(""), parity(N))      # "parity(N)%2:0,1"
    conserve_string(Electron(""), (Ntot, 2Sz))   # "2Sz:0,1,-1,0;Ntot:0,1,1,2"
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
        # the characters the recorded form and the relabelling of a strong symmetry read: a
        # name ending in %2 was read back as a charge modulo 2
        if endswith(name, '!') || endswith(name, '*') || any(in(name), (':', ';', '%'))
            error("cannot conserve $name: a name holding :, ; or %, or ending in ! or *, " *
                  "cannot be told from how a site records its charges")
        end
        # ITensors refuses a longer charge name when the index is built, far from here, and
        # the bra of a strong one takes a star (`ITensors.SmallStrings.smallLength`, internal)
        limit = ITensors.SmallStrings.smallLength - (spec isa Strong ? 1 : 0)
        if length(name) > limit
            error("cannot conserve $name: ITensors takes names of at most $limit characters " *
                  "here, give it a shorter one with named")
        end
        head = modulus == 1 ? name : "$name%$modulus"
        # a strong symmetry is marked on the quantity and not on the site, so that one site
        # may hold both kinds, and at the end of the head so that the name and the modulus
        # are read exactly as before
        return (spec isa Strong ? head * "!" : head) * ":" * join(q, ",")
    end
    # sorted by name, see `Conserved`
    return sorted_conserve(join(parts, ";"))
end
