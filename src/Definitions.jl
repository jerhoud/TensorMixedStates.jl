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

refuse a matrix that does not commute with `F` on each of its sites: its one site factors are
placed with no Jordan-Wigner string and crossed by strings as if even, which is wrong for an
operator moving a fermion from one site to another
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
others, with an axis for each of the former: the component along the identity is projected
out on the former, the latter are traced out and divided by their dimension
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
from a vectorised matrix. The sites are split off one by one by singular value
decompositions, made charge by charge so that every factor has a definite one.
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
    # real factors for a real operator, which give a real MPO, cheaper to contract
    if eltype(m) <: Complex && all(x -> iszero(imag(x)), m)
        m = real(m)
    end
    m = float(m)
    tol = rounding_tol * norm(m)
    check_even(name, m, sites, tol)
    # one axis per site, holding the vectorised matrix of a one site operator: a matrix
    # reshaped puts the last site first and every output before every input, so the two axes
    # of each site are brought together, the output varying faster
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
its sites, or an expression such as `exp(X ⊗ X)`), split once and for all into a sum of tensor
products of one site operators, which `simplify` substitutes as it does for `Swap`: this lets
it into a Hamiltonian, a Lindbladian or `expect`. Created without its sites, such an operator
can only be applied as a gate. On a single site, the definition is replaced by its matrix.

The sites come in the order of the indices, a single one standing for `N` identical ones, and
a matrix is written in their basis with the last site varying fastest, as `matrix` gives it.
Every factor carries a definite charge of what the sites conserve. A matrix is taken with no
Jordan-Wigner string: one not commuting with `F` on each of its fermionic sites is refused, to
be written with `C` and `dag(C)`.

The factors are named after the operator: `P2¹₂` is its second factor on its first site, the
subscript 0 marking the part acting on that site alone. No term spans more sites than it acts
on.

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
    # its type, fermionic included, says what simplify does with it
    return Operator{1}(name, m, type)
end


############### Declaring and naming operators ###############

"""
    check_declared(site, declared)

check the operators `@def_operators` has just declared for `site`, given as `(name, type)`
pairs, against their types on that site, and `F` against being an involution, a refusal taking
the definitions back
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

define operators for `site`, grouped by `OpType` as in the example below. A definition is a
matrix, an expression of operators of one site, or a function of the site giving either; an
operator of several sites is defined with `named` or `Operator{N}`. `F` declares the
Jordan-Wigner operator of a fermionic site and binds no name of its own.

Each name becomes a `const` of the calling module the first time it is seen, so the macro is
used at the top level of a module or a script, not inside a function, a `let` or a
`@testset`. A name already in scope is registered for the new site, and refused if it stands
for something else or with another `OpType`.

An operator neither fermionic nor self adjoint is `plain_op`. On a fermionic site, a
`fermionic_op` has to anticommute with `F` and any other type to commute with it. The types
are checked on the site given and on each site an operator is placed on. A declaration refused
by its check is taken back and can be made again once corrected, but a name it bound keeps its
type until a new session. One interrupted by an error in its block is not taken back.

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
                # `F` is shared by every fermionic site and is not an `Operator{1}`: the site
                # is registered and the name left alone
                push!(e.args,
                quote
                    add_operator($__module__, $(esc(site)), $nsym, $(esc(val)), $(esc(type)))
                end)
            elseif isdefined(__module__, sym)
                # not bound again, which would rebind it for every site using it, and is an
                # error up to Julia 1.11 for a name brought in by `using`: decided at
                # expansion time, so that no binding is emitted
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
    # checked at once, simplify relying on the type; parity needs the F of a site, where it
    # is checked
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
    ss = named_sites(name, def, (site, sites...))
    n = length(ss)
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
denominators, as for `Zd` on a `Qudit`. On ±1, where both readings apply, the integer one
wins, and `parity` asks for the other.
"""
function site_charges(op::SimpleOp, site::AbstractSite; tol::Float64 = charge_tol)
    # a renamed operator, or one declared for the site, reads its charges off its definition,
    # which may carry a modulus, as parity(N) does
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
    # relative to the largest charge, whose rounding grows with it
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

# a power or the adjoint of a charge modulo m is one modulo m, whatever its eigenvalues: those
# of Zd^2 on Qudit(4), all ±1, would be read as integers
function site_charges(a::Union{IntPowOp{Pure, 1}, GenPowOp{Pure, 1}, DagOp{1}}, site::AbstractSite;
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
    flux(op, site)

the charge the one site operator `op` carries on `site`, as a `QN`: the difference between the
charges of the states it connects, taken modulo the modulus of the charge when it has one. It
is the zero charge for an operator commuting with everything the site conserves, and `QN()` on
a site conserving nothing.

An operator connecting states whose charges differ in more than one way has no flux and is
refused: `X` raises and lowers `N` at once, while under `parity(N)` both differences are 1
modulo 2.

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
modulus when that is not 1, and the charge of every basis state. `spec` is one operator or a
tuple of them, each possibly marked `strong`.

A site type of your own calls it to fill its `conserve` field, on the site built bare, as
below.

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
    names = [ obs_name(o isa Strong ? o.arg : o) for o in ops ]
    for n in names
        if count(==(n), names) > 1
            error("$n is given twice in conserve")
        end
    end
    parts = map(ops) do spec
        op = spec isa Strong ? spec.arg : spec
        modulus, q = site_charges(op, site)
        name = obs_name(op)
        if endswith(name, '!') || endswith(name, '*') || any(in(name), (':', ';', '%'))
            error("cannot conserve $name: a name holding :, ; or %, or ending in ! or *, " *
                  "cannot be told from how a site records its charges")
        end
        # the bra of a strong quantity takes a star in its name
        limit = ITensors.SmallStrings.smallLength - (spec isa Strong ? 1 : 0)
        if length(name) > limit
            error("cannot conserve $name: ITensors takes names of at most $limit characters " *
                  "here, give it a shorter one with named")
        end
        return encode_quantity(name, modulus, q, spec isa Strong)
    end
    # sorted by name, see `Conserved`
    return sorted_conserve(join(parts, ";"))
end
