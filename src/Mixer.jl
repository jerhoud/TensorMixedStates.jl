export tensor, matrix

"""
    combinerto(i::Index, j::Index...)

return a combiner that combines the indices `j` into the index `i`
"""
function combinerto(i::Index, j::Index...)
    c = combiner(j...; tags="")
    x = combinedind(c)
    replaceind(c, x, i)
end

"""
    mixer(j::Index, k::Index)

return the combiner that turns the ket index `j` and the bra index `j'` into the mixed
index `k`.

The bra is daggered, so that the charge `k` carries is the difference of the two and not
their sum. This is what makes a density matrix a tensor of zero flux, and it is the reason
everything that crosses between the two representations goes through this one combiner:
the pairs are then flattened in the same order everywhere. Without charges `dag` is the
identity and this is the plain combiner it always was.
"""
function mixer(j::Index, k::Index, site::AbstractSite)
    b = bra_index(j, site)
    return b, combinerto(k, j, dag(b'))
end


"""
    ket_bra(system, i, m, n)

the element ``|m\\rangle\\langle n|`` of site `i`, vectorised on its mixed index
"""
function ket_bra(system::System, i::Int, m::Int, n::Int)
    j = SysIndex{Pure}(system, i)
    k = SysIndex{Mixed}(system, i)
    b, c = mixer(j, k, system[i])
    e = zeros(dim(j), dim(j))
    e[m, n] = 1.
    return op_on_sites(e, [j], [dag(b')]) * c
end

"""
    weak_map(strong, weak, i, relab)

the tensor carrying site `i` of a strongly conserving system onto the same site of its
weakened one, `relab` being the relabelling of `weak_index`.

It pairs the element ``|m\\rangle\\langle n|`` of one index with the same element of the other,
and the relabelling has already put the two under the same charge, so every term has a flux
of zero and their sum has one too. What it does beyond renaming is to gather the blocks the
relabelling left apart, which the two indices order differently.
"""
function weak_map(strong::System, weak::System, i::Int, relab)
    d = dim(SysIndex{Pure}(strong, i))
    return sum( ket_bra(weak, i, m, n) *
                dag(relabel(ket_bra(strong, i, m, n), relab))
                for m in 1:d, n in 1:d )
end

"""
    adj_map(system, i, relab)

the tensor sending ``|x\\rangle\\langle y|`` to ``|y\\rangle\\langle x|`` on site `i`, from
the index relabelled by `adjoint_index` to the index of the system, `relab` being that
relabelling.

On its own the exchange has no definite flux under a strong symmetry, the two charges it
swaps being different ones. The relabelling has already given each element the charge of the
one it is sent to, so every term here has a flux of zero, as in `weak_map`.

This relies on the relabelled index being a new one. It always is on a charged system, but
on a plain one there is nothing to relabel and the two sides would contract into a scalar,
which is one reason, besides the cost, why a system without a strong symmetry keeps
`tensor_dag`.
"""
function adj_map(system::System, i::Int, relab)
    d = dim(SysIndex{Pure}(system, i))
    return sum( ket_bra(system, i, y, x) * dag(relabel(ket_bra(system, i, x, y), relab))
                for x in 1:d, y in 1:d )
end

"""
    dense_map(charged, dense, i)

the tensor carrying site `i` of a charged system onto the same site of one conserving
nothing.

The last rung of the ladder cannot be a relabelling: an ITensor holds either charged indices
or plain ones, never both, so the state is densified first. Densifying lays the blocks out in
the order the charges put them, which is not the order a plain combiner gives — measured, not
a single basis element of a mixed index lands in the same place. This tensor is the
permutation between the two, and having no charges it has no flux to respect.
"""
function dense_map(charged::System, plain::System, i::Int)
    d = dim(SysIndex{Pure}(charged, i))
    return sum( ket_bra(plain, i, m, n) * dag(dense(ket_bra(charged, i, m, n)))
                for m in 1:d, n in 1:d )
end

"""
    relabel(t, f)

the tensor `t` with each of its indices passed through `f`, the storage left as it is. This
leans on `ITensors.setinds`, which is not part of the public ITensors interface.
"""
relabel(t::ITensor, f) = ITensors.setinds(t, map(f, inds(t)))

"""
    relabeller(f)

a function relabelling an index by `f`, giving the same one back for the same identity and
prime level, and in the direction asked for.

Keyed on identity and prime level rather than on the index itself, because a link appears
daggered on one of the two sites it joins, and the two have to relabel to the same index or
the tensors stop contracting.
"""
function relabeller(f)
    seen = Dict{Tuple{ITensors.IDType, Int}, Index}()
    return function (i::Index)
        x = get!(() -> f(i), seen, (id(i), plev(i)))
        return dir(x) == dir(i) ? x : dag(x)
    end
end

"""
    all_sites(a, site)

the sites an operator acts on, a single site standing for as many identical ones as the
operator needs.
"""
function all_sites(a::GenericOp{R, N}, site) where {R, N}
    sites = length(site) == 1 ? fill(site[1], N) : collect(site)
    if length(sites) ≠ N
        error("$a acts on $N sites and was given $(length(site))")
    end
    return sites
end

"""
    matrix(a::GenericOp, site::AbstractSite...)

return the matrix of a generic operator for the given sites. If sites are all identical, you may give only one

The basis is the one of the sites, whatever they conserve: a matrix knows no charge. The last
site varies fastest and, for an operator acting on a density matrix, the ket of a site varies
faster than its bra.

# Examples

    matrix(X, Qubit())
    matrix(Swap, Qubit())
    matrix(X⊗A, Qubit(), Boson(2))
    matrix(Left(X), Qubit())
"""
function matrix(a::Union{TensorOp, Left, Right, SetState}, site::AbstractSite...)
    sites = all_sites(a, site)
    # plain indices: nothing is there to reorder the basis, and no charge to check
    js = [ Index(dim(s)) for s in sites ]
    if a isa GenericOp{Pure}
        t = legs(a, sites, js)
        outs, ins = [ j' for j in js ], js
    else
        bs = [ sim(j) for j in js ]
        t = legs(a, sites, js, bs)
        outs, ins = mixed_sides(js, bs)
    end
    d = prod(dim, outs)
    # the reverse of `op_on_sites`
    return reshape(Array(t, reverse(outs)..., reverse(ins)...), d, d)
end

"""
    tensor(a::GenericOp, site::AbstractSite...)

return the ITensor of a generic operator for the given sites. If sites are all identical, you may give only one

For several sites it lives on a single index combining theirs.

# Examples

    tensor(X, Qubit())
    tensor(Swap, Qubit())
    tensor(X⊗A, Qubit(), Boson(2))
"""
function tensor(a::GenericOp, site::AbstractSite...)
    sites = all_sites(a, site)
    # a site conserving nothing takes a trivial charge when another one conserves, as it does
    # in a system
    charged = any(s -> !isempty(conserved(s)), sites)
    js = [ site_index(s, charged) for s in sites ]
    if a isa GenericOp{Pure}
        return combine_sites(legs(a, sites, js), js)
    end
    bs = fresh_bras(js, sites)
    ks = [ mixed_index(j, s) for (j, s) in zip(js, sites) ]
    return combine_sites(onto_mixed(legs(a, sites, js, bs), js, bs, ks), ks)
end

"""
    combine_sites(t, is)

the operator `t`, laid on one pair of indices per site, gathered on a single pair combining
them, which is the form `tensor` gives for several sites. An operator placed on a system never
goes through this: it keeps one pair per site.
"""
function combine_sites(t::ITensor, is)
    if length(is) == 1
        return t
    end
    # built on the daggered indices, the ones the operator takes in, so that the combiner
    # carries them the other way and meets them; its primed dagger meets the outputs
    c = combiner(reverse([ dag(i) for i in is ])...; tags = "")
    return t * c * dag(c')
end

function tensor(a::Matrix, site::AbstractSite, sites::AbstractSite...)
    n, _ = size(a)
    if n ≠ prod(dim, (site, sites...))
        # the shorthand of one site standing for several identical ones, which names no
        # system and therefore carries no charge
        i = Index(n)
        return ITensor(a, i', dag(i))
    end
    ss = [site, sites...]
    charged = any(s -> !isempty(conserved(s)), ss)
    js = [ site_index(s, charged) for s in ss ]
    return combine_sites(op_on_sites(a, [ j' for j in js ], [ dag(j) for j in js ]), js)
end

matrix(a::Matrix, ::AbstractSite, ::AbstractSite...) = a

matrix(a::Function, site::AbstractSite, sites::AbstractSite...) =
    matrix(a(site, sites...), site, sites...)

matrix(a::String, site::AbstractSite, ::AbstractSite...) =
    matrix(operator_info(site, a), site)

matrix(::Identity, site::AbstractSite) =
    identity_operator(site)

matrix(::JW_F, site::AbstractSite) =
    matrix(F_info(site), site)

matrix(a::Operator, site::AbstractSite...) =
    if isnothing(a.expr)
        matrix(a.name, site...)
    else
        matrix(a.expr, site...)
    end

function matrix(a::Proj, site::AbstractSite, ::AbstractSite...)
    st = state(site, a.state)
    if !(st isa Vector)
        error("Proj can only project on a pure state, \"$(a.state)\" is a mixed state")
    end
    return st * adjoint(st)
end

matrix(a::JW, site::AbstractSite) =
    matrix(a.arg, site)

matrix(a::ScalarOp, site::AbstractSite...) =
    a.coef * matrix(a.arg, site...)

matrix(a::ProdOp, site::AbstractSite...) =
    prod(matrix(o, site...) for o in a.subs)

matrix(a::SumOp, site::AbstractSite...) =
    sum(matrix(o, site...) for o in a.subs)

matrix(a::ExpOp, site::AbstractSite...) =
    exp(matrix(a.arg, site...))

matrix(a::ModOp, site::AbstractSite...) =
    exp(2im * π * matrix(a.arg, site...) / a.modulus)

function matrix(a::PowOp, site::AbstractSite...)
    m = matrix(a.arg, site...)
    # Julia 1.10 takes a non integer power of a real diagonal matrix entry by entry and
    # refuses a negative entry, where later versions go complex
    if !isinteger(a.expo) && eltype(m) <: Real && isdiag(m) && any(<(0), diag(m))
        return complex(m) ^ a.expo
    end
    return m ^ a.expo
end

matrix(a::DagOp, site::AbstractSite...) =
    collect(adjoint(matrix(a.arg, site...)))

function matrix(a::Dissipator, site::AbstractSite...)
    aa = dag(a.arg) * a.arg
    return matrix(Gate(a.arg), site...) - 0.5 * (matrix(Left(aa), site...) + matrix(Right(aa), site...))
end

matrix(a::Gate, site::AbstractSite...) =
    matrix(Left(a.arg), site...) * matrix(Right(a.arg), site...)


############### Operators of several sites split into one site factors ###############

"""
    Operator{N}(name, def, type, sites...)

an operator of several sites whose definition `simplify` cannot develop, a matrix, a function
of its sites or an expression such as `exp(X ⊗ X)`, split once and for all into a sum of
products of one site operators. That sum becomes its definition, which `simplify` replaces it
with as it does for `Swap`, and this is what lets it into a hamiltonian, a lindbladian or
`expect`. Given without its sites, such an operator can only be applied as a gate.

On a single site there is nothing to split, and the definition is only replaced by its
matrix on that site, computed once rather than each time a tensor is built.

The sites come in the order of the indices, a single one standing for as many identical ones,
and a matrix is written in their basis with the last site varying fastest, as `matrix` gives
it. The factors are computed from them:

- the dimension of each site, which the size of a matrix does not give when the sites differ;
- what each site conserves: every factor carries a definite charge, so that the operator acts
  on these sites, and on the same sites once weakened;
- whether a site is fermionic. A matrix is taken as it is, with no Jordan-Wigner string, and
  that is only right for an operator commuting with `F` on each of its sites. One that does
  not is refused, and is to be written as an expression of `C` and `dag(C)`, into which
  `simplify` inserts the strings.

The factors are named after the operator, `P2¹₂` being its second factor on its first site
and the index 0 the part acting on that site alone. What acts as the identity on a site is
taken out first, so that a term acting on a single site takes no channel of an MPO, and the
rest is split by singular value decompositions, which give the fewest terms across each link.

# Examples

    P2 = Operator{2}("P2", m, selfadjoint_op, Spin(1))
    K = Operator{2}("K", mk, plain_op, Spin(1, conserve = 2Sz), Qubit(conserve = 2Sz))
    R = Operator{2}("R", exp(-0.3im * (X ⊗ X)), plain_op, Qubit())
"""
function Operator{N}(name::String, def::Union{Matrix, Function, GenericOp{Pure, N}},
                     type::OpType, site::AbstractSite, sites::AbstractSite...) where N
    ss = isempty(sites) ? fill(site, N) : AbstractSite[site, sites...]
    if length(ss) ≠ N
        error("$name acts on $N sites and was given $(length(ss))")
    end
    m = matrix(def, ss...)
    if N > 1
        return Operator{N}(name, split_matrix(name, m, ss), type)
    end
    # one site has nothing to be split: its matrix is only computed once and for all, and
    # its type, fermionic included, says what simplify does with it
    d = dim(only(ss))
    if size(m) ≠ (d, d)
        error("$name is given by a $(size(m, 1))×$(size(m, 2)) matrix, but $(only(ss)) " *
              "has dimension $d")
    end
    return Operator{1}(name, m, type)
end

const superscripts = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")
const subscripts = collect("₀₁₂₃₄₅₆₇₈₉")

# a number written with the given digits
script(digits, n::Int) = join(digits[c - '0' + 1] for c in string(n))

"""
    split_matrix(name, m, sites)

the definition `Operator{N}(name, def, type, sites...)` gives its operator, `m` being the
matrix of `def` on `sites`: a sum of tensor products of one site operators, each of a definite
charge.

Every part acting on some of the sites and as the identity on the others is split on its own,
which is what keeps a term from spanning more sites than it acts on.
"""
function split_matrix(name::String, m::AbstractMatrix, sites::Vector)
    n = length(sites)
    d = [ dim(s) for s in sites ]
    D = prod(d)
    if size(m) ≠ (D, D)
        error("$name is given by a $(size(m, 1))×$(size(m, 2)) matrix, but " *
              "$(join(sites, " ⊗ ")) has dimension $D")
    end
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
                push!(terms, y[] * TensorOp{n}(fill(Id, n)))
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
    check_even(name, m, sites, tol)

refuse a matrix that does not commute with `F` on each of its sites. A matrix is taken as it
is, with no Jordan-Wigner string between its sites, and the strings of the other factors of a
product cross its factors as if they were even: that is right for a density, a spin or a pair,
all even on each site, and not for an operator moving a fermion from one site to another.
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
            error("$name does not commute with F on its site $j, $s, so it moves a fermion " *
                  "there. A matrix is taken as it is, with no Jordan-Wigner string, which " *
                  "is only right for an operator even on each of its sites: write it as an " *
                  "expression of C and dag(C) instead, into which simplify inserts the strings")
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
    qs = decode_conserve(conserved(site))
    c = [ QN([ (name, ch[k], modulus) for (name, modulus, ch, _) in qs ]...) for k in 1:dim(site) ]
    return vec([ c[a] - c[b] for a in eachindex(c), b in eachindex(c) ])
end

"""
    total_charge(name, x, charges, sites, tol)

the charge an operator carries, read off its elements, refused when they do not agree: such an
operator could act on these sites neither whole nor split
"""
function total_charge(name, x, charges, sites, tol)
    q = nothing
    for i in CartesianIndices(x)
        if abs(x[i]) ≤ tol
            continue
        end
        c = sum(charges[j][i[j]] for j in eachindex(charges))
        if isnothing(q)
            q = c
        elseif c ≠ q
            no_definite_charge(name, sites)
        end
    end
    return isnothing(q) ? QN() : q
end

"""
    part_on(x, d, on)

the part of an operator acting on the sites where `on` is true and as the identity on the
others, with an axis for each of the former: what acts as the identity is taken out on each of
them, and the others are traced out, divided by their dimension so that what the part stands
for there is the identity itself.
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

# the linear map f applied along axis j of y
function along(f::AbstractMatrix, y::AbstractArray, j::Int)
    p = [j; setdiff(1:ndims(y), j)]
    z = permutedims(y, p)
    w = reshape(f * reshape(z, size(z, 1), :), size(f, 1), size(z)[2:end]...)
    return permutedims(w, invperm(p))
end

"""
    svd_terms(y, charges, q, tol, make)

the terms of an operator acting on each of its sites, as pairs of a coefficient and of one
factor per site: `y` has an axis per site, `q` is its charge and `make(k, v)` builds the factor
of its `k`-th site from a vectorised matrix.

The first site is split from the others by a singular value decomposition, made charge by
charge so that every factor has a definite one, and what it leaves on the others is split in
the same way.
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
    op_on_sites(m, outs, ins)

a matrix of the combined space of several sites, put back on their indices one by one.

The shape of the array is read off the indices themselves rather than from the dimensions of
the sites, which is the only way the two cannot disagree, and the order is the reversed one
the package uses everywhere it combines sites.
"""
function op_on_sites(m::Matrix, outs, ins)
    idx = [ reverse(outs) ; reverse(ins) ]
    return charged_itensor(reshape(m, ntuple(k -> dim(idx[k]), length(idx))), idx)
end

"""
    legs(a, sites, js)
    legs(a, sites, js, bs)

the tensor of an operator on the given sites, with one pair of legs per site: `j'` and
`dag(j)` for an operator acting on a pure state, `js` being the indices of the sites.

An operator acting on a density matrix needs four legs per site, while a mixed index and its
primed form offer only three distinct ones, so its bra lives on indices of its own, `bs`, drawn
by `fresh_bras`. Its legs are `j'` and `dag(b'')` on the way out, `dag(j)` and `b'` on the way
in, and `onto_mixed` gathers each pair onto the mixed index of its site.

Nothing here combines two sites. Every matrix is laid on the indices of the sites one by one,
and the index of a site keeps one block per basis state in the order of the basis, so no charge
reorders anything. The one thing that knows how charged sectors are sorted is the combiner of
`mixer`, and it is only ever handed a ket and its bra.
"""
function legs(a::GenericOp{Pure}, sites, js)
    t = lay(checked_matrix(a, sites, js), [ j' for j in js ], [ dag(j) for j in js ])
    if isnothing(t)
        no_definite_charge(a, sites)
    end
    return t
end

legs(a::TensorOp{N}, sites, js) where N =
    prod(tensor_apply((o, p...) -> legs(o, sites[[p...]], js[[p...]]), a, (1:N)...))

# the identities are dense blocked because ITensors has no outer product of two charged deltas
legs(a::Left, sites, js, bs) =
    legs(a.arg, sites, js) * prod(denseblocks(delta(b', dag(b''))) for b in bs)

# the operator laid on the bras, `b'` out and `dag(b)` in, daggered and primed: its conjugate,
# `dag(b'')` out and `b'` in. It goes through `legs` as the one of `Left` does, so that a
# factor of no definite charge is refused by name rather than by ITensors
legs(a::Right, sites, js, bs) =
    prod(denseblocks(delta(j', dag(j))) for j in js) * dag(legs(a.arg, sites, bs))'

function legs(a::SetState, sites, js, bs)
    site, j, b = only(sites), only(js), only(bs)
    v = state(site, a.state)
    m = v isa Matrix ? v : v * v'
    return denseblocks(delta(dag(j), b')) *
           charged_state(() -> op_on_sites(m, [j'], [dag(b'')]), j, a.state, site)
end

function legs(a::GenericOp{Mixed}, sites, js, bs)
    m = matrix(a, sites...)
    t = lay(m, mixed_sides(js, bs)...)
    if !isnothing(t)
        return t
    end
    # the weak pairing, whose bra carries the charges of the ket, tells a jump that only the
    # strong symmetry forbids, one moving the charge, from one no conservation allows
    st = unique(reduce(vcat, [ strong_names(s) for s in sites ]))
    if !isempty(st) && !isnothing(lay(m, mixed_sides(js, [ sim(j) for j in js ])...))
        error("$a changes $(join(st, ", ")) between its two sides, which conserving it " *
              "strongly forbids: drop `strong` to allow a jump that moves the charge")
    end
    no_definite_charge(a, sites)
end

"""
    lay(m, outs, ins)
    no_definite_charge(a, sites)

the matrix `m` laid on the given legs, or `nothing` when their charges cannot carry it, and
the refusal of the operator it came from. ITensors refuses such a matrix with `Fluxes not
all equal`, from a place where neither the operator nor its sites are in sight, so the
question is asked here and the answer given in terms of the operator.
"""
function lay(m::Matrix, outs, ins)
    try
        return op_on_sites(m, outs, ins)
    catch e
        if !(e isa ErrorException)
            rethrow()
        end
        return nothing
    end
end

no_definite_charge(a, sites) =
    error("$a carries no definite charge of " *
          "$(join(unique(q[1] for s in sites for q in decode_conserve(conserved(s))), ", ")), " *
          "so it cannot act on sites that conserve it")

"""
    checked_matrix(a, sites, js)

the matrix of an operator, refused by a message naming it when it carries no definite flux on
a single site whose index is charged, rather than deep inside ITensors. The check is made only
there: a matrix knows no charge, and `matrix` lays it on plain indices.

A matrix whose size is not the dimension of its sites is refused here as well. A matrix knows
no site either, and it failed on a `DimensionMismatch` from `reshape`, which named neither the
operator nor the site.
"""
function checked_matrix(a, sites, js)
    m = matrix(a, sites...)
    n = prod(dim, sites)
    if size(m) ≠ (n, n)
        error("$a is given by a $(size(m, 1))×$(size(m, 2)) matrix and cannot act on " *
              "$(join(sites, " ⊗ ")), whose dimension is $n")
    end
    if length(sites) == 1 && hasqns(only(js))
        charge_flux(m, a, only(sites))
    end
    return m
end

"""
    fresh_bras(js, sites)

the indices the bra of an operator on a density matrix lives on, one per site: starred, so that
a site conserving something strongly keeps its bra apart from its ket, and drawn with `sim` for
the fourth slot. Starring nothing gives the index back, and this is the `sim` it always was.
"""
fresh_bras(js, sites) = [ sim(star(j, strong_names(s))) for (j, s) in zip(js, sites) ]

"""
    mixed_sides(js, bs)

the outgoing and the incoming legs of an operator on a density matrix, in the order a matrix of
the mixed space reads them: the last site varying fastest and, inside a site, the ket faster
than the bra, which is how `mixer` pairs them.
"""
mixed_sides(js, bs) =
    (reduce(vcat, [ [dag(b''), j'] for (j, b) in zip(js, bs) ]),
     reduce(vcat, [ [b', dag(j)] for (j, b) in zip(js, bs) ]))

"""
    onto_mixed(t, js, bs, ks)

the operator `t`, laid on the ket and bra legs of each site, carried onto `ks`, the mixed
indices of the sites, through the combiner of `mixer`.
"""
function onto_mixed(t::ITensor, js, bs, ks)
    for (j, b, k) in zip(js, bs, ks)
        x = combinerto(k, j, dag(b'))
        t = t * dag(x) * x'
    end
    return t
end

"""
    tensor(::System, ::AtIndex)

returns a tensor representing the given simple indexed operator acting on this system

It is built on the indices of the system, one pair per site, so there is nothing to carry over
afterwards. See `legs`.
"""
function tensor(system::System, a::AtIndex{R}) where R
    sites = [ system[i] for i in a.index ]
    js = [ SysIndex{Pure}(system, i) for i in a.index ]
    if R === Pure
        return legs(a.op, sites, js)
    end
    bs = fresh_bras(js, sites)
    return onto_mixed(legs(a.op, sites, js, bs), js, bs,
                      [ SysIndex{Mixed}(system, i) for i in a.index ])
end

tensor_next(f, o::GenericOp{Pure, N}, site::Vararg{Union{AbstractSite, Int}, M}; kwargs...) where {N, M} =
    (f(o, site[1:N]...; kwargs...), site[N+1:M])

function tensor_apply(f, a::TensorOp{N}, idx::Vararg{Union{AbstractSite, Int}, N}; kwargs...) where N
    rest = idx
    r = map(a.subs) do o
        t, rest = tensor_next(f, o, rest...; kwargs...)
        t
    end
end

tensor_apply(::Any, ::TensorOp; kwargs...) = error("bug: tensor_apply")
tensor_next(::Any, ::GenericOp{Pure, N}; kwargs...) where N = error("bug: tensor_next")

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
    # a renamed operator reads its charges off what it renames, which is where a modulus is
    # carried: named(parity(N), "P") read its ±1 as integers and conserved their sum instead
    # of a parity
    if op isa Operator && op.expr isa GenericOp
        return site_charges(op.expr, site; tol)
    end
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
flux(op::SimpleOp, site::AbstractSite; tol::Float64 = rounding_tol) =
    charge_flux(matrix(op, site), op, site; tol)

flux(op::GenericOp{Pure}, site::AbstractSite) =
    error("flux is only defined for one site operators, and $op acts on several")

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
    return join(parts, ";")
end
