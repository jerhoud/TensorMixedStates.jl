# From an operator to its matrix and its ITensor: the pairing of a ket and its bra into a mixed
# index, the laying of a matrix on the indices of its sites with its charges checked, and the
# tensor of an operator, generic or placed on a system.

export tensor, matrix


############### Between the pure and the mixed representations ###############

"""
    combinerto(i::Index, j::Index...)

a combiner gathering the indices `j` into the index `i`
"""
function combinerto(i::Index, j::Index...)
    c = combiner(j...; tags="")
    x = combinedind(c)
    replaceind(c, x, i)
end

"""
    mixer(j::Index, k::Index, site::AbstractSite)

the pair `(b, c)`: `b` the index carrying the bra of the ket index `j` on `site`, see
`bra_index`, and `c` the combiner gathering `j` and `dag(b')` into the mixed index `k`, the
bra daggered so that a density matrix has zero flux. Everything crossing between the pure and
the mixed representations goes through this combiner, so that the pairs are flattened in the
same order everywhere.
"""
function mixer(j::Index, k::Index, site::AbstractSite)
    b = bra_index(j, site)
    return b, combinerto(k, j, dag(b'))
end

"""
    ket_bras(system, i)

the function giving, for `(m, n)`, the element ``|m\\rangle\\langle n|`` of site `i`,
vectorised on its mixed index
"""
function ket_bras(system::System, i::Int)
    j = SysIndex{Pure}(system, i)
    b, c = mixer(j, SysIndex{Mixed}(system, i), system[i])
    return (m, n) -> onehot(j => m, dag(b') => n) * c
end

"""
    relabel(t, f)

the tensor `t` with each of its indices passed through `f`, its storage untouched
"""
relabel(t::ITensor, f) = ITensors.setinds(t, map(f, inds(t)))

"""
    relabeller(f)

a function relabelling an index by `f`, memoised on its identity and prime level, and giving
the result in the direction of the index it is given: the two ends of a link, one daggered,
must be relabelled to the same index for the tensors to keep contracting
"""
function relabeller(f)
    seen = Dict{Tuple{ITensors.IDType, Int}, Index}()
    return function (i::Index)
        x = get!(() -> f(i), seen, (id(i), plev(i)))
        return dir(x) == dir(i) ? x : dag(x)
    end
end

"""
    element_map(to, from, i, f; swap = false)

the tensor pairing each element ``|m\\rangle\\langle n|`` of site `i` of the system `from`,
passed through `f`, with the same element of `to`, or with ``|n\\rangle\\langle m|`` when
`swap` is true
"""
function element_map(to::System, from::System, i::Int, f; swap::Bool = false)
    d = dim(SysIndex{Pure}(from, i))
    to_element, from_element = ket_bras(to, i), ket_bras(from, i)
    return sum( to_element((swap ? (n, m) : (m, n))...) * dag(f(from_element(m, n)))
                for m in 1:d, n in 1:d )
end

"""
    element_maps(to, from, f; swap = false)

`element_map` on every site, built once for the sites alike and carried onto the indices of
the others: two sites whose pure indices have the same space in both systems, and which
conserve the same quantities strongly, give the same tensor but for its two indices
"""
function element_maps(to::System, from::System, f; swap::Bool = false)
    # the mixed index of site `i` in `to`, and the one `f` puts in place of that of `from`
    ends(i) = (SysIndex{Mixed}(to, i), only(inds(f(onehot(SysIndex{Mixed}(from, i) => 1)))))
    first_alike = Dict{Any, Int}()
    maps = ITensor[]
    for i in 1:length(to)
        key = (space(SysIndex{Pure}(to, i)), space(SysIndex{Pure}(from, i)),
               strong_names(to[i]), strong_names(from[i]))
        k = get!(first_alike, key, i)
        push!(maps, k == i ? element_map(to, from, i, f; swap) :
                             replaceinds(maps[k], ends(k), ends(i)))
    end
    return maps
end

"""
    weak_maps(strong, weak, relab)

the tensors carrying each site of the system `strong` onto the same site of `weak`, the
charged system it is weakened to, `relab` being the relabelling by `weak_index`: beyond
renaming, they regroup the blocks, which the two mixed indices order differently
"""
weak_maps(strong::System, weak::System, relab) =
    element_maps(weak, strong, t -> relabel(t, relab))

"""
    pure_map(old, new)

the tensor carrying the pure index `old` of a site onto `new`, basis state by basis state: the
pure counterpart of `weak_maps`, for two indices holding the basis in the same order and the
same charges but cut into other blocks
"""
pure_map(old::Index, new::Index) = sum(onehot(new => m) * dag(onehot(old => m)) for m in 1:dim(old))

"""
    adj_maps(system, relab)

the tensors sending ``|x\\rangle\\langle y|`` to ``|y\\rangle\\langle x|`` on each site,
from the index relabelled by `adjoint_index`, `relab` being that relabelling, to the index of
`system`. The relabelled index must be a new one, as it is on a charged system: on a plain
one the two sides would contract into a scalar, and a system without a strong symmetry uses
`tensor_dag` instead.
"""
adj_maps(system::System, relab) =
    element_maps(system, system, t -> relabel(t, relab); swap = true)

"""
    dense_maps(charged, plain)

the tensors carrying each site of the system `charged`, once densified, onto the same site of
`plain`, the same sites conserving nothing: densifying lays the blocks out in the order of
the charges, and each tensor is the permutation to the order of a plain combiner
"""
dense_maps(charged::System, plain::System) = element_maps(plain, charged, dense)


############### Laying a matrix on indices ###############

"""
    expand_sites(what, n, sites)
    expand_sites(a::GenericOp, sites)

the `n` sites `what` acts on, given as `sites`, `n` being that of the generic operator `a` in
the second form: a single site stands for `n` identical ones, and any other number of sites
than `n` is refused
"""
function expand_sites(what, n::Int, sites)
    ss = length(sites) == 1 ? fill(only(sites), n) : collect(AbstractSite, sites)
    if length(ss) ≠ n
        error("$what acts on $n sites and was given $(length(sites))")
    end
    return ss
end

expand_sites(a::GenericOp, sites) = expand_sites(a, nsites(a), sites)

"""
    named_sites(what, def, sites)

the sites of `def`, the definition of `what`, given with the sites `sites`, see `expand_sites`:
an expression acts on as many as it does, a function on those given, and a matrix given a
single site on as many copies of it as its size asks for
"""
named_sites(what, ::GenericOp{Pure, N}, sites) where N = expand_sites(what, N, sites)
named_sites(what, ::Function, sites) = expand_sites(what, length(sites), sites)
function named_sites(what, m::Matrix, sites)
    if length(sites) > 1
        return expand_sites(what, length(sites), sites)
    end
    site = only(sites)
    d, n, k = dim(site), size(m, 1), 1
    while d > 1 && d^k < n
        k += 1
    end
    if d^k ≠ n
        error("a $(size(m, 1))×$(size(m, 2)) matrix acts on no number of $site")
    end
    return expand_sites(what, k, sites)
end

"""
    site_indices(sites)

the indices of `sites` on their own, as `tensor` lays an operator on them: charged as soon as
one of them conserves something, a site conserving nothing then taking a trivial charge, as it
does in a system
"""
function site_indices(sites)
    charged = is_charged(sites)
    return [ site_index(s, charged) for s in sites ]
end

"""
    check_size(what, m, sites)

the matrix `m` of `what`, refused by a message naming `what` when its size is not the dimension
of `sites`
"""
function check_size(what, m::AbstractMatrix, sites)
    n = prod(dim, sites)
    if size(m) ≠ (n, n)
        error("$what is given by a $(size(m, 1))×$(size(m, 2)) matrix and cannot act on " *
              "$(join(sites, " ⊗ ")), whose dimension is $n")
    end
    return m
end

"""
    on_legs(m, outs, ins)

the matrix `m` of the combined space of several sites reshaped with one axis per index, and
those indices, `outs` then `ins`, each reversed: the array and the indices a tensor is built
from
"""
function on_legs(m::Matrix, outs, ins)
    idx = [ reverse(outs) ; reverse(ins) ]
    return (reshape(m, ntuple(k -> dim(idx[k]), length(idx))), idx)
end

"""
    lay(m, outs, ins)

the matrix `m` laid on the legs `outs` and `ins`, or `nothing` when their charges cannot
carry it, see `has_definite_flux`, for the caller to refuse it with `no_definite_charge`
"""
function lay(m::Matrix, outs, ins)
    a, idx = on_legs(m, outs, ins)
    return has_definite_flux(a, idx) ? charged_itensor(a, idx) : nothing
end

"""
    no_definite_charge(a, sites)

refuse the operator `a`, which carries no definite charge of what `sites` conserve
"""
no_definite_charge(a, sites) =
    error("$a carries no definite charge of " *
          "$(join(unique(q[1] for s in sites for q in decode_conserve(conserved(s))), ", ")), " *
          "so it cannot act on sites that conserve it")

"""
    pure_sides(js)

the outgoing and the incoming legs of an operator on a pure state, `j'` and `dag(j)` for each
index `j`
"""
pure_sides(js) = ([ j' for j in js ], [ dag(j) for j in js ])

"""
    mixed_sides(js, bs)

the outgoing and the incoming legs of an operator on a density matrix, in the order a matrix of
the mixed space reads them: the last site varying fastest and, inside a site, the ket faster
than the bra, as `mixer` pairs them
"""
mixed_sides(js, bs) =
    (reduce(vcat, [ [dag(b''), j'] for (j, b) in zip(js, bs) ]),
     reduce(vcat, [ [b', dag(j)] for (j, b) in zip(js, bs) ]))

"""
    fresh_bras(js, sites)

the new indices the bra of an operator on a density matrix lives on, one per site, see `legs`,
starred so that a site conserving something strongly keeps its bra apart from its ket
"""
fresh_bras(js, sites) = [ sim(star(j, strong_names(s))) for (j, s) in zip(js, sites) ]

"""
    onto_mixed(t, js, bs, ks)

the operator `t`, laid on the ket and bra legs of each site, carried onto `ks`, the mixed
indices of the sites, through the combiner of `mixer`
"""
function onto_mixed(t::ITensor, js, bs, ks)
    for (j, b, k) in zip(js, bs, ks)
        x = combinerto(k, j, dag(b'))
        t = t * dag(x) * x'
    end
    return t
end

"""
    combine_sites(t, is)

the operator `t`, laid on one pair of indices per site, gathered on a single pair combining
them: the form `tensor` gives for several sites
"""
function combine_sites(t::ITensor, is)
    if length(is) == 1
        return t
    end
    # on the daggered indices, so that the combiner meets the inputs and its primed dagger the
    # outputs
    c = combiner(reverse([ dag(i) for i in is ])...; tags = "")
    return t * c * dag(c')
end


############### The type of a matrix ###############

"""
    nearly(a, b)

whether `a` and `b` are equal up to rounding, relative to the larger of the two, see
`rounding_tol`
"""
nearly(a, b) = norm(a - b) ≤ rounding_tol * max(norm(a), norm(b))

"""
    violation(type, m, f)

how the matrix `m` fails to be of `type`, see `OpType`, as a phrase for an error message, or
`nothing` when it does not. `f` is the `F` of its site, or `nothing` when parity is not to be
checked.
"""
function violation(type::OpType, m::AbstractMatrix, f)
    if type in (selfadjoint_op, involution_op) && !nearly(m', m)
        return "is not self adjoint"
    elseif type == involution_op && !nearly(m * m, Matrix(I, size(m)...))
        return "its square is not the identity"
    elseif isnothing(f)
        return nothing
    elseif type == fermionic_op
        return nearly(f * m, -m * f) ? nothing : "does not anticommute with F"
    else
        return nearly(f * m, m * f) ? nothing : "does not commute with F"
    end
end


############### The matrix and the legs of an operator, which call each other ###############

"""
    site_F(sites)

the `F` of the site when `sites` holds a single one, and `nothing` otherwise: parity is only
checked on a single site, see `checked_type`
"""
site_F(sites) = length(sites) == 1 ? matrix(F, only(sites)) : nothing

"""
    checked_type(what, type, m, sites)

the matrix `m` of the operator `what`, declared of `type`, refused when it is not of that type,
see `OpType`, parity being checked only on a single site, see `split_matrix` for several. A
matrix whose size does not fit the sites is let through, to be refused where it is laid.
"""
function checked_type(what, type::OpType, m::AbstractMatrix, sites)
    n = prod(dim, sites)
    if size(m) ≠ (n, n)
        return m
    end
    v = violation(type, m, site_F(sites))
    if !isnothing(v)
        hint = v == "does not commute with F" ?
            ": declare it fermionic_op, or write it with C and dag(C)" : ""
        error("$what is declared $type but $v on $(join(sites, " ⊗ "))$hint")
    end
    return m
end

"""
    matrix(a::GenericOp, site::AbstractSite...)

the matrix of the generic operator `a` on the given sites, a single one standing for several
identical ones. It is written in the basis of the sites, whatever they conserve, the last site
varying fastest and, on a density matrix, the ket of a site faster than its bra. A tensor
product of fermionic factors holds the Jordan-Wigner strings between them: the matrix of
`C ⊗ dag(C)` is that of `C(1) * dag(C)(2)`, not the Kronecker product of the two matrices.

# Examples

    matrix(X, Qubit())
    matrix(Swap, Qubit())
    matrix(X⊗A, Qubit(), Boson(2))
    matrix(Left(X), Qubit())
"""
function matrix(a::Union{TensorOp, Left, Right}, site::AbstractSite...)
    sites = expand_sites(a, site)
    # plain indices: nothing is there to reorder the basis, and no charge to check
    js = [ Index(dim(s)) for s in sites ]
    if a isa GenericOp{Pure}
        t = legs(a, sites, js)
        outs, ins = pure_sides(js)
    else
        bs = [ sim(j) for j in js ]
        t = legs(a, sites, js, bs)
        outs, ins = mixed_sides(js, bs)
    end
    d = prod(dim, outs)
    # the reverse of `on_legs`
    return reshape(Array(t, reverse(outs)..., reverse(ins)...), d, d)
end

# a copy: changing it would change the declared operator
matrix(a::Matrix, ::AbstractSite, ::AbstractSite...) = copy(a)

matrix(a::Function, site::AbstractSite, sites::AbstractSite...) =
    matrix(a(site, sites...), site, sites...)

matrix(a::String, site::AbstractSite, ::AbstractSite...) =
    matrix(operator_info(site, a), site)

matrix(a::IdentityOp{Pure, Generic}, site::AbstractSite...) =
    identity_operator(prod(dim, expand_sites(a, site)))

matrix(a::IdentityOp{Mixed, Generic}, site::AbstractSite...) =
    identity_operator(prod(s -> dim(s)^2, expand_sites(a, site)))

function matrix(::JW_F, site::AbstractSite)
    m = matrix(F_info(site), site)
    if !isnothing(violation(involution_op, m, nothing))
        error("F is not an involution on $site: it has to be self adjoint and square to the identity")
    end
    return m
end

function matrix(a::Operator, site::AbstractSite...)
    # every site, which a definition by a function of the sites needs
    sites = expand_sites(a, site)
    m = isnothing(a.expr) ? matrix(a.name, site...) : matrix(a.expr, sites...)
    return checked_type(a, a.type, m, sites)
end

"""
    spectral_tol(values)

the distance below which two eigenvalues of an operator whose eigenvalues are `values` are one,
see `rounding_tol`
"""
spectral_tol(values) = rounding_tol * max(1, maximum(abs, values))

"""
    eigenspaces(a, site, what)
    eigenspaces(m, f, what, a, on = "")

the eigenvalues of the Hermitian operator `a` on `site`, or of its matrix `m`, increasing, and
the projectors on their eigenspaces, eigenvalues equal to rounding gathered into one, given by
their mean; `a` is refused when it is not Hermitian or does not commute with the fermionic
parity, `f` the matrix of the parity, `what` naming the caller in the refusals
"""
function eigenspaces(m::AbstractMatrix, f::AbstractMatrix, what::String, a, on::String = "")
    if !nearly(m, m')
        error("$what needs a Hermitian operator, which $a is not$on")
    end
    if !nearly(f * m, m * f)
        error("$what needs an operator commuting with the fermionic parity, which $a does " *
              "not$on")
    end
    e = eigen(Hermitian(complex(m)))
    tol = spectral_tol(e.values)
    d = size(m, 1)
    values = Float64[]
    ps = Matrix{ComplexF64}[]
    start = 1
    for k in 1:d
        if k == d || e.values[k+1] - e.values[start] > tol
            v = e.vectors[:, start:k]
            push!(values, sum(e.values[start:k]) / (k - start + 1))
            push!(ps, v * v')
            start = k + 1
        end
    end
    return values, ps
end

eigenspaces(a::SimpleOp, site::AbstractSite, what::String) =
    eigenspaces(matrix(a, site), matrix(F, site), what, a, " on $site")

"""
    eigenprojector(a, λ, sites...)

the projector on the eigenspace of eigenvalue `λ` of the Hermitian operator `a` on `sites`,
`λ` matching an eigenvalue to rounding, see `eigenspaces`
"""
function eigenprojector(a::GenericOp{Pure}, λ::Real, sites::AbstractSite...)
    on = join(sites, " ⊗ ")
    f = foldl(kron, [ matrix(F, s) for s in sites ])
    values, ps = eigenspaces(matrix(a, sites...), f, "Proj", a, " on $on")
    k = findfirst(v -> abs(v - λ) ≤ spectral_tol(values), values)
    if isnothing(k)
        error("$λ is not an eigenvalue of $a on $on, whose eigenvalues are $(join(values, ", "))")
    end
    return ps[k]
end

function matrix(a::Proj, given::AbstractSite...)
    sites = expand_sites(a, given)
    if a.state isa Pair
        return eigenprojector(a.state.first, a.state.second, sites...)
    end
    site = only(sites)
    st = state(site, a.state)
    if !(st isa Vector)
        error("Proj can only project on a pure state, \"$(a.state)\" is a mixed state")
    end
    # `simplify` takes a projector as even, commuting the F of the strings across it
    f = matrix(F, site)
    if f != I && !(nearly(f * st, st) || nearly(f * st, -st))
        error("$(repr(a.state)) has no definite fermionic parity on $site, which a projector " *
              "of a fermionic site needs")
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

matrix(a::IntPowOp, site::AbstractSite...) =
    matrix(a.arg, site...) ^ a.expo

function matrix(a::GenPowOp, site::AbstractSite...)
    m = matrix(a.arg, site...)
    r = rank(m)
    if real(a.expo) ≤ 0 && r < size(m, 1)
        error("$a does not exist: $(a.arg) is not invertible")
    elseif real(a.expo) > 0 && r ≠ rank(m * m)
        error("$a does not exist: the eigenvalue 0 of $(a.arg) is defective, as for C^0.5")
    end
    # Julia 1.10 refuses a non integer power of a real diagonal matrix with a negative entry,
    # and `Matrix` unwraps the Symmetric or Hermitian it may return
    if eltype(m) <: Real && isdiag(m) && any(<(0), diag(m))
        return Matrix(complex(m) ^ a.expo)
    end
    return Matrix(m ^ a.expo)
end

matrix(a::DagOp, site::AbstractSite...) =
    collect(adjoint(matrix(a.arg, site...)))

matrix(a::Gate, site::AbstractSite...) =
    matrix(Left(a.arg), site...) * matrix(Right(a.arg), site...)

function matrix(a::Dissipator, site::AbstractSite...)
    aa = dag(a.arg) * a.arg
    return matrix(Gate(a.arg), site...) - 0.5 * (matrix(Left(aa), site...) + matrix(Right(aa), site...))
end

# ρ ↦ m tr(ρ), the ket varying fastest; under a strong symmetry, where resetting a site moves
# the charge of one side only, `lay` refuses it
function matrix(a::SetState, site::AbstractSite)
    v = state(site, a.state)
    m = v isa Matrix ? v : v * v'
    return vec(m) * transpose(vec(identity_operator(site)))
end

# ρ ↦ Σ P ρ P over the eigenprojectors of the operator, see `eigenspaces`
matrix(a::Dephase, site::AbstractSite) =
    sum(kron(conj(p), p) for p in last(eigenspaces(a.arg, site, "Dephase")))

"""
    checked_matrix(a, sites, js)

the matrix of the operator `a` on `sites`, refused by a message naming `a` when its size is not
the dimension of the sites, or when it has no definite flux on a single site whose index is
charged
"""
function checked_matrix(a, sites, js)
    m = check_size(a, matrix(a, sites...), sites)
    if length(sites) == 1 && hasqns(only(js))
        charge_flux(m, a, only(sites))
    end
    return m
end

"""
    legs(a, sites, js)
    legs(a, sites, js, bs)

the tensor of the operator `a` on `sites`, `js` being their indices. An operator acting on a
pure state has one pair of legs per site, `j'` out and `dag(j)` in. One acting on a density
matrix has its bra on indices of its own, `bs`, drawn by `fresh_bras`: its legs are `j'` and
`dag(b'')` out, `dag(j)` and `b'` in, and `onto_mixed` then gathers each pair onto the mixed
index of its site.

No two sites are combined here, so no charge reorders the basis: only the combiner of `mixer`
sorts charged sectors.
"""
function legs(a::GenericOp{Pure}, sites, js)
    t = lay(checked_matrix(a, sites, js), pure_sides(js)...)
    if isnothing(t)
        no_definite_charge(a, sites)
    end
    return t
end

# (A₁ ⊗ … ⊗ Aₙ) on consecutive sites is A₁(1)…Aₙ(n): each factor takes an F per odd factor
# that follows, and a factor of no definite parity after a fermionic site is refused. Placed on
# sites apart, it misses the strings in between: `hasfermionic` sends it through `simplify`
function legs(a::TensorOp{Pure}, sites, js)
    pos = factor_sites(a)
    ps = map(jw_parity, a.subs)
    fs = GenericOp{Pure}[]
    for (k, o) in enumerate(a.subs)
        rest = ps[k+1:end]
        j = findfirst(isnothing, rest)
        if all(i -> matrix(F, sites[i]) == I, pos[k])
            push!(fs, o)
        elseif !isnothing(j)
            error("$a has a factor of no definite fermionic parity, $(a.subs[k+j]), after a " *
                  "fermionic site: write it as a sum of tensor products")
        elseif isodd(sum(rest; init = 0))
            push!(fs, o * reduce(⊗, fill(F, length(pos[k]))))
        else
            push!(fs, o)
        end
    end
    return prod(legs(o, sites[p], js[p]) for (o, p) in zip(fs, pos))
end

# each factor on its sites: no string, a tensor product holding a fermionic operator being
# expanded by `prepare_gate` and `simplify` before it is laid
legs(a::TensorOp{Mixed}, sites, js, bs) =
    prod(legs(o, sites[p], js[p], bs[p]) for (o, p) in zip(a.subs, factor_sites(a)))

# the identities are dense blocked because ITensors has no outer product of two charged deltas
legs(a::Left, sites, js, bs) =
    legs(a.arg, sites, js) * prod(denseblocks(delta(b', dag(b''))) for b in bs)

# the operator laid on the bras, `b'` out and `dag(b)` in, daggered and primed: its conjugate,
# `dag(b'')` out and `b'` in
legs(a::Right, sites, js, bs) =
    prod(denseblocks(delta(j', dag(j))) for j in js) * dag(legs(a.arg, sites, bs))'

function legs(a::GenericOp{Mixed}, sites, js, bs)
    m = matrix(a, sites...)
    t = lay(m, mixed_sides(js, bs)...)
    if !isnothing(t)
        return t
    end
    # the weak pairing, whose bra carries the charges of the ket, tells a jump that only the
    # strong symmetry forbids, one moving the charge, from one no conservation allows
    st = strong_names(sites)
    if !isempty(st) && !isnothing(lay(m, mixed_sides(js, [ sim(j) for j in js ])...))
        error("$a changes $(join(st, ", ")) between its two sides, which conserving it " *
              "strongly forbids: drop `strong` to allow a jump that moves the charge")
    end
    no_definite_charge(a, sites)
end


############### The tensor of an operator ###############

"""
    legs_on(a, sites, js, ks)

the tensor of the operator `a` on `sites`, whose pure indices are `js` and mixed ones `ks`: its
legs on `js` for an operator on pure states, gathered onto `ks` for one on density matrices,
see `legs`
"""
legs_on(a::GenericOp{Pure}, sites, js, _) = legs(a, sites, js)

function legs_on(a::GenericOp{Mixed}, sites, js, ks)
    bs = fresh_bras(js, sites)
    return onto_mixed(legs(a, sites, js, bs), js, bs, ks)
end

"""
    tensor(a::GenericOp, site::AbstractSite...)
    tensor(m::Matrix, site::AbstractSite...)

the ITensor of the generic operator `a`, or of the matrix `m`, on the given sites, a single
one standing for several identical ones, as many as the size of a matrix asks for. Several
sites are gathered on a single index combining theirs, and its primed form. The indices are
new ones, which contract with no state: `tensor(system, X(1))` lays a placed operator on the
indices of a system.

When a site conserves something, an operator or a matrix of no definite charge is refused. A
tensor product of fermionic factors holds the Jordan-Wigner strings between them, as in
`matrix`.

`tensor` of two operators is their tensor product, `tensor(X, Y)` being `X ⊗ Y`, and of two
systems their product, see `⊗`.

# Examples

    tensor(X, Qubit())
    tensor(Swap, Qubit())
    tensor(X⊗A, Qubit(), Boson(2))
    tensor([0. 1. ; 1. 0.], Qubit())
"""
function tensor(a::GenericOp, site::AbstractSite...)
    sites = expand_sites(a, site)
    js = site_indices(sites)
    ks = a isa GenericOp{Pure} ? js : [ mixed_index(j, s) for (j, s) in zip(js, sites) ]
    return combine_sites(legs_on(a, sites, js, ks), ks)
end

function tensor(a::Matrix, site::AbstractSite, sites::AbstractSite...)
    ss = named_sites("the matrix", a, (site, sites...))
    check_size("the operator", a, ss)
    js = site_indices(ss)
    t = lay(a, pure_sides(js)...)
    if isnothing(t)
        no_definite_charge("the matrix", ss)
    end
    return combine_sites(t, js)
end

"""
    tensor(system::System, a::AtIndex)

the ITensor of the operator `a` placed on sites, as `X(1)` or `(X ⊗ Y)(2, 3)`, on the indices
of `system`, with one pair of indices per site it acts on. It takes a single placed generic
operator, not a coefficient, as in `2X(1)`, nor a product, as in `X(1) * Y(2)`. Jordan-Wigner
strings are inserted only between the factors of a tensor product on consecutive sites: `C(3)`
gets the bare matrix of `C`, see `hasfermionic`.

# Examples

    tensor(System(3, Qubit()), X(1))
"""
function tensor(system::System, a::AtIndex{R}) where R
    sites = [ system[i] for i in a.index ]
    js = [ SysIndex{Pure}(system, i) for i in a.index ]
    ks = R === Pure ? js : [ SysIndex{Mixed}(system, i) for i in a.index ]
    return legs_on(a.op, sites, js, ks)
end
