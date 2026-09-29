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
`bra_index`, and `c` the combiner gathering `j` and `dag(b')` into the mixed index `k`.

The bra is daggered, so that without a strong symmetry the charge `k` carries is the
difference of those of the ket and the bra, which gives a density matrix zero flux. Under a
strong symmetry `b` carries its charges under starred names, and `k` keeps both apart.
Everything crossing between the pure and the mixed representations goes through this
combiner, so that the pairs are flattened in the same order everywhere.
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
    b, c = mixer(j, SysIndex{Mixed}(system, i), system[i])
    return onehot(j => m, dag(b') => n) * c
end

"""
    relabel(t, f)

the tensor `t` with each of its indices passed through `f`, its storage untouched, through
the internal `ITensors.setinds`
"""
relabel(t::ITensor, f) = ITensors.setinds(t, map(f, inds(t)))

"""
    relabeller(f)

a function relabelling an index by `f`, memoised on its identity and prime level, and giving
the result in the direction of the index it is given.

The key leaves the direction out because a link appears daggered on one of the two sites it
joins, and both ends have to be relabelled to the same index or the tensors stop contracting.
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
`swap` is true: the form `weak_map`, `adj_map` and `dense_map` share
"""
function element_map(to::System, from::System, i::Int, f; swap::Bool = false)
    d = dim(SysIndex{Pure}(from, i))
    return sum( ket_bra(to, i, (swap ? (n, m) : (m, n))...) * dag(f(ket_bra(from, i, m, n)))
                for m in 1:d, n in 1:d )
end

"""
    weak_map(strong, weak, i, relab)

the tensor carrying site `i` of the system `strong` onto the same site of `weak`, the charged
system it is weakened to, `relab` being the relabelling by `weak_index`.

It pairs each element ``|m\\rangle\\langle n|`` of one mixed index with the same element of
the other. The relabelling has already put both under the same charge, so every term has zero
flux; what the tensor does beyond renaming is to regroup the blocks, which the two indices
order differently.
"""
weak_map(strong::System, weak::System, i::Int, relab) =
    element_map(weak, strong, i, t -> relabel(t, relab))

"""
    adj_map(system, i, relab)

the tensor sending ``|x\\rangle\\langle y|`` to ``|y\\rangle\\langle x|`` on site `i`, from
the index relabelled by `adjoint_index`, `relab` being that relabelling, to the index of
`system`.

Under a strong symmetry the exchange alone has no definite flux, but the relabelling has
already given each element the charge of the one it is sent to, so every term has zero flux.
This needs the relabelled index to be a new one, which it always is on a charged system: on a
plain one the two sides would contract into a scalar, one reason, besides the cost, why a
system without a strong symmetry uses `tensor_dag` instead.
"""
adj_map(system::System, i::Int, relab) =
    element_map(system, system, i, t -> relabel(t, relab); swap = true)

"""
    dense_map(charged, plain, i)

the tensor carrying site `i` of the system `charged`, once densified, onto the same site of
`plain`, the same sites conserving nothing.

An ITensor holds either charged indices or plain ones, never both, so this last step of
weakening cannot be a relabelling: the state is densified first. Densifying lays the blocks
out in the order of the charges, not in the order a plain combiner gives, and this tensor is
the permutation between the two. Having no charges, it has no flux to respect.
"""
dense_map(charged::System, plain::System, i::Int) = element_map(plain, charged, i, dense)


############### Laying a matrix on indices ###############

"""
    expand_sites(what, n, sites)

the `n` sites `what` acts on, given as `sites`: a single site stands for `n` identical ones, and
any other number of sites than `n` is refused
"""
function expand_sites(what, n::Int, sites)
    ss = length(sites) == 1 ? fill(only(sites), n) : collect(AbstractSite, sites)
    if length(ss) ≠ n
        error("$what acts on $n sites and was given $(length(sites))")
    end
    return ss
end

"""
    all_sites(a, site)

the sites the operator `a` acts on, see `expand_sites`
"""
all_sites(a::GenericOp{R, N}, site) where {R, N} = expand_sites(a, N, site)

"""
    named_sites(def, site, sites)

the number of sites of a definition given with its sites: an expression has its own, a
function acts on those given, and a matrix given a single site acts on as many copies of it as
its size asks for
"""
named_sites(::GenericOp{Pure, N}, _, _) where N = N
named_sites(::Function, _, sites) = 1 + length(sites)
function named_sites(m::Matrix, site, sites)
    if !isempty(sites)
        return 1 + length(sites)
    end
    d, n, k = dim(site), size(m, 1), 1
    while d > 1 && d^k < n
        k += 1
    end
    if d^k ≠ n
        error("a $(size(m, 1))×$(size(m, 2)) matrix acts on no number of $site")
    end
    return k
end

"""
    site_indices(sites)

the indices of `sites` on their own, as `tensor` lays an operator on them: charged as soon as
one of them conserves something, a site conserving nothing then taking a trivial charge, as it
does in a system
"""
function site_indices(sites)
    charged = any(s -> !isempty(conserved(s)), sites)
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
carry it, see `has_definite_flux`. The caller then refuses the operator by name with
`no_definite_charge`, where ITensors would refuse the matrix with `Fluxes not all equal`, with
neither the operator nor its sites in sight.
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
than the bra, which is how `mixer` pairs them.
"""
mixed_sides(js, bs) =
    (reduce(vcat, [ [dag(b''), j'] for (j, b) in zip(js, bs) ]),
     reduce(vcat, [ [b', dag(j)] for (j, b) in zip(js, bs) ]))

"""
    fresh_bras(js, sites)

the indices the bra of an operator on a density matrix lives on, one per site, see `legs`:
starred, so that a site conserving something strongly keeps its bra apart from its ket, and
drawn with `sim`, so that they are new ones
"""
fresh_bras(js, sites) = [ sim(star(j, strong_names(s))) for (j, s) in zip(js, sites) ]

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
    combine_sites(t, is)

the operator `t`, laid on one pair of indices per site, gathered on a single pair combining
them: the form `tensor` gives for several sites. An operator placed on a system keeps one pair
per site and never goes through this.
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
see `OpType`. Parity is checked only on a single site; an operator of several sites has it
checked where it is split, see `split_matrix`. A matrix whose size does not fit the sites is
let through, to be refused where it is laid.
"""
function checked_type(what, type::OpType, m::AbstractMatrix, sites)
    n = prod(dim, sites)
    # a size the sites do not fit is refused by whoever lays the matrix, naming the operator
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

the matrix of the generic operator `a` on the given sites. When the sites are all identical,
a single one may be given for all of them.

The matrix is written in the basis of the sites, whatever they conserve. The last site varies
fastest and, for an operator acting on a density matrix, the ket of a site varies faster than
its bra.

# Examples

    matrix(X, Qubit())
    matrix(Swap, Qubit())
    matrix(X⊗A, Qubit(), Boson(2))
    matrix(Left(X), Qubit())
"""
function matrix(a::Union{TensorOp, Left, Right}, site::AbstractSite...)
    sites = all_sites(a, site)
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

matrix(a::Matrix, ::AbstractSite, ::AbstractSite...) = a

matrix(a::Function, site::AbstractSite, sites::AbstractSite...) =
    matrix(a(site, sites...), site, sites...)

matrix(a::String, site::AbstractSite, ::AbstractSite...) =
    matrix(operator_info(site, a), site)

matrix(a::IdentityOp{Pure, Generic}, site::AbstractSite...) =
    identity_operator(prod(dim, all_sites(a, site)))

# on a density matrix, the identity of the ket and the bra of each site
matrix(a::IdentityOp{Mixed, Generic}, site::AbstractSite...) =
    identity_operator(prod(s -> dim(s)^2, all_sites(a, site)))

function matrix(::JW_F, site::AbstractSite)
    m = matrix(F_info(site), site)
    if !isnothing(violation(involution_op, m, nothing))
        error("F is not an involution on $site: it has to be self adjoint and square to the identity")
    end
    return m
end

function matrix(a::Operator, site::AbstractSite...)
    # one site given for identical ones, which a function of the sites cannot take
    sites = all_sites(a, site)
    m = isnothing(a.expr) ? matrix(a.name, site...) : matrix(a.expr, sites...)
    return checked_type(a, a.type, m, sites)
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
    # Julia 1.10 takes a non integer power of a real diagonal matrix entry by entry and
    # refuses a negative entry, where later versions go complex; it also returns a
    # Symmetric or Hermitian wrapper for a non integer power of such a matrix, which the
    # rest of the package, laying matrices on indices, does not take
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

# the density matrix of the state times the trace: ρ ↦ m tr(ρ), the ket varying fastest. It is
# laid through `lay` as any operator of a density matrix, so that under a strong symmetry, where
# resetting a site moves the charge of one side only, it is refused rather than left with the
# block of charge zero alone, which gave a state of trace zero
function matrix(a::SetState, site::AbstractSite)
    v = state(site, a.state)
    m = v isa Matrix ? v : v * v'
    return vec(m) * transpose(vec(identity_operator(site)))
end

"""
    checked_matrix(a, sites, js)

the matrix of the operator `a` on `sites`, refused by a message naming `a` when its size is not
the dimension of the sites, or when it has no definite flux on a single site whose index is
charged. Otherwise the first would fail on a `DimensionMismatch` from `reshape` and the second
deep inside ITensors, naming neither the operator nor the site. The flux is only checked on a
charged index: `matrix` lays operators on plain indices, where there is none to check.
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
pure state has one pair of legs per site, `j'` out and `dag(j)` in.

An operator acting on a density matrix needs four legs per site, while a mixed index and its
primed form offer only three distinct ones, so its bra lives on indices of its own, `bs`,
drawn by `fresh_bras`. Its legs are `j'` and `dag(b'')` out, `dag(j)` and `b'` in, and
`onto_mixed` then gathers each pair onto the mixed index of its site.

Nothing here combines two sites: every matrix is laid on the indices of the sites one by one,
each keeping one block per basis state in the order of the basis, so no charge reorders
anything. The only thing that sorts charged sectors is the combiner of `mixer`, which is only
ever given a ket and its bra.
"""
function legs(a::GenericOp{Pure}, sites, js)
    t = lay(checked_matrix(a, sites, js), pure_sides(js)...)
    if isnothing(t)
        no_definite_charge(a, sites)
    end
    return t
end

# (A₁ ⊗ … ⊗ Aₙ) on consecutive sites is A₁(1)…Aₙ(n): the string of each odd factor, moved left
# through those before it, leaves on each of them an F per odd factor that follows, and a factor
# of no definite parity after a fermionic site makes the product a sum, which is refused rather
# than laid without its string. Sites apart, as a gate places it, still miss the strings of the
# sites in between, which is why `has_fermionic` sends such an operator through `simplify`
function legs(a::TensorOp, sites, js)
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

# the identities are dense blocked because ITensors has no outer product of two charged deltas
legs(a::Left, sites, js, bs) =
    legs(a.arg, sites, js) * prod(denseblocks(delta(b', dag(b''))) for b in bs)

# the operator laid on the bras, `b'` out and `dag(b)` in, daggered and primed: its conjugate,
# `dag(b'')` out and `b'` in. It goes through `legs` as the one of `Left` does, so that a
# factor of no definite charge is refused by name rather than by ITensors
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
    st = unique(reduce(vcat, [ strong_names(s) for s in sites ]))
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

the ITensor of the generic operator `a`, or of the matrix `m`, on the given sites. When the
sites are all identical, a single one may be given for all of them; for a matrix, their number
is then read off its size. Several sites are gathered on a single index combining theirs, and
its primed form.

When any of the sites conserves something, the indices carry charges, and an operator or a
matrix of no definite charge is refused.

# Examples

    tensor(X, Qubit())
    tensor(Swap, Qubit())
    tensor(X⊗A, Qubit(), Boson(2))
    tensor([0. 1. ; 1. 0.], Qubit())
"""
function tensor(a::GenericOp, site::AbstractSite...)
    sites = all_sites(a, site)
    js = site_indices(sites)
    ks = a isa GenericOp{Pure} ? js : [ mixed_index(j, s) for (j, s) in zip(js, sites) ]
    return combine_sites(legs_on(a, sites, js, ks), ks)
end

function tensor(a::Matrix, site::AbstractSite, sites::AbstractSite...)
    ss = expand_sites("the matrix", named_sites(a, site, sites), (site, sites...))
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
of `system`, with one pair of indices per site it acts on

# Examples

    tensor(System(3, Qubit()), X(1))
"""
function tensor(system::System, a::AtIndex{R}) where R
    sites = [ system[i] for i in a.index ]
    js = [ SysIndex{Pure}(system, i) for i in a.index ]
    ks = R === Pure ? js : [ SysIndex{Mixed}(system, i) for i in a.index ]
    return legs_on(a.op, sites, js, ks)
end
