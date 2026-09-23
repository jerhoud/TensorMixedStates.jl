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

# with a single index there is nothing to combine, and going through a combiner would not be
# harmless: it sorts the sectors it is given, while a site index keeps them in the order of
# the basis, one block per state. The two orders agree only when the charges of the basis
# happen to increase, which is why this was invisible until a site whose charges do not, an
# electron or a t-J, was asked for an observable.
#
# both are daggered because a combiner carries the indices it combines daggered and its own
# the other way, and `from` reaches here as it appears in the operator, which is already
# daggered
combinerto(to::Index, from::Index) = delta(dag(from), dag(to))

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
    vec_pieces(system, i, m)
    diag_elements(system, i)

the vectorised form of the one site matrix `m`, cut into pieces of definite charge, each
with the charge it carries. `diag_elements` is the identity, whose pieces are the diagonal
elements ``|x\\rangle\\langle x|`` and whose charges are all zero unless the site conserves
something strongly.

These are what a trace is made of, the trace of a density matrix being the sum of its
diagonal. When the site keeps its ket and its bra apart, each of them carries a charge of
its own and the sum has no single flux, which is why a trace then needs a chain rather than
one vector per site.
"""
function vec_pieces(system, i::Int, m::AbstractMatrix)
    j = SysIndex{Pure}(system, i)
    k = SysIndex{Mixed}(system, i)
    b, c = mixer(j, k, system[i])
    d = dim(j)
    # transposed, because the mixed index pairs the ket with the bra while a matrix is read
    # row by column: this is the convention `tensor_obs` has always used
    mt = transpose(m)
    acc = Dict{QN, ITensor}()
    for x in 1:d, y in 1:d
        if iszero(mt[x, y])
            continue
        end
        e = zeros(eltype(mt), d, d)
        e[x, y] = mt[x, y]
        t = dag(op_on_sites(e, [j], [dag(b')]) * c)
        g = flux(t)
        acc[g] = haskey(acc, g) ? acc[g] + t : t
    end
    return [ (t, g) for (g, t) in acc ]
end

diag_elements(system, i::Int) =
    vec_pieces(system, i, Matrix(1.0I, dim(system[i]), dim(system[i])))

"""
    adj_pieces(system, i)

the map sending ``|x\\rangle\\langle y|`` to ``|y\\rangle\\langle x|`` on site `i`, cut into
pieces of definite charge, each with the charge it carries.

Taking the adjoint of a density matrix exchanges its ket and its bra, which on a site that
keeps their charges apart swaps the two. That is a permutation of the blocks, and no single
tensor does it at a definite flux: the piece exchanging a pair whose charges differ by `d`
carries `d` on both sides. A chain along the running difference does, and since that
difference is zero over the whole state, it closes on nothing.
"""
function adj_pieces(system, i::Int)
    j = SysIndex{Pure}(system, i)
    k = SysIndex{Mixed}(system, i)
    b, c = mixer(j, k, system[i])
    d = dim(j)
    acc = Dict{QN, ITensor}()
    for x in 1:d, y in 1:d
        t = prime(ket_bra(system, i, y, x)) * dag(ket_bra(system, i, x, y))
        g = flux(t)
        acc[g] = haskey(acc, g) ? acc[g] + t : t
    end
    return [ (t, g) for (g, t) in acc ]
end

"""
    qn_list(i::Index)

the charges an index carries, one per block
"""
qn_list(i::Index) = [ first(p) for p in space(i) ]

"""
    charge_links(gs, total, tag)

the links a chain runs along, given the charges `gs[i]` its pieces may carry at each site and
the `total` they have to add up to.

Only what the sites on the left can have accumulated and what the ones on the right can still
bring is kept. That intersection is what makes a link small: most running totals cannot be
completed into the one the chain has to reach.
"""
function charge_links(gs, total::QN, tag::String)
    n = length(gs)
    forward = [[QN()]]
    for i in 1:n
        push!(forward, unique([ u - g for u in forward[i] for g in gs[i] ]))
    end
    backward = Vector{Vector{QN}}(undef, n + 1)
    backward[n+1] = [total]
    for i in n:-1:1
        backward[i] = unique([ v + g for v in backward[i+1] for g in gs[i] ])
    end
    return [ Index([ u => 1 for u in forward[i+1] if u in backward[i+1] ]...;
                   tags = "$tag,l=$i") for i in 1:n-1 ]
end

"""
    chain_at(links, ts, i, n, total)

the element at site `i` of a chain of `n` sites running along `links`, built from the pieces
`ts` and the charge each one carries. The last site closes the chain on `total`.
"""
function chain_at(links, ts, i::Int, n::Int, total::QN)
    ins = i == 1 ? [QN()] : qn_list(links[i-1])
    parts = ITensor[]
    for (t, g) in ts, (p, u) in enumerate(ins)
        left = i == 1 ? ITensor(1.) : onehot(dag(links[i-1]) => p)
        if i == n
            if u - g ≠ total
                continue
            end
            push!(parts, t * left)
        else
            r = findfirst(==(u - g), qn_list(links[i]))
            if isnothing(r)
                continue
            end
            push!(parts, t * left * onehot(links[i] => r))
        end
    end
    return sum(parts)
end

"""
    ket_bra(system, i, m, n)

the element ``|m\\rangle\\langle n|`` of site `i`, vectorised on its mixed index
"""
function ket_bra(system, i::Int, m::Int, n::Int)
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

Unlike the trace and the adjoint, this one is a plain local tensor and needs no chain. It
pairs the element ``|m\\rangle\\langle n|`` of one index with the same element of the other,
and the relabelling has already put the two under the same charge, so every term has a flux
of zero and their sum has one too. What it does beyond renaming is to gather the blocks the
relabelling left apart, which the two indices order differently.
"""
# the systems are left unannotated because this file is read before the one defining them
function weak_map(strong, weak, i::Int, relab)
    d = dim(SysIndex{Pure}(strong, i))
    return sum( ket_bra(weak, i, m, n) *
                dag(relabel(ket_bra(strong, i, m, n), relab))
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
    site_qns(system, i)

every charge the mixed index of site `i` can put on a vectorised one site tensor, as the
pieces of `vec_pieces` carry them. This is what the charge links of a trace have to be wide
enough for, an observable of non zero flux moving the running charge as it passes.
"""
site_qns(system, i::Int) = [ -first(p) for p in space(SysIndex{Mixed}(system, i)) ]

"""
    tensor_index(t::ITensor)

return the first index of an ITensor that is not primed
"""
tensor_index(t::ITensor) = getfirst(i->hasplev(i, 0), inds(t))

"""
    matrix(a::GenericOp, site::AbstractSite...)

return the matrix of a generic operator for the given sites. If sites are all identical, you may give only one

# Examples

    matrix(X, Qubit())
    matrix(Swap, Qubit())
    matrix(X⊗A, Qubit(), Boson(2))
"""
function matrix(a::GenericOp, site::AbstractSite...)
    t = tensor(a, site...)
    i = tensor_index(t)
    Matrix(t, i', i)
end

"""
    tensor(a::GenericOp, site::AbstractSite...)

return the ITensor of a generic operator for the given sites. If sites are all identical, you may give only one

# Examples

    tensor(X, Qubit())
    tensor(Swap, Qubit())
    tensor(X⊗A, Qubit(), Boson(2))
"""
function tensor(a::GenericOp, site::AbstractSite...; charged::Bool = false)
    m = matrix(a, site...)
    if length(site) == 1
        charge_flux(m, a, site[1])
    end
    return tensor(m, site...; charged)
end

"""
    op_index(sites, charged)

the index an operator of those sites lives on: theirs when there is one, the combination of
theirs otherwise.

Every index of the package comes from a site, so that two of them built for the same site are
interchangeable and a tensor computed here can be put on a system without a translation.
`charged` is the system's, a site conserving nothing taking a trivial index inside a system
where another one conserves. A matrix given for several identical sites at once, which is a
convenience of `matrix`, has no system behind it and stays dense.
"""
function op_index(sites, charged::Bool)
    is = [ site_index(s, charged) for s in sites ]
    return length(is) == 1 ? is[1] : combinedind(combiner(reverse(is)...; tags = ""))
end

function tensor(a::Matrix, site::AbstractSite, sites::AbstractSite...; charged::Bool = false)
    n, _ = size(a)
    if n ≠ dim(site) ^ (1 + length(sites))
        # the shorthand of one site standing for several identical ones, which names no
        # system and therefore carries no charge
        i = Index(n)
        return ITensor(a, i', dag(i))
    end
    # laid on the indices of the sites and only then combined, never written straight onto
    # the combination: combining charged indices sorts and merges their sectors, so the flat
    # order of the combined basis is not the order of the matrix
    is = [ site_index(s, charged) for s in (site, sites...) ]
    t = op_on_sites(a, [ i' for i in is ], [ dag(i) for i in is ])
    if length(is) == 1
        return t
    end
    c = combiner(reverse(is)...; tags = "")
    return t * c * dag(c')
end

matrix(a::Matrix, ::AbstractSite, ::AbstractSite...) = a

matrix(a::Function, site::AbstractSite, sites::AbstractSite...) =
    matrix(a(site, sites...), site, sites...)

matrix(a::String, site::AbstractSite, ::AbstractSite...) =
    matrix(operator_info(site, a), site)

matrix(::Identity, site::AbstractSite) =
    identity_operator(site)

matrix(::JW_F, site::AbstractSite) =
    matrix(F_info, site)    

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

matrix(a::PowOp, site::AbstractSite...) =
    matrix(a.arg, site...) ^ a.expo

matrix(a::DagOp, site::AbstractSite...) =
    collect(adjoint(matrix(a.arg, site...)))

function matrix(a::Dissipator, site::AbstractSite...)
    aa = dag(a.arg) * a.arg
    return matrix(Gate(a.arg), site...) - 0.5 * (matrix(Left(aa), site...) + matrix(Right(aa), site...))
end

matrix(a::Gate, site::AbstractSite...) =
    matrix(Left(a.arg), site...) * matrix(Right(a.arg), site...)

"""
    op_on_sites(m, outs, ins)

a matrix of the combined space of several sites, put back on their indices one by one.

The shape of the array is read off the indices themselves rather than from the dimensions of
the sites, which is the only way the two cannot disagree, and the order is the reversed one
the package uses everywhere it combines sites.
"""
function op_on_sites(m::Matrix, outs, ins)
    idx = [ reverse(outs) ; reverse(ins) ]
    return ITensor(reshape(m, ntuple(k -> dim(idx[k]), length(idx))), idx...)
end

"""
    super_tensor(m, is, left)

the tensor of the superoperator acting on one side of the density matrix by the matrix `m`,
on the left when `left` and on the right otherwise.

The density matrix carries a ket and a bra, and the mixed index of a site pairs them, so a
superoperator needs four slots while a mixed index and its primed form offer only three
distinct ones. The bra therefore lives on indices of its own, starred and drawn with `sim`, which leaves
room for the operator on one side and the identity on the other. ``\\rho \\mapsto A\\rho``
acts on the ket and leaves the bra alone, and its mirror ``\\rho \\mapsto \\rho
A^\\dagger`` does the reverse and conjugates.

Nothing here combines the sites into one index before splitting them again: a charged index
carries a direction, and going through a combined index is what no arrangement of `dag` could
be made to survive.
"""
function super_tensor(m::Matrix, is::Vector{<:Index}, sites, left::Bool)
    # starred first, so that a site conserving something strongly keeps its bra apart from
    # its ket, and `sim` then makes the fresh copy the four slots need. Starring nothing
    # gives the index back, and this is the `sim` it always was
    bs = [ sim(star(is[k], strong_names(sites[k]))) for k in eachindex(is) ]
    if left
        tk = op_on_sites(m, [ i' for i in is ], [ dag(i) for i in is ])
        tb = prod(delta(b', dag(b'')) for b in bs)
    else
        tk = prod(delta(i', dag(i)) for i in is)
        tb = op_on_sites(conj(m), [ dag(b'') for b in bs ], [ b' for b in bs ])
    end
    # the mixed index of a site pairs its ket with its bra, and the sites are combined
    # afterwards, which is the order a system builds its own indices in
    cs = [ combiner(is[k], dag(bs[k]'); tags = "") for k in eachindex(is) ]
    c = combiner(reverse(combinedind.(cs))...; tags = "")
    t = tk * tb
    for x in cs
        t = t * dag(x) * x'
    end
    return t * dag(c) * c'
end

"""
    tensor(a::GenericOp{Mixed}, site...)

the tensor of a superoperator whose matrix is already known, such as a `Gate` or a
`Dissipator`, both of which are written in terms of `Left` and `Right`. Only the index it
lives on has to be found, and that is the mixed one of its sites.
"""
function tensor(a::GenericOp{Mixed}, site::AbstractSite...; charged::Bool = false)
    is = [ site_index(s, charged) for s in site ]
    ms = [ mixed_index(is[k], site[k]) for k in eachindex(is) ]
    j = combinedind(combiner(reverse(ms)...; tags = ""))
    m = matrix(a, site...)
    if !hasqns(j) || all(s -> isempty(strong_names(s)), site)
        return ITensor(m, j', dag(j))
    end
    try
        return ITensor(m, j', dag(j))
    catch e
        if !(e isa ErrorException)
            rethrow()
        end
        # a superoperator of definite flux on the weak pairing may have none on the strong
        # one, which is exactly the case of a jump operator that moves the charge
        error("$a changes $(join(strong_names(site[1]), ", ")) between its two sides, which " *
              "conserving it strongly forbids: drop `strong` to allow a jump that moves " *
              "the charge")
    end
end

"""
    side_matrix(a, site)

the matrix of the operator a superoperator acts by, checked against the charges of its site
so that one carrying no flux is named here rather than deep inside ITensors
"""
function side_matrix(a, site::AbstractSite...)
    m = matrix(a, site...)
    if length(site) == 1
        charge_flux(m, a, site[1])
    end
    return m
end

tensor(a::Left, site::AbstractSite...; charged::Bool = false) =
    super_tensor(side_matrix(a.arg, site...), [ site_index(s, charged) for s in site ],
                 site, true)

tensor(a::Right, site::AbstractSite...; charged::Bool = false) =
    super_tensor(side_matrix(a.arg, site...), [ site_index(s, charged) for s in site ],
                 site, false)

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

function tensor(a::TensorOp{N}, site::AbstractSite...; charged::Bool = false) where N
    if length(site) == 1
        ts = [tensor(o, site...; charged) for o in a.subs]
    elseif length(site) ≠ N
        error("number of sites does not match operator")
    else
        ts = tensor_apply(tensor, a, site...; charged)
    end
    c = combiner((tensor_index(t) for t in reverse(ts))...; tags="")
    c * prod(ts) * c'
end

function tensor(a::SetState, site::AbstractSite; charged::Bool = false)
    i = site_index(site, charged)
    j = mixed_index(i, site)
    v = state(site, a.state)
    if v isa Matrix
        m = v
    else
        m = v * v'
    end
    # the target is laid on the two site indices and only then gathered, never written
    # straight onto the mixed one: combining charged indices merges and sorts their
    # sectors, so the flat order of the mixed basis is not the order of the matrix
    b, c = mixer(i, j, site)
    tr = denseblocks(delta(dag(i), b')) * dag(c)
    t = charged_state(() -> op_on_sites(m, [i'], [dag(b'')]), i, a.state, site)
    return tr * t * last(mixer(i', j', site))
end

