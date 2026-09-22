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
    to_bra(t, from, to)

the tensor `t` with its index `from` carried by `to` instead.

The two indices hold the same blocks in the same order and differ only in the names their
charges go under, so the storage is valid as it stands and only the index set is swapped.
This leans on `ITensors.setinds`, which is not part of the public ITensors interface: an
upgrade that moves it breaks here. Nothing to do when the site conserves nothing strongly,
the bra being the ket index itself.
"""
to_bra(t::ITensor, from::Index, to::Index) =
    from === to ? t : ITensors.setinds(t, map(i -> i == from ? to : i, inds(t)))

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
    if n == dim(site) ^ (1 + length(sites))
        i = op_index((site, sites...), charged)
    else
        # the shorthand of one site standing for several identical ones, which names no
        # system and therefore carries no charge
        i = Index(n)
    end
    ITensor(a, i', dag(i))
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
    ms = [ mix(is[k], site[k]) for k in eachindex(is) ]
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
    j = mix(i, site)
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

