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
mixer(j::Index, k::Index) = combinerto(k, j, dag(j'))

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
distinct ones. The bra therefore lives on indices of its own, drawn with `sim`, which leaves
room for the operator on one side and the identity on the other. ``\\rho \\mapsto A\\rho``
acts on the ket and leaves the bra alone, and its mirror ``\\rho \\mapsto \\rho
A^\\dagger`` does the reverse and conjugates.

Nothing here combines the sites into one index before splitting them again: a charged index
carries a direction, and going through a combined index is what no arrangement of `dag` could
be made to survive.
"""
function super_tensor(m::Matrix, is::Vector{<:Index}, left::Bool)
    bs = [ sim(i) for i in is ]
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
    j = combinedind(combiner(reverse([ mix(i) for i in is ])...; tags = ""))
    return ITensor(matrix(a, site...), j', dag(j))
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
    super_tensor(side_matrix(a.arg, site...), [ site_index(s, charged) for s in site ], true)

tensor(a::Right, site::AbstractSite...; charged::Bool = false) =
    super_tensor(side_matrix(a.arg, site...), [ site_index(s, charged) for s in site ], false)

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
    j = mix(i)
    v = state(site, a.state)
    if v isa Matrix
        m = v
    else
        m = v * v'
    end
    # the target is laid on the two site indices and only then gathered, never written
    # straight onto the mixed one: combining charged indices merges and sorts their
    # sectors, so the flat order of the mixed basis is not the order of the matrix
    tr = denseblocks(delta(dag(i), i')) * dag(mixer(i, j))
    t = charged_state(() -> op_on_sites(m, [i'], [dag(i'')]), i, a.state, site)
    return tr * t * mixer(i', j')
end

