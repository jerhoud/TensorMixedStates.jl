export tensor, matrix, fermionic, left_tensor, right_tensor

function combinerto(i::Index, j::Index...)
    c = combiner(j...; tags="")
    x = combinedind(c)
    replaceind(c, x, i)
end

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
function tensor(a::GenericOp, site::AbstractSite...)
    m = matrix(a, site...)
    n, _ = size(m)
    i = Index(n)
    ITensor(m, i', i)
end

matrix(a::Matrix, ::AbstractSite, ::AbstractSite...) = a

matrix(a::Function, site::AbstractSite) =
    matrix(a(site), site)

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
    return st * adjoint(st)
end

matrix(a::ProdGenOp, site::AbstractSite...) =
    a.coef * prod(matrix(o, site...) for o in a.subs)

matrix(a::SumGenOp, site::AbstractSite...) =
    sum(matrix(o, site...) for o in a.subs)

matrix(a::ExpOp, site::AbstractSite...) =
    exp(matrix(a.arg, site...))

matrix(a::PowOp, site::AbstractSite...) =
    matrix(a.arg, site...) ^ a.expo

matrix(a::DagOp, site::AbstractSite...) =
    adjoint(matrix(a.arg, site...))

matrix(::Dissipator, ::AbstractSite...) = error("cannot give matrix nor tensor for Dissipator")

function tensor(a::Gate, site::AbstractSite...)
    ti = tensor(a.arg, site...)
    i = tensor_index(ti)
    is = map(Index, site)
    j = sim(i)
    js = sim.(is)
    tj = replaceinds(ti, (i, i'), (j, j'))
    ci = combinerto(i, reverse(is)...)
    cj = combinerto(j, reverse(js)...)
    ijs = Iterators.flatten(zip(is, js))
    c = combiner(ijs...; tags="")
    return (ti * ci * ci') * (dag(tj) * cj * cj') * c * c'
end

function left_tensor(a::GenericOp, site::AbstractSite...)
    ti = tensor(a, site...)
    i = tensor_index(ti)
    j = sim(i)
    c = combiner(i, j; tags="")
    return ti * delta(j, j') * c * c'
end

function right_tensor(a::GenericOp, site::AbstractSite...)
    ti = tensor(a, site...)
    i = tensor_index(ti)
    j = sim(i)
    c = combiner(j, i; tags="")
    return dag(ti) * delta(j, j') * c * c'
end

tensor_next(f, o::GenericOp{Pure, N}, site::Vararg{AbstractSite, M}) where {N, M} =
    (f(o, site[1:N]...), site[N+1:M])

function tensor_apply(f, a::TensorOp{N}, idx::Vararg{AbstractSite, N}) where N
    rest = idx
    r = map(a.subs) do o 
        t, rest = tensor_next(f, o, rest...)
        t
    end
end

tensor_next(f, o::GenericOp{Pure, N}, site::Vararg{Int, M}) where {N, M} =
    (f(o, site[1:N]...), site[N+1:M])

function tensor_apply(f, a::TensorOp{N}, idx::Vararg{Int, N}) where N
    rest = idx
    r = map(a.subs) do o 
        t, rest = tensor_next(f, o, rest...)
        t
    end
end

function tensor(a::TensorOp{N}, site::AbstractSite...) where N
    if length(site) == 1
        ts = [tensor(o, site...) for o in a.subs]
    elseif length(site) ≠ N
        error("number of sites does not match operator")
    else
        ts = tensor_apply(tensor, a, site...)
    end
    c = combiner((tensor_index(t) for t in reverse(ts))...; tags="")
    c * prod(ts) * c'
end

