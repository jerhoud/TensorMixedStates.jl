export tensor, matrix

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
tensor(a::GenericOp, site::AbstractSite...) =
    tensor(matrix(a, site...), site...)

function tensor(a::Matrix, ::AbstractSite, ::AbstractSite...)
    n, _ = size(a)
    i = Index(n)
    ITensor(a, i', i)   
end

matrix(a::Matrix, ::AbstractSite, ::AbstractSite...) = a

# a multi site operator may be defined by a function of its sites, which lets it adapt to
# them, as in a gate whose dimension is read from the site it acts on
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

function tensor(a::Left, site::AbstractSite...)
    ti = tensor(a.arg, site...)
    i = tensor_index(ti)
    is = Index.(site)
    j = sim(i)
    js = sim.(is)
    ci = combinerto(i, reverse(is)...)
    cj = combinerto(j, reverse(js)...)
    ijs = Iterators.flatten(zip(reverse(is), reverse(js)))
    c = combiner(ijs...; tags="")
    return (ti * ci * ci') * (delta(j, j') * cj * cj') * c * c'
end

function tensor(a::Right, site::AbstractSite...)
    ti = tensor(a.arg, site...)
    i = tensor_index(ti)
    is = Index.(site)
    j = sim(i)
    js = sim.(is)
    ci = combinerto(i, reverse(is)...)
    cj = combinerto(j, reverse(js)...)
    jis = Iterators.flatten(zip(reverse(js), reverse(is)))
    c = combiner(jis...; tags="")
    return (delta(j, j') * cj * cj') * (dag(ti) * ci * ci') * c * c'
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

function tensor(a::SetState, site::AbstractSite)
    i = Index(site)
    j = mix(i)
    c = combinerto(j, i, i')
    v = state(site, a.state)
    if v isa Matrix
        m = v
    else
        m = v * v'
    end
    return dense(delta(i, i')) * c * ITensor(m, j')
end

