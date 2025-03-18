export tensor, matrix, fermionic, has_fermionic, has_multiop

function combinerto(i::Index, j::Index...)
    c = combiner(j...)
    x = combinedind(c)
    replaceind(c, x, i)
end

tensor_index(t::ITensor) = getfirst(i->hasplev(i, 0), inds(t))

function matrix(a, site::AbstractSite...)
    t = tensor(a, site...)
    i = tensor_index(t)
    Matrix(t, i', i)
end

function tensor(a, site::AbstractSite...)
    m = matrix(a, site...)
    n, _ = size(m)
    i = Index(n)
    ITensor(m, i', i)
end

matrix(a::Matrix, ::AbstractSite) = a

matrix(a::Function, site::AbstractSite) = a(site)

matrix(a::String, site::AbstractSite) =
    matrix(operator_info(site, a), site)

matrix(a::Operator, site::AbstractSite...) =
    if isnothing(a.expr)
        matrix(a.name, site...)
    else
        matrix(a.expr, site...)
    end

matrix(a::ProdOp, site::AbstractSite...) =
    a.coef * prod(matrix(o, site...) for o in a.subs)

matrix(a::SumOp, site::AbstractSite...) =
    sum(matrix(o, site...) for o in a.subs)

matrix(a::ExpOp, site::AbstractSite...) =
    exp(matrix(a.arg, site...))

matrix(a::SqrtOp, site::AbstractSite...) =
    sqrt(matrix(a.arg, site...))

matrix(a::PowOp, site::AbstractSite...) =
    matrix(a.arg, site...) ^ a.expo

matrix(a::DagOp, site::AbstractSite...) =
    adjoint(matrix(a.arg, site...))

function tensor(a::Gate, site::AbstractSite...)
    ti = tensor(a.arg, site...)
    i = tensor_index(ti)
    j = sim(i)
    tj = replaceind(ti, (i, i'), (j, j'))
    c = combiner(j, i)
    return ti * dag(tj) * c * c'
end

function tensor(a::Left, site::AbstractSite...)
    ti = tensor(a.arg, site...)
    i = tensor_index(ti)
    j = sim(i)
    c = combiner(j, i)
    return ti * delta(j, j') * c * c'
end

function tensor(a::Right, site::AbstractSite...)
    ti = tensor(a.arg, site...)
    i = tensor_index(ti)
    j = sim(i)
    c = combiner(i, j)
    return dag(ti) * delta(j, j') * c * c'
end

tensor_next(o::ExprOp{N}, site::Vararg{AbstractSite, M}) where {N, M} =
    (tensor(o, site[1:N]...), site[N+1:M])

function tensor(a::TensorOp{N}, site::AbstractSite...) where N
    if length(site) == 1
        ts = [tensor(o, site...) for o in a.subs]
    elseif length(site) ≠ N
        error("number of sites does not match operator")
    else
        rest = site
        ts = map(a.subs) do o
            t, rest = tensor_next(o, rest...)
            t
        end 
    end
    c = combiner((tensor_index(t) for t in reverse(ts))...)
    c * prod(ts) * c'
end


function fermionic(a::SumOp{1})
    n = length(a.subs)
    nf = count(fermionic, a.subs)
    if nf == n
        return true
    elseif n == 0
        return false
    else
        error("cannot sum fermionic with non fermionic operators")
    end
end

fermionic(a::Operator{1}) = a.fermionic
fermionic(a::ProdOp{1}) = isodd(count(fermionic, a.subs))
fermionic(a::Union{Gate{1}, Dissipator{1}, Indexed{1}}) = fermionic(a.arg)

has_fermionic(a::Union{ProdOp, SumOp, TensorOp}) = any(has_fermionic, a.subs)
has_fermionic(a::Union{ExpOp, SqrtOp, PowOp, Gate, Dissipator}) = has_fermionic(a.arg)
has_fermionic(a::Indexed) = has_fermionic(a.op)
has_fermionic(a::Operator) = a.fermionic

has_multiop(a::Union{ProdOp, SumOp, TensorOp}) = any(has_multiop, a.subs)
has_multiop(a::Union{ExpOp, SqrtOp, PowOp}) = has_multiop(a.arg)
has_multiop(a::Indexed) = has_multiop(a.op)
has_multiop(a::Union{Operator{1}, Dissipator{1}, Gate{1}}) = false
has_multiop(a::Union{Operator, Dissipator, Gate}) = true

function signature(p)
    n = length(p)
    t = fill(false, n)
    r = 1
    for i in 1:n
        if t[i] continue end
        j = p[i]
        while j ≠ i
            r = -r
            t[j] = true
            j = p[j]
        end
    end
    return r
end
