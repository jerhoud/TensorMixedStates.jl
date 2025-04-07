export tensor, matrix, fermionic

function combinerto(i::Index, j::Index...)
    c = combiner(j...; tags="")
    x = combinedind(c)
    replaceind(c, x, i)
end

tensor_index(t::ITensor) = getfirst(i->hasplev(i, 0), inds(t))

function matrix(a::ExprOp, site::AbstractSite...)
    t = tensor(a, site...)
    i = tensor_index(t)
    Matrix(t, i', i)
end

function tensor(a::ExprOp, site::AbstractSite...)
    m = matrix(a, site...)
    n, _ = size(m)
    i = Index(n)
    ITensor(m, i', i)
end

matrix(::ExprIndexed, ::AbstractSite...) = error("cannot give matrix nor tensor for indexed expressions")

matrix(a::Matrix, ::AbstractSite, ::AbstractSite...) = a

matrix(a::Function, site::AbstractSite, ::AbstractSite...) = a(site)

matrix(a::String, site::AbstractSite, ::AbstractSite...) =
    matrix(operator_info(site, a), site)

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

matrix(::Union{Dissipator, Evolver}, ::AbstractSite...) = error("cannot give matrix nor tensor for Dissipator or Evolver")

function tensor(a::Gate, site::AbstractSite...)
    ti = tensor(a.arg, site...)
    i = tensor_index(ti)
    j = sim(i)
    tj = replaceinds(ti, (i, i'), (j, j'))
    c = combiner(j, i; tags="")
    return ti * dag(tj) * c * c'
end

function tensor(a::Left, site::AbstractSite...)
    ti = tensor(a.arg, site...)
    i = tensor_index(ti)
    j = sim(i)
    c = combiner(j, i; tags="")
    return ti * delta(j, j') * c * c'
end

function tensor(a::Right, site::AbstractSite...)
    ti = tensor(a.arg, site...)
    i = tensor_index(ti)
    j = sim(i)
    c = combiner(i, j; tags="")
    return dag(ti) * delta(j, j') * c * c'
end

tensor_next(f, o::ExprOp{T, N}, site::Vararg{U, M}) where {T, N, U, M} =
    (f(o, site[1:N]...), site[N+1:M])

function tensor_apply(f, a::TensorOp{T, N}, idx::Vararg{U, N}) where {T, U, N}
    rest = idx
    r = map(a.subs) do o 
        t, rest = tensor_next(f, o, rest...)
        t
    end
end

function tensor(a::TensorOp{T, N}, site::AbstractSite...) where {T, N}
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


function fermionic(a::SumOp)
    n = length(a.subs)
    nf = count(fermionic, a.subs)
    if nf == n
        return true
    elseif nf == 0
        return false
    else
        error("cannot sum fermionic with non fermionic operators: $a")
    end
end

fermionic(a::Operator) = a.fermionic
fermionic(a::ProdOp) = isodd(count(fermionic, a.subs))
fermionic(a::Indexed) = fermionic(a.op)
fermionic(a::Union{Gate, Left, Right, DagOp}) = fermionic(a.arg)
fermionic(a::Union{ExpOp, SqrtOp, PowOp}) =
    if fermionic(a.arg)
        error("cannot take functionals of fermionic operators: $a")
    else
        false
    end
fermionic(::Proj) = false
fermionic(a) = error("bug: fermionic($a)")

has_dissipator(a) = eval_expr(has_dissipator, a)
function has_dissipator(a::Union{ProdOp, TensorOp})
    d = any(has_dissipator, a.subs)
    if d && length(a.subs) >= 2
        error("cannot multiply dissipators with other operators ($a)")
    end
    return d
end
has_dissipator(a::Dissipator) = true
has_dissipator(a::Operator) =
    if a.expr isa ExprOp
        has_dissipator(a.expr)
    else
        false
    end
has_dissipator(::Proj) = false
has_dissipator(a::Union{ExpOp, SqrtOp, PowOp, DagOp}) =
    if has_dissipator(a.arg)
        error("cannot take functional of dissipators ($a)")
    else
        false
    end