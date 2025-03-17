export tensor, matrix, fermionic

function matrix(a, site::AbstractSite...)
    t = tensor(a, site...)
    i = getfirst(j->hasplev(j, 0), inds(t))
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
    matrix(first(operator_info(site, a)), site)

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
    c = combiner((getfirst(j->hasplev(j, 0), inds(t)) for t in reverse(ts))...)
    c * prod(ts) * c'
end

fermionic(a::Operator{1}) = a.fermionic

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

fermionic(a::ProdOp{1}) = isodd(count(fermionic, a.subs))

fermionic(a::Gate{1}) = fermionic(a.arg)

fermionic(a::Dissipator{1}) = fermionic(a.arg)

fermionic(a::Indexed{1}) = fermionic(a.op)