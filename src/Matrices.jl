# The matrix of an operator on sites and the vector or density matrix of a local state, which
# call each other: a projector on a named state needs its vector, and an eigenstate `A => λ`
# the matrix of `A`. The later layers also use `expand_sites`, `check_size`, `eigenspaces`,
# `nearly`, `rounding_tol` and `spectral_tol`.

export matrix

"""
    rounding_tol

the size, relative to the norm of a matrix, below which an element, a singular value or a
difference is rounding: a tolerance, not a setting
"""
const rounding_tol = 1e-13

"""
    nearly(a, b)

whether `a` and `b` are equal up to rounding, relative to the larger of the two
"""
nearly(a, b) = norm(a - b) ≤ rounding_tol * max(norm(a), norm(b))

"""
    spectral_tol(values)

the distance below which two eigenvalues among `values` are one
"""
spectral_tol(values) = rounding_tol * max(1, maximum(abs, values))

"""
    expand_sites(what, n, sites)
    expand_sites(a::GenericOp, sites)

the `n` sites `what` acts on, a single site standing for `n` identical ones
"""
function expand_sites(what, n::Int, sites)
    ss = length(sites) == 1 ? fill(only(sites), n) : collect(AbstractSite, sites)
    if length(ss) ≠ n
        error("$what acts on $n site$(n == 1 ? "" : "s") and was given $(length(sites))")
    end
    return ss
end

expand_sites(a::GenericOp, sites) = expand_sites(a, nsites(a), sites)

"""
    check_size(what, m, sites)

the matrix `m` of `what`, refused when its size is not the dimension of `sites`
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
    parity_of(m, f)

the fermionic parity of the matrix `m`, `f` being the parity of its sites: 0 if they commute,
1 if they anticommute, `nothing` otherwise
"""
parity_of(m::AbstractMatrix, f::AbstractMatrix) =
    if nearly(f * m * f, m)
        0
    elseif nearly(f * m * f, -m)
        1
    else
        nothing
    end

"""
    matrix(op, site, sites...)

the matrix of the generic operator `op` on the given sites, a single one standing for several
identical ones, the last site varying fastest. A tensor product of fermionic factors holds the
Jordan-Wigner strings between them: the matrix of `C ⊗ dag(C)` is that of `C(1) * dag(C)(2)`.

# Examples

    matrix(X, Qubit())
    matrix(X ⊗ A, Qubit(), Boson(2))
"""
matrix(a::Matrix, ::AbstractSite, ::AbstractSite...) = copy(a)

matrix(a::Function, site::AbstractSite, sites::AbstractSite...) =
    matrix(a(site, sites...), site, sites...)

function matrix(a::IdentityOp{IsOperator, Generic}, site::AbstractSite, sites::AbstractSite...)
    d = prod(dim, expand_sites(a, (site, sites...)))
    return Matrix{Float64}(I, d, d)
end

function matrix(a::Operator, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    if isnothing(a.expr)
        def = operator_definition(only(sites), Val(Symbol(a.name)))
        if isnothing(def)
            error("operator $(a.name) is not defined for site $(typeof(only(sites)))")
        end
        return check_size("operator $(a.name)", matrix(def, sites...), sites)
    end
    return check_size("operator $(a.name)", matrix(a.expr, sites...), sites)
end

function matrix(a::Proj, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    if a.state isa Pair
        return eigenprojector(a.state.first, a.state.second, sites...)
    end
    v = local_state(only(sites), a.state)
    if !(v isa Vector)
        error("Proj projects on a pure state, and $(repr(a.state)) is mixed on $(only(sites))")
    end
    return v * v'
end

matrix(a::ScalarOp, site::AbstractSite, sites::AbstractSite...) =
    a.coef * matrix(a.arg, site, sites...)

matrix(a::ProdOp, site::AbstractSite, sites::AbstractSite...) =
    prod(matrix(o, site, sites...) for o in a.subs)

matrix(a::SumOp, site::AbstractSite, sites::AbstractSite...) =
    sum(matrix(o, site, sites...) for o in a.subs)

matrix(a::ExpOp, site::AbstractSite, sites::AbstractSite...) = exp(matrix(a.arg, site, sites...))

matrix(a::ModOp, site::AbstractSite, sites::AbstractSite...) =
    exp(2im * π * matrix(a.arg, site, sites...) / a.modulus)

matrix(a::IntPowOp, site::AbstractSite, sites::AbstractSite...) =
    matrix(a.arg, site, sites...) ^ a.expo

function matrix(a::GenPowOp, site::AbstractSite, sites::AbstractSite...)
    m = matrix(a.arg, site, sites...)
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

matrix(a::DagOp, site::AbstractSite, sites::AbstractSite...) =
    collect(adjoint(matrix(a.arg, site, sites...)))

# (A₁ ⊗ … ⊗ Aₙ) on consecutive sites is A₁(1)…Aₙ(n): a factor on fermionic sites takes their
# parity for each odd factor after it, and a factor of no definite parity there is refused
function matrix(a::TensorOp{IsOperator}, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    pos = factor_sites(a)
    fs = [ foldl(kron, [ matrix(F, s) for s in sites[p] ]) for p in pos ]
    ms = [ matrix(o, sites[p]...) for (o, p) in zip(a.subs, pos) ]
    ps = [ parity_of(m, f) for (m, f) in zip(ms, fs) ]
    strung = map(eachindex(ms)) do k
        rest = ps[k+1:end]
        j = findfirst(isnothing, rest)
        if fs[k] == I
            ms[k]
        elseif !isnothing(j)
            error("$a has a factor of no definite fermionic parity, $(a.subs[k+j]), after a " *
                  "fermionic site: write it as a sum of tensor products")
        elseif isodd(sum(rest; init = 0))
            ms[k] * fs[k]
        else
            ms[k]
        end
    end
    return foldl(kron, strung)
end

"""
    eigenspaces(m, what, a, on = "")

the eigenvalues of the Hermitian matrix `m` of `a`, increasing, and the projectors on their
eigenspaces, eigenvalues equal to rounding gathered into one, `what` naming the caller in the
refusal of a matrix that is not Hermitian
"""
function eigenspaces(m::AbstractMatrix, what::String, a, on::String = "")
    if !nearly(m, m')
        error("$what needs a Hermitian operator, which $a is not$on")
    end
    # real projectors for a real operator, which give real tensors, cheaper to contract
    e = eltype(m) <: Real ? eigen(Symmetric(m)) : eigen(Hermitian(m))
    tol = spectral_tol(e.values)
    d = size(m, 1)
    values = Float64[]
    ps = Matrix{eltype(e.vectors)}[]
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

"""
    eigenprojector(a, λ, sites...)

the projector on the eigenspace of eigenvalue `λ` of the Hermitian operator `a` on `sites`
"""
function eigenprojector(a::GenericOp{IsOperator}, λ::Real, sites::AbstractSite...)
    on = " on " * join(sites, " ⊗ ")
    values, ps = eigenspaces(matrix(a, sites...), "Proj", a, on)
    k = findfirst(v -> abs(v - λ) ≤ spectral_tol(values), values)
    if isnothing(k)
        error("$λ is not an eigenvalue of $a$on, whose eigenvalues are $(join(values, ", "))")
    end
    return ps[k]
end

"""
    local_state(site, state)

the vector, or the density matrix, of the local state `state` on `site`: a name, `A => λ`, a
vector of amplitudes, a density matrix, or a function of the site returning one of those. A
name is looked for among the states the site declares, then in its rule, see
`string_state`, then among the states every site has.
"""
function local_state(site::AbstractSite, name::String)
    declared = state_definition(site, Val(Symbol(name)))
    if !isnothing(declared)
        return local_state(site, declared)
    end
    ruled = string_state(site, name)
    if !isnothing(ruled)
        return local_state(site, ruled)
    end
    common = get(common_states, name, nothing)
    if !isnothing(common)
        return local_state(site, common(site))
    end
    error("state $name is not defined for site $(typeof(site))")
end

# an eigenvector of a diagonal operator for a simple eigenvalue is a basis vector, which a pure
# state takes; any other eigenstate is the uniform mixture of its eigenspace
function local_state(site::AbstractSite, p::Pair{<:SimpleOp, <:Real})
    m = matrix(first(p), site)
    pr = eigenprojector(first(p), last(p), site)
    k = findfirst(x -> abs(x - 1) ≤ rounding_tol, real(diag(pr)))
    if isdiag(m) && !isnothing(k) && nearly(pr, pr[:, k:k] * pr[k:k, :])
        v = fill(0., dim(site))
        v[k] = 1.
        return v
    end
    return pr / tr(pr)
end

function local_state(site::AbstractSite, v::Vector{<:Number})
    if length(v) ≠ dim(site)
        error("a state of $(length(v)) components cannot be one of $site, whose dimension is $(dim(site))")
    end
    return copy(v)
end

function local_state(site::AbstractSite, m::Matrix{<:Number})
    d = dim(site)
    if size(m) ≠ (d, d)
        error("a $(size(m, 1))×$(size(m, 2)) density matrix cannot be one of $site, whose " *
              "dimension is $d")
    end
    return copy(m)
end

local_state(site::AbstractSite, f::Function) = local_state(site, f(site))

local_state(::AbstractSite, n::Number) =
    error("a number does not denote a state: write \"$n\" or Level => $n")

local_state(site::AbstractSite, x) =
    error("$(repr(x)) is not a local state of $site: give a name, A => λ, a vector, a density " *
          "matrix or a function of the site")
