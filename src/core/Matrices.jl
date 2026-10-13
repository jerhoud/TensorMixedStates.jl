# The matrix of an operator, a channel or a generator on sites and the vector or density matrix
# of a local state, which call each other: a projector on a named state needs its vector, and an
# eigenstate `A => λ` the matrix of `A`. The later layers also use `expand_sites`,
# `eigenspaces`, `nearly`, `rounding_tol` and `spectral_tol`.

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

the distance below which two eigenvalues among `values` are one, relative to the largest
"""
spectral_tol(values) = rounding_tol * maximum(abs, values)

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
    graded_parts(m, f)

the parts of each fermionic parity of the matrix `m` of an operator of one site, `f` the parity
of the site: pairs `(p, part)`, a single one for an operator of definite parity
"""
function graded_parts(m::AbstractMatrix, f::AbstractMatrix)
    p = parity_of(m, f)
    if !isnothing(p)
        return [ (p, m) ]
    end
    return [ (0, (m + f * m * f) / 2), (1, (m - f * m * f) / 2) ]
end

"""
    dense_max

the largest dimension of the dense matrix of a function of an operator, or of a superoperator,
built on sites: beyond, it is refused rather than exhausting memory
"""
const dense_max = 4096

"""
    check_dense(a, sites)

refuse the matrix of the function `a` on `sites` when its dimension, that of the states of the
sites or, for a channel or a generator, of their density matrices, is beyond `dense_max`
"""
function check_dense(a::GenericOp{K}, sites) where K
    ss = expand_sites(a, sites)
    d = prod(dim, ss)^(K === IsOperator ? 1 : 2)
    if d > dense_max
        error("$a on $(length(ss)) sites needs a dense matrix of dimension $d, beyond the " *
              "$dense_max TMS builds")
    end
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
    mixed_order(dims)

the positions in vec(ρ), ρ a density matrix of sites of dimensions `dims`, of its elements in
the order of the matrices of channels and generators, see `matrix`
"""
function mixed_order(dims)
    n = length(dims)
    ds = reverse(dims)
    # vec(ρ) as the kets, then the bras, each from the last site to the first
    axes = reshape(1:prod(dims)^2, ds..., ds...)
    return vec(permutedims(axes, vec(permutedims([1:n (n + 1):2n]))))
end

"""
    sided(a, b, sites)

the matrix of the map ρ ↦ a ρ b of the density matrices of `sites`, see `matrix`
"""
function sided(a::AbstractMatrix, b::AbstractMatrix, sites)
    o = mixed_order(map(dim, sites))
    return kron(transpose(b), a)[o, o]
end

"""
    kraus_sum(m, sites)

the sum of ``K^\\dagger K`` over the Kraus operators ``K`` of the channel of matrix `m` on
`sites`: the ``M`` such that ``\\mathrm{tr}\\,\\Phi(\\rho) = \\mathrm{tr}(M\\rho)``
"""
function kraus_sum(m::AbstractMatrix, sites)
    d = prod(dim, sites)
    o = mixed_order(map(dim, sites))
    # the row of the trace after the channel, which holds the transpose of M
    r = transpose(m) * vec(Matrix{Float64}(I, d, d))[o]
    return Matrix(transpose(reshape(r[invperm(o)], d, d)))
end

"""
    definite_parts(m, strings)

the parts of the matrix `m`, each of definite parity in every grading of `strings`, whose sum
is `m`, see `graded_parts`
"""
definite_parts(m::AbstractMatrix, strings::Vector) =
    foldl((ps, f) -> [ y for x in ps for (_, y) in graded_parts(x, f) ], strings; init = [m])

"""
    strung_kron(ms, strings)

the matrix of a tensor product of factors of matrices `ms`: `(A₁ ⊗ … ⊗ Aₙ)` on consecutive
sites is `A₁(1)…Aₙ(n)`, so a factor on fermionic sites takes each of its `strings`, one per
grading of its sites, for each factor after it odd in that grading. A factor after the first of
no definite parity is the sum of its parts that have one, see `definite_parts`, and the
product the sum of the products of these parts
"""
function strung_kron(ms::Vector, strings::Vector)
    parts = [ k == 1 ? [m] : definite_parts(m, fs) for (k, (m, fs)) in enumerate(zip(ms, strings)) ]
    return sum(Iterators.product(parts...)) do c
        ps = [ [ parity_of(m, f) for f in fs ] for (m, fs) in zip(c, strings) ]
        strung = map(eachindex(c)) do k
            m = c[k]
            for (g, f) in enumerate(strings[k])
                if f ≠ I && isodd(sum(p[g] for p in ps[k+1:end]; init = 0))
                    m = m * f
                end
            end
            m
        end
        foldl(kron, strung)
    end
end

"""
    matrix(op, site, sites...)

the matrix of the generic operator `op` on the given sites, a single one standing for several
identical ones, the last site varying fastest. A tensor product of fermionic factors holds the
Jordan-Wigner strings between them: the matrix of `C ⊗ dag(C)` is that of `C(1) * dag(C)(2)`.
A channel or a generator acts on the density matrices ρ of the sites: its matrix is that of
the map of vec(ρ), the last site varying fastest and, on each site, the ket faster than the
bra.

# Examples

    matrix(X, Qubit())
    matrix(X ⊗ A, Qubit(), Boson(2))
    matrix(Gate(X), Qubit())
"""
matrix(a::Matrix, ::AbstractSite, ::AbstractSite...) = copy(a)

matrix(a::Function, site::AbstractSite, sites::AbstractSite...) =
    matrix(a(site, sites...), site, sites...)

function matrix(a::IdentityOp{IsOperator, Generic}, site::AbstractSite, sites::AbstractSite...)
    d = prod(dim, expand_sites(a, (site, sites...)))
    return Matrix{Float64}(I, d, d)
end

function matrix(a::IdentityOp{<:Union{IsChannel, IsGenerator}, Generic}, site::AbstractSite,
                sites::AbstractSite...)
    d = prod(dim, expand_sites(a, (site, sites...)))^2
    return Matrix{Float64}(I, d, d)
end

function matrix(a::Operator, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    def = a.expr
    if isnothing(def)
        def = operator_definition(only(sites), Val(Symbol(a.name)))
        if isnothing(def)
            error("operator $(a.name) is not defined for site $(typeof(only(sites)))")
        end
    end
    m = matrix(def, sites...)
    n = prod(dim, sites)
    if size(m) ≠ (n, n)
        error("operator $(a.name) is given by a $(size(m, 1))×$(size(m, 2)) matrix and cannot " *
              "act on $(join(sites, " ⊗ ")), whose dimension is $n")
    end
    return m
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

function matrix(a::ExpOp, site::AbstractSite, sites::AbstractSite...)
    check_dense(a, (site, sites...))
    return exp(matrix(a.arg, site, sites...))
end

function matrix(a::ModOp, site::AbstractSite, sites::AbstractSite...)
    check_dense(a, (site, sites...))
    return exp(2im * π * matrix(a.arg, site, sites...) / a.modulus)
end

function matrix(a::IntPowOp, site::AbstractSite, sites::AbstractSite...)
    check_dense(a, (site, sites...))
    return matrix(a.arg, site, sites...) ^ a.expo
end

function matrix(a::GenPowOp, site::AbstractSite, sites::AbstractSite...)
    check_dense(a, (site, sites...))
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

# the string of a factor is the parity of its sites, see `strung_kron`
function matrix(a::TensorOp{IsOperator}, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    pos = factor_sites(a)
    ms = [ matrix(o, sites[p]...) for (o, p) in zip(a.subs, pos) ]
    return strung_kron(ms, [ [ foldl(kron, [ matrix(F, s) for s in sites[p] ]) ] for p in pos ])
end

# the strings of a factor are the parity of its sites on the kets and on the bras
function matrix(a::TensorOp{IsChannel}, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    pos = factor_sites(a)
    ms = [ matrix(o, sites[p]...) for (o, p) in zip(a.subs, pos) ]
    strings = map(pos) do p
        f = foldl(kron, [ matrix(F, s) for s in sites[p] ])
        [ sided(f, one(f), sites[p]), sided(one(f), f, sites[p]) ]
    end
    return strung_kron(ms, strings)
end

function matrix(a::Gate{Generic}, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    k = matrix(a.arg, sites...)
    return sided(k, k', sites)
end

function matrix(a::Reset, site::AbstractSite, others::AbstractSite...)
    s = only(expand_sites(a, (site, others...)))
    v = local_state(s, a.state)
    σ = v isa Vector ? v * v' : v
    return vec(σ) * transpose(vec(Matrix{Float64}(I, dim(s), dim(s))))
end

function matrix(a::Dephase, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    check_dense(a, sites)
    _, ps = eigenspaces(matrix(a.arg, sites...), "Dephase", a.arg, " on " * join(sites, " ⊗ "))
    return sum(sided(p, p, sites) for p in ps)
end

# the dissipator of a channel is the sum of those of its Kraus operators K, which take the sum
# of the K†K
function matrix(a::Dissipator{Generic}, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    if a.arg isa GenericOp{IsOperator}
        l = matrix(a.arg, sites...)
        g, m = sided(l, l', sites), l' * l
    else
        g = matrix(a.arg, sites...)
        m = kraus_sum(g, sites)
    end
    return g - (sided(m, one(m), sites) + sided(one(m), m, sites)) / 2
end

function matrix(a::Evolver{Generic}, site::AbstractSite, others::AbstractSite...)
    sites = expand_sites(a, (site, others...))
    k = matrix(a.arg, sites...)
    return sided(k, one(k), sites) + sided(one(k), k', sites)
end

"""
    eigenprojector(a, λ, sites...)

the projector on the eigenspace of eigenvalue `λ` of the Hermitian operator `a` on `sites`
"""
function eigenprojector(a::GenericOp{IsOperator}, λ::Real, sites::AbstractSite...)
    check_dense(a, sites)
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
