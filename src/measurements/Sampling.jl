# The projective measurements of a state: the results of measuring a site or a Hermitian
# operator with their probabilities, a result drawn, and the state projected on a result.

export probabilities, sample, collapse

"""
    measured_spectrum(system, op, what)

the sites of `system` the Hermitian operator `op` acts on, increasing, its eigenvalues and the
projectors on its eigenspaces, matrices on those sites in the Jordan-Wigner basis of those sites
alone, see `eigenspaces`; `what` names the caller in the refusals
"""
function measured_spectrum(system::System, op::IndexedOp{IsOperator}, what::String)
    sites = sort(sites_of(op))
    if isempty(sites)
        error("$what needs an operator placed on sites, which $op is not")
    end
    # a Float64, which does not wrap around as an Int does on many sites
    d = prod(float(dim(system[k])) for k in sites)
    if d > dense_max
        error("$what of $op needs its matrix on its $(length(sites)) sites, of dimension " *
              "$(Printf.@sprintf("%.0f", d)), beyond the $dense_max TMS builds")
    end
    on = " on " * join([ system[k] for k in sites ], " ⊗ ")
    values, ps = eigenspaces(matrix_on(system, op, sites), what, op, on)
    return sites, values, ps
end

"""
    even_on(system, i, m)

whether the matrix `m` of an operator on site `i` of `system` is even, so that it commutes with
the operators of other sites
"""
even_on(system::System, i::Int, m::AbstractMatrix) = parity_of(m, matrix(F, system[i])) == 0

"""
    involution(system, op)

`(c, P)` when `op` is `c P`, `c` real and not zero and `P` a product of Hermitian involutions of
one site on distinct sites, two at least, each even, so that they commute, and diagonal on a
system conserving something; `nothing` otherwise
"""
function involution(system::System, op::IndexedOp{IsOperator})
    c = scalarcoef(op)
    fs = IndexedOp{IsOperator}[]
    for x in prodsubs(scalararg(op))
        c *= scalarcoef(x)
        a = scalararg(x)
        if !(a isa AtIndex && length(a.index) == 1)
            return nothing
        end
        i = only(a.index)
        m = matrix(a.op, system[i])
        if !(nearly(m, m') && nearly(m * m, one(m)) && even_on(system, i, m)) ||
           (is_charged(system) && !nearly(m, Diagonal(diag(m))))
            return nothing
        end
        push!(fs, a)
    end
    sites = [ only(a.index) for a in fs ]
    if length(sites) < 2 || !allunique(sites) || !isreal(c) || iszero(c)
        return nothing
    end
    return (float(real(c)), ProdOp(fs))
end

"""
    counting(system, op, what)

the sites of `op` and, for each, twice the eigenvalues of its part on that site, as integers,
with the projectors on their eigenspaces, when `op` is a sum of operators of one site on two
sites at least, each even, so that they commute, of half integer eigenvalues and diagonal on a
system conserving something; `nothing` otherwise. `what` names the caller in the refusals.
"""
function counting(system::System, op::IndexedOp{IsOperator}, what::String)
    c = scalarcoef(op)
    parts = Dict{Int, Vector{Any}}()
    for x in sumsubs(scalararg(op))
        a = scalararg(x)
        if !(a isa AtIndex && length(a.index) == 1)
            return nothing
        end
        push!(get!(parts, only(a.index), []), c * scalarcoef(x) * a.op)
    end
    sites = sort(collect(keys(parts)))
    if length(sites) < 2
        return nothing
    end
    locals = map(sites) do i
        m = matrix(sum(parts[i]), system[i])
        if !even_on(system, i, m) || (is_charged(system) && !nearly(m, Diagonal(diag(m))))
            return nothing
        end
        values, ps = eigenspaces(m, what, op, " on $(system[i])")
        if any(λ -> abs(2λ - round(2λ)) > 2 * spectral_tol(values), values)
            return nothing
        end
        return (round.(Int, 2 .* values), ps)
    end
    return any(isnothing, locals) ? nothing : (sites, locals)
end

"""
    counting_probabilities(state, sites, locals)

the totals, twice the values of the sum that `counting` describes, with their probabilities, in
increasing order: the totals reached lie within `K` successive integers, whose probabilities are
read off the mean values of the `K` products over the sites of ``e^{2i\\pi k t/K}``, `t` twice
the value on the site
"""
function counting_probabilities(state::State, sites, locals)
    totals = foldl((s, l) -> Set(x + t for x in s for t in first(l)), locals; init = Set(0))
    t0, t1 = extrema(totals)
    K = t1 - t0 + 1
    φ = map(0:K-1) do k
        phase(t, p) = Operator{1}("Count$k", sum(cis(2π * k * u / K) * q for (u, q) in zip(t, p)))
        expect(state, prod(phase(l...)(i) for (i, l) in zip(sites, locals)))
    end
    return [ T => real(sum(φ[k+1] * cis(-2π * k * T / K) for k in 0:K-1)) / K
             for T in sort(collect(totals)) ]
end

"""
    count_mpo(system, sites, locals, T)

the MPO of the projector on the total `T` of the sum that `counting` describes, on the pure
indices of `system`: its bond carries the totals of the sites on its left that can still reach
`T`
"""
function count_mpo(system::System, sites, locals, T::Int)
    n = length(system)
    at = Dict(i => k for (k, i) in enumerate(sites))
    steps(i) = haskey(at, i) ? first(locals[at[i]]) : [0]
    reached = [ Set(0) ]
    for i in 1:n
        push!(reached, Set(s + t for s in reached[i] for t in steps(i)))
    end
    bonds = Vector{Vector{Int}}(undef, n + 1)
    bonds[n + 1] = [T]
    for i in n:-1:1
        bonds[i] = sort([ s for s in reached[i] if any(t -> (s + t) in bonds[i + 1], steps(i)) ])
    end
    links = [ is_charged(system) ? Index(QN() => length(b); tags = "Count,l=$(i - 1)") :
                                   Index(length(b); tags = "Count,l=$(i - 1)")
              for (i, b) in enumerate(bonds) ]
    # real for real projectors, which keep a real state real
    T = promote_type(Float64, (eltype(p) for (_, ps) in locals for p in ps)...)
    ts = map(1:n) do i
        x = siteind(Pure, system, i)
        d = dim(x)
        right = Dict(s => q for (q, s) in enumerate(bonds[i + 1]))
        blocks = haskey(at, i) ? collect(zip(locals[at[i]]...)) : [ (0, Matrix{Float64}(I, d, d)) ]
        a = zeros(T, d, d, length(bonds[i]), length(bonds[i + 1]))
        for (p, s) in enumerate(bonds[i]), (t, m) in blocks
            q = get(right, s + t, nothing)
            if !isnothing(q)
                a[:, :, p, q] .+= m
            end
        end
        return charged_itensor(a, [x', dag(x), dag(links[i]), links[i + 1]])
    end
    ts[1] *= onehot(links[1] => 1)
    ts[n] *= onehot(dag(links[n + 1]) => 1)
    return MPO(ts)
end

"""
    project_parity(state, P, s; limits)

the state projected by ``(1 + sP)/2``, `P` an involution and `s` its eigenvalue ±1, and
normalized: a density matrix through the gate of that projector
"""
project_parity(state::State{Pure}, P, s; limits::Limits) =
    normalize(+(state, s * apply(P, state); limits))

project_parity(state::State{Mixed}, P, s; limits::Limits) =
    normalize(apply(Gate((Id(1) + s * P) / 2), state; limits))

"""
    project_count(state, sites, locals, T; limits)

the state projected on the total `T` of the sum that `counting` describes, and normalized: a
density matrix by the MPO of the projector on its kets, then on its bras, see `lifted`
"""
project_count(state::State{Pure}, sites, locals, T; limits::Limits) =
    normalize(apply(count_mpo(state.system, sites, locals, T), state; limits))

function project_count(state::State{Mixed}, sites, locals, T; limits::Limits)
    w = count_mpo(state.system, sites, locals, T)
    bra = relabeller(i -> star(i, strong_names(state.system)))
    for left in (true, false)
        state = apply(MPO([ lifted(state.system, w[k], k, left, bra) for k in eachindex(w) ]),
                      state; limits)
    end
    return normalize(state)
end

"""
    measurement(state, op, what)

the results of measuring the Hermitian operator `op`, in increasing order, their probabilities,
and the function `(k, limits) -> state` projecting the state on the `k`-th of them and
normalizing it: as a product of involutions, as a sum of operators of one site, or through the
projectors on the eigenspaces of its matrix, placed on its sites, whose Jordan-Wigner strings
the descent lays; `what` names the caller in the refusals
"""
function measurement(state::State, op::IndexedOp{IsOperator}, what::String)
    nonzero_trace(state, what)
    # a probability below rounding, of either sign, as a sum of Fourier leaves for an impossible
    # result, is zero
    zeroed(probs) = [ abs(p) ≤ rounding_tol ? zero(p) : p for p in probs ]
    system = state.system
    check_indices(system, op)
    inv = involution(system, op)
    if !isnothing(inv)
        c, P = inv
        m = real(expect(state, P))
        signs = c > 0 ? [-1, 1] : [1, -1]
        return (c .* signs, zeroed((1 .+ signs .* m) ./ 2),
                (k, limits) -> project_parity(state, P, signs[k]; limits))
    end
    cnt = counting(system, op, what)
    if !isnothing(cnt)
        ts = counting_probabilities(state, cnt...)
        return (first.(ts) ./ 2, zeroed(last.(ts)),
                (k, limits) -> project_count(state, cnt..., first(ts[k]); limits))
    end
    sites, values, ps = measured_spectrum(system, op, what)
    projectors = [ Operator{length(sites)}("Proj($op => $λ)", p)(sites...)
                   for (λ, p) in zip(values, ps) ]
    probs = [ real(expect(state, p)) for p in projectors ]
    function project(k, limits)
        if !definite_on(system, ps[k], sites)
            error("the result $(values[k]) of $op leaves the state without definite charges: " *
                  "collapse weaken(state, ()) instead")
        end
        return normalize(apply(projectors[k], state; limits))
    end
    return (values, zeroed(probs), project)
end

"""
    probabilities(state, pos)
    probabilities(state, op)

the results of measuring site `pos`, or the Hermitian operator `op` placed on a few sites, with
their probabilities, as pairs `result => probability` in increasing order, those of zero
probability included, a probability below rounding being zero: the levels of the site, counted from 0, or the eigenvalues of `op`, to
rounding. The matrix of `op` on its sites is diagonalized, its size being the product of their
dimensions, except for two kinds of operators measured on any number of sites: a product of
involutions, as a string of Pauli operators, by its mean value, and a sum of operators of one
site with half integer eigenvalues, as the number of particles of a region, by the mean values
of products of phases.

# Examples

    probabilities(state, 3)
    probabilities(state, Z(1) * Z(2))
    probabilities(state, sum(N(i) for i in 1:10))
"""
function probabilities(state::State, op::IndexedOp{IsOperator})
    values, probs, _ = measurement(state, op, "probabilities")
    return [ λ => p for (λ, p) in zip(values, probs) ]
end

function probabilities(state::State, pos::Int)
    check_positions(state, [pos], "probabilities")
    return [ round(Int, x) => p for (x, p) in probabilities(state, Level(pos)) ]
end

"""
    draw(p, d, rnd)

the first outcome among `0:d-1` whose cumulated probability exceeds `rnd`, `p(x)` being the
probability of `x`, the last one taking what the others leave
"""
function draw(p, d::Int, rnd)
    ptot = 0.
    for x in 0:d - 2
        ptot += p(x)
        if rnd < ptot
            return x
        end
    end
    return d - 1
end

"""
    level_projector(d, x)

the matrix of the projector on the level `x`, counted from 0, of a site of dimension `d`
"""
level_projector(d::Int, x::Int) = Matrix{Float64}(Diagonal([ k == x ? 1. : 0. for k in 0:d-1 ]))

"""
    sample(::State [; rng])
    sample(::State, pos [; rng])
    sample(::State, op [; rng])

a random outcome of measuring the state in the basis of its sites, numbered from 0: a vector
with one outcome per site, or, given `pos`, the outcome of that site alone; or, given a
Hermitian operator placed on a few sites, one of its eigenvalues, see `probabilities`. The state
is left as it is, see `collapse`. `rng` is the random number generator, the global one by
default.

# Examples

    sample(state)
    sample(state, 3)
    sample(state, Z(1) * Z(2))
"""
function sample(state::State{Pure}; rng = Random.default_rng())
    nonzero_trace(state, "sample")
    st = orthogonalize(state.state, 1)
    st[1] /= norm(st[1])
    return sample(rng, st) .- 1
end

function sample(state::State{Pure}, pos::Int; rng = Random.default_rng())
    check_positions(state, [pos], "sample")
    nonzero_trace(state, "sample")
    st = orthogonalize(state.state, pos)
    t = st[pos] / norm(st[pos])
    ind = siteind(Pure, state.system, pos)
    # the conjugate leg, the direction the tensor takes on a charged site
    return draw(dim(ind), rand(rng)) do x
        tx = t * onehot(dag(ind) => x + 1)
        real(scalar(tx * dag(tx)))
    end
end

function sample(state::State{Mixed}, pos::Int; rng = Random.default_rng())
    check_positions(state, [pos], "sample")
    state = weak_form(state)
    l = get_left(state, pos) / nonzero_trace(state, "sample")
    r = get_right(state, pos)
    d = dim(state.system[pos])
    return draw(x -> real(scalar(l * obs_tensor(state, level_projector(d, x), pos) * r)), d,
                rand(rng))
end

function sample(state::State{Mixed}; rng = Random.default_rng())
    nonzero_trace(state, "sample")
    state = weak_form(state)
    n = length(state)
    result = Vector{Int}(undef, n)
    l = ITensor(1.)
    for pos in 1:n
        a = l * state.state[pos]
        r = get_right(state, pos)
        d = dim(state.system[pos])
        tot = real(scalar(a * identity_at(state, pos) * r))
        x = draw(i -> real(scalar(a * obs_tensor(state, level_projector(d, i), pos) * r)), d,
                 rand(rng) * tot)
        result[pos] = x
        l = a * obs_tensor(state, level_projector(d, x), pos)
        # rescaled, the probabilities to come being ratios: it would underflow on long chains
        l /= norm(l)
    end
    return result
end

function sample(state::State, op::IndexedOp{IsOperator}; rng = Random.default_rng())
    ps = probabilities(state, op)
    return first(ps[draw(x -> last(ps[x + 1]), length(ps), rand(rng)) + 1])
end

"""
    collapse(state, pos [; rng])
    collapse(state, op [; rng, limits])
    collapse(state, pos => x)
    collapse(state, op => λ [; limits])

the result of measuring site `pos`, or the Hermitian operator `op` placed on a few sites, its
probability, and the state it leaves, `(result, probability, state)`: the result is drawn with
its probability, see `probabilities`, or given as a pair, `pos => x` or `op => λ`, a result that
is not possible or of zero probability being refused, and the state is projected onto it and
normalized. The projection on several sites is truncated by `limits`, `Limits()` by default:
for a product of involutions, the state plus its image by the product, which at most doubles
the bond dimension, or quadruples it on a density matrix; for a sum of operators of one site,
an MPO whose bond carries the partial sums, as many as its dimension; otherwise the projector,
with the Jordan-Wigner strings of the fermionic sites between its sites, refused when it
leaves the state without definite charges, on a system conserving them, see `weaken`. `rng` is
the random number generator, the global one by default.

# Examples

    x, p, state = collapse(state, 3)
    λ, p, state = collapse(state, X(3))
    s, p, state = collapse(state, Z(1) * Z(2); limits = Limits(maxdim = 100))
    _, p, state = collapse(state, Z(1) * Z(2) => 1)
"""
function collapse(state::State, op::IndexedOp{IsOperator}; rng = Random.default_rng(),
                  limits::Limits = Limits())
    values, probs, project = measurement(state, op, "collapse")
    k = draw(x -> probs[x + 1], length(probs), rand(rng)) + 1
    return (result = values[k], probability = probs[k], state = project(k, limits))
end

function collapse(state::State, pos::Int; rng = Random.default_rng())
    check_positions(state, [pos], "collapse")
    r = collapse(state, Level(pos); rng)
    return (result = round(Int, r.result), probability = r.probability, state = r.state)
end

function collapse(state::State, result::Pair{<:IndexedOp{IsOperator}, <:Real};
                  limits::Limits = Limits())
    op, λ = result
    values, probs, project = measurement(state, op, "collapse")
    k = findfirst(v -> abs(v - λ) ≤ spectral_tol(values), values)
    if isnothing(k)
        error("$λ is not a result of measuring $op, whose results are $(join(values, ", "))")
    end
    if probs[k] ≤ 0
        error("collapse cannot project on the result $λ of $op, of probability $(probs[k])")
    end
    return (result = values[k], probability = probs[k], state = project(k, limits))
end

function collapse(state::State, result::Pair{Int, <:Real})
    pos, x = result
    check_positions(state, [pos], "collapse")
    r = collapse(state, Level(pos) => x)
    return (result = round(Int, r.result), probability = r.probability, state = r.state)
end
