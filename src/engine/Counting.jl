# The operators whose measurement or dephasing needs no dense matrix of their sites: products of
# involutions of one site and sums of operators of one site, each commuting with the operators
# of other sites, and the products of phases and the partial totals a sum of operators of one
# site is measured and projected with.

"""
    commuting_on(system, i, m)

whether the matrix `m` of an operator on site `i` of `system` is even, so that it commutes with
the operators of other sites, and diagonal on a system conserving something, so that its
eigenspaces keep the charges
"""
commuting_on(system::System, i::Int, m::AbstractMatrix) =
    parity_of(m, matrix(F, system[i])) == 0 &&
    (!is_charged(system) || nearly(m, Diagonal(diag(m))))

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
        if !(nearly(m, m') && nearly(m * m, one(m)) && commuting_on(system, i, m))
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
        if !commuting_on(system, i, m)
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
    phase_operators(sites, locals)

the totals of the sum that `counting` describes, twice its values, in increasing order, and the
`K` products over its sites of ``e^{2i\\pi k t/K}``, `t` twice the value on the site, for `k`
in `0:K-1`, the totals lying within `K` successive integers: their mean values give the
probabilities of the totals, their mean as gates the dephasing in the eigenbasis of the sum
"""
function phase_operators(sites, locals)
    totals = foldl((s, l) -> Set(x + t for x in s for t in first(l)), locals; init = Set(0))
    t0, t1 = extrema(totals)
    K = t1 - t0 + 1
    us = map(0:K-1) do k
        phase(t, p) = Operator{1}("Count$k", sum(cis(2π * k * u / K) * q for (u, q) in zip(t, p)))
        prod(phase(l...)(i) for (i, l) in zip(sites, locals))
    end
    return sort(collect(totals)), us
end

"""
    count_bonds(n, sites, locals, T)

for each bond of a chain of `n` sites, from the left of site 1 to the right of site `n`, the
partial totals of the sum that `counting` describes, twice its values, that the sites on its
left reach and from which those on its right can still reach `T`, in increasing order
"""
function count_bonds(n::Int, sites, locals, T::Int)
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
    return bonds
end
