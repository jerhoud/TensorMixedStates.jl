# The Kraus operators of a channel, each of definite charge and parity, for the trajectories,
# which draw one of them where a density matrix takes their sum.

"""
    kraus_operators(system, a)

the Kraus operators `K` of the channel `a` on `system`, operators placed on its sites, `a` being
``\\rho \\mapsto \\sum_K K\\rho K^\\dagger``, each of definite charge and parity but for a gate
placed as a whole, as `Gate(X(1) + X(2))`
"""
function kraus_operators(system::System, a::IndexedOp{IsChannel})
    if a isa AtIndex
        return site_kraus(system, a)
    elseif a isa ScalarOp
        return [ scale(sqrt(a.coef), k) for k in kraus_operators(system, a.arg) ]
    elseif a isa Union{SumOp, ProdOp}
        g = gathered_channel(a)
        if !isnothing(g)
            return site_kraus(system, g)
        end
        parts = [ kraus_operators(system, x) for x in a.subs ]
        return a isa SumOp ? reduce(vcat, parts) :
               vec([ prod(ks) for ks in Iterators.product(parts...) ])
    elseif is_pure(a)
        return [kraus(a)]
    end
    error("cannot find the Kraus operators of $(obs_name(a))")
end

"""
    gathered_channel(a)

the sum or the product `a` of channels placed on the same sites as one channel placed there,
`nothing` when they lie on different ones
"""
function gathered_channel(a::Union{SumOp, ProdOp})
    xs = map(scalararg, a.subs)
    if !all(x -> x isa AtIndex && x.index == xs[1].index, xs)
        return nothing
    end
    ops = [ scale(scalarcoef(x), scalararg(x).op) for x in a.subs ]
    return (a isa SumOp ? sum(ops) : prod(ops))(xs[1].index...)
end

"""
    site_kraus(system, a)

the Kraus operators of the channel `a` placed on sites, from its structure when they have a
definite charge and parity, from its Choi matrix otherwise
"""
function site_kraus(system::System, a::AtIndex{IsChannel})
    sites = [ system[i] for i in a.index ]
    ks = generic_kraus(a.op, sites)
    if !all(k -> definite_grade(k, sites), ks)
        ks = choi_kraus(a.op, sites)
    end
    return [ k(a.index...) for k in ks ]
end

"""
    generic_kraus(a, sites)

the Kraus operators of the channel `a` of `sites` from its structure, from the Choi matrix of a
part with none or with more of them than the dimension of its superoperator
"""
function generic_kraus(a::GenericOp{IsChannel}, sites)
    most = prod(dim, sites)^2
    if is_pure(a)
        return [kraus(a)]
    elseif a isa ScalarOp
        return [ scale(sqrt(a.coef), k) for k in generic_kraus(a.arg, sites) ]
    elseif a isa Dephase
        values, _ = eigenspaces(matrix(a.arg, sites...), "Dephase", a.arg)
        return [ Proj(a.arg => λ) for λ in values ]
    elseif a isa SumOp
        parts = [ generic_kraus(x, sites) for x in a.subs ]
        return sum(length, parts) > most ? choi_kraus(a, sites) : reduce(vcat, parts)
    elseif a isa Union{ProdOp, IntPowOp, TensorOp}
        parts = a isa ProdOp ? [ generic_kraus(x, sites) for x in a.subs ] :
                a isa IntPowOp ? fill(generic_kraus(a.arg, sites), a.expo) :
                [ generic_kraus(x, sites[r]) for (x, r) in zip(a.subs, factor_sites(a)) ]
        if prod(float ∘ length, parts) > most
            return choi_kraus(a, sites)
        end
        combine = a isa TensorOp ? (ks -> reduce(⊗, ks)) : prod
        return vec([ combine(ks) for ks in Iterators.product(parts...) ])
    end
    return choi_kraus(a, sites)
end

"""
    basis_grades(sites)

the charge and the parity of each basis state of `sites`, in the order of `kron`
"""
basis_grades(sites) =
    foldl((ga, gb) -> vec([ (x[1] + y[1], x[2] + y[2]) for y in gb, x in ga ]),
          [ collect(zip(basis_charges(s), fermion_parities(s))) for s in sites ])

# the charge and the parity an element (a, b) of an operator changes
element_grade(g, a, b) = (g[a][1] - g[b][1], mod(g[a][2] + g[b][2], 2))

"""
    definite_grade(k, sites)

whether the operator `k` of `sites` has a definite charge and a definite parity
"""
function definite_grade(k::GenericOp{IsOperator}, sites)
    m = matrix(k, sites...)
    g = basis_grades(sites)
    small = rounding_tol * norm(m)
    return length(unique(element_grade(g, a, b) for a in axes(m, 1), b in axes(m, 2)
                         if abs(m[a, b]) > small)) ≤ 1
end

"""
    choi_kraus(a, sites)

the Kraus operators of the channel `a` of `sites`, its Choi matrix diagonalized block by block of
charge and parity, those of weight below rounding dropped
"""
function choi_kraus(a::GenericOp{IsChannel}, sites)
    check_dense(a, sites)
    dims = map(dim, sites)
    d = prod(dims)
    o = invperm(mixed_order(dims))
    # the matrix on vec(ρ), then J[(k, i), (l, j)] = a(|i⟩⟨j|)[k, l]
    s = matrix(a, sites...)[o, o]
    j = reshape(permutedims(reshape(s, d, d, d, d), (1, 3, 2, 4)), d^2, d^2)
    on = join(sites, " ⊗ ")
    if !nearly(j, j')
        error("$(obs_name(a)) does not keep density matrices Hermitian on $on, so it has no " *
              "Kraus operators")
    end
    g = basis_grades(sites)
    grades = vec([ element_grade(g, k, i) for k in 1:d, i in 1:d ])
    labels = indexin(grades, unique(grades))
    if any(abs.(j[labels .≠ permutedims(labels)]) .> rounding_tol * norm(j))
        error("$(obs_name(a)) has no definite charge on $on, so it has no Kraus operators of " *
              "definite charge")
    end
    blocks = [ findall(==(l), labels) for l in unique(labels) ]
    es = [ eigen(Hermitian(j[b, b])) for b in blocks ]
    top = maximum(e -> maximum(e.values), es)
    if any(e -> minimum(e.values) < -rounding_tol * top, es)
        error("$(obs_name(a)) is not completely positive on $on, so it has no Kraus operators")
    end
    ks = Operator{length(sites)}[]
    for (b, e) in zip(blocks, es), (λ, u) in zip(e.values, eachcol(e.vectors))
        if λ > rounding_tol * top
            v = zeros(eltype(u), d^2)
            v[b] = sqrt(λ) * u
            push!(ks, Operator{length(sites)}("K$(length(ks) + 1)[$(obs_name(a))]",
                                               reshape(v, d, d)))
        end
    end
    return ks
end
