# The Kraus operators of a channel, each of definite charge, for the trajectories,
# which draw one of them where a density matrix takes their sum.

"""
    kraus_operators(system, a)

the Kraus operators `K` of the channel `a` on `system`, operators placed on its sites, `a` being
``\\rho \\mapsto \\sum_K K\\rho K^\\dagger``, each of definite charge but for a gate
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

the Kraus operators of the channel `a` placed on sites, from its structure, refused when they
have no definite charge, as a trajectory draws among the operators written
"""
function site_kraus(system::System, a::AtIndex{IsChannel})
    sites = [ system[i] for i in a.index ]
    ks = generic_kraus(a.op, sites)
    if !all(k -> definite_charge(k, sites), ks)
        error("the Kraus operators of $(obs_name(a.op)) as written have no definite charge on " *
              "$(join(sites, " ⊗ ")): write the channel as a sum of gates of definite charge, " *
              "among which its trajectories are drawn")
    end
    return [ k(a.index...) for k in ks ]
end

"""
    generic_kraus(a, sites)

the Kraus operators of the channel `a` of `sites` from its structure, as written, from the Choi
matrix of a part with none
"""
function generic_kraus(a::GenericOp{IsChannel}, sites)
    if is_pure(a)
        return [kraus(a)]
    elseif a isa ScalarOp
        return [ scale(sqrt(a.coef), k) for k in generic_kraus(a.arg, sites) ]
    elseif a isa Dephase
        values, _ = dephase_projectors(a, sites)
        ps = [ Proj(a.arg => λ) for λ in values ]
        if a.kraus === :projectors
            return ps
        end
        d = length(ps)
        # the identity written as such, which a Dissipator of the channel drops
        unit = scale(1 / sqrt(d), IdentityOp{IsOperator, Generic, length(sites)}())
        return [ [unit] ; [ sum(cis(2π * j * k / d) / sqrt(d) * p for (j, p) in enumerate(ps))
                           for k in 1:d - 1 ] ]
    elseif a isa Reset
        return reset_kraus(a, only(sites))
    elseif a isa SumOp
        return reduce(vcat, [ generic_kraus(x, sites) for x in a.subs ])
    elseif a isa Union{ProdOp, IntPowOp, TensorOp}
        parts = a isa ProdOp ? [ generic_kraus(x, sites) for x in a.subs ] :
                a isa IntPowOp ? fill(generic_kraus(a.arg, sites), a.expo) :
                [ generic_kraus(x, sites[r]) for (x, r) in zip(a.subs, factor_sites(a)) ]
        combine = a isa TensorOp ? (ks -> reduce(⊗, ks)) : prod
        return vec([ combine(ks) for ks in Iterators.product(parts...) ])
    end
    return choi_kraus(a, sites)
end

"""
    reset_kraus(a, site)

the Kraus operators of the reset `a` of `site` to σ, ``\\sqrt{λ}|s⟩⟨n|`` for each eigenvector `s`
of σ, see `site_eigen`, and each basis state `n`
"""
function reset_kraus(a::Reset, site::AbstractSite)
    x = local_state(site, a.state)
    d = dim(site)
    basis(k) = [ i == k ? 1. : 0. for i in 1:d ]
    if x isa AbstractMatrix && !nearly(x, x')
        error("$(obs_name(a)) does not keep density matrices Hermitian on $site, so it has no " *
              "Kraus operators")
    end
    σ = x isa AbstractMatrix ? x / real(tr(x)) : x * x' / real(x' * x)
    es = site_eigen(σ, site)
    if isnothing(es)
        error("$(obs_name(a)) resets $site to a state of no definite charge, which a trajectory " *
              "cannot hold")
    elseif any(e -> first(e) < -rounding_tol, es)
        error("$(obs_name(a)) is not completely positive on $site, so it has no Kraus operators")
    end
    pairs = filter(e -> first(e) > rounding_tol, es)
    return vec([ Operator{1}("K$k[$(obs_name(a))]", sqrt(λ) * s * basis(n)')
                 for (k, ((λ, s), n)) in enumerate(Iterators.product(pairs, 1:d)) ])
end

# the flux, in the strongly conserved quantities `names`, of the operator `k` of definite charge on
# `sites`: its nonzero components
function strong_flux(k::GenericOp{IsOperator}, sites, names)
    m = matrix(k, sites...)
    g = basis_grades(sites)
    i = findfirst(x -> abs(x) > rounding_tol * norm(m), m)
    if isnothing(i)
        return Pair{String, Int}[]
    end
    return strong_part((g[i[1]] - g[i[2]]).q, names)
end

# the sum of the fluxes `fs`, as `strong_flux` gives them
function added_fluxes(fs)
    total = Dict{String, Int}()
    for f in fs, (nm, v) in f
        total[nm] = get(total, nm, 0) + v
    end
    return sort!([ nm => v for (nm, v) in total if v ≠ 0 ])
end

"""
    basis_grades(sites)

the charge and the parity of each basis state of `sites`, as a `Grade`, in the order of `kron`:
an element `(a, b)` of an operator changes them by `g[a] - g[b]`
"""
basis_grades(sites) =
    foldl((ga, gb) -> vec([ x + y for y in gb, x in ga ]),
          [ [ Grade(q, p) for (q, p) in zip(basis_charges(s), fermion_parities(s)) ]
            for s in sites ])

"""
    definite_charge(k, sites)

whether the operator `k` of `sites` has a definite charge, in the quantities they conserve
"""
function definite_charge(k::GenericOp{IsOperator}, sites)
    m = matrix(k, sites...)
    g = [ x.q for x in basis_grades(sites) ]
    small = rounding_tol * norm(m)
    return length(unique(g[a] - g[b] for a in axes(m, 1), b in axes(m, 2)
                         if abs(m[a, b]) > small)) ≤ 1
end

"""
    choi_kraus(a, sites)

the Kraus operators of the channel `a` of `sites`, its Choi matrix diagonalized block by block of
the charges they conserve and, unless it has coherences between them, of the fermionic parities,
as `site_eigen` does, those of weight below rounding dropped
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
    function labels_of(g)
        gs = vec([ g[k] - g[i] for k in 1:d, i in 1:d ])
        return indexin(gs, unique(gs))
    end
    coherent(l) = any(abs.(j[l .≠ permutedims(l)]) .> rounding_tol * norm(j))
    graded = basis_grades(sites)
    labels = labels_of([ x.q for x in graded ])
    if coherent(labels)
        error("$(obs_name(a)) has no definite charge on $on, so it has no Kraus operators of " *
              "definite charge")
    end
    parted = labels_of(graded)
    if !coherent(parted)
        labels = parted
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
