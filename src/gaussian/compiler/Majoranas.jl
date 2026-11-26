# Operators of a system as polynomials in Majoranas, for the Gaussian backend. The modes are those
# of the fermionic sites, in the order of the sites and, on a site, of its species (↑ before ↓ on an
# Electron), the order of the Jordan-Wigner strings of TMS. A mode of the Gaussian backend comes
# first, the ancilla: Majoranas 1 and 2, the system's from 3 on, and an odd monomial c(α) of the
# system becomes i γ₁ c(α), which keeps every state even (local/gaussian/experiences.md, « Mode
# ancillaire »).

import .GaussianKernels: MajoranaPoly, monomial_product

"""
    MajoranaLayout(system)

the Majoranas of the fermionic sites of `system`: for each site, its local Majoranas as matrices
on the site and the index of its first one in the system, the ancilla excluded; an error for a site
whose annihilation operators do not make fermionic modes.
"""
struct MajoranaLayout
    first::Vector{Int}
    majoranas::Vector{Vector{Matrix{ComplexF64}}}
    count::Int
end

# the Majoranas c + c†, i(c† - c) of the species of a site, which must satisfy {γ_a, γ_b} = 2δ_ab
# to rounding, as TMS checks its species, and span its operators; none on a site that is not fermionic
function local_majoranas(site::AbstractSite)
    species = fermion_species(site)
    if isnothing(species)
        if any(isodd, fermion_parities(site))
            throw(ArgumentError("the Gaussian backend does not accept sites $(typeof(site)), whose " *
                                "annihilation operators do not make fermionic modes"))
        end
        return Matrix{ComplexF64}[]
    end
    γs = Matrix{ComplexF64}[]
    for c in species
        m = matrix(c, site)
        push!(γs, m + m', im * (m' - m))
    end
    d = dim(site)
    clifford = all(norm(γs[a] * γs[b] + γs[b] * γs[a] - (a == b ? 2I : 0I)) <= rounding_tol * dim(site) for a in eachindex(γs), b in eachindex(γs))
    if !clifford || 2^length(species) != d
        throw(ArgumentError("the Gaussian backend does not accept sites $(typeof(site)), whose " *
                            "annihilation operators do not make fermionic modes"))
    end
    γs
end

function MajoranaLayout(system::System)
    first = Int[]
    majoranas = Vector{Matrix{ComplexF64}}[]
    count = 0
    for site in system.sites
        γs = local_majoranas(site)
        push!(first, count + 1)
        push!(majoranas, γs)
        count += length(γs)
    end
    MajoranaLayout(first, majoranas, count)
end

# the local monomials of a site with q modes, c(β) for β ⊆ 1:2q, as (β, matrix)
function local_monomials(γs::Vector{Matrix{ComplexF64}}, d::Int)
    out = Tuple{Vector{Int}, Matrix{ComplexF64}}[]
    for mask in 0:(2^length(γs) - 1)
        β = [r for r in eachindex(γs) if (mask >> (r - 1)) & 1 == 1]
        push!(out, (β, foldl(*, γs[β]; init = Matrix{ComplexF64}(I, d, d))))
    end
    out
end

# the matrix m of a site as Σ_β c_β c(β) in its local Majoranas, c_β = tr(c(β)† m)/d, the global
# indices of the Majoranas; the monomials are orthonormal for tr(A†B)/d
function local_expansion(layout::MajoranaLayout, k::Int, m::AbstractMatrix)
    γs = layout.majoranas[k]
    d = size(m, 1)
    out = Pair{Vector{Int}, ComplexF64}[]
    for (β, c) in local_monomials(γs, d)
        coef = tr(c' * m) / d
        if coef != 0
            push!(out, (β .+ (layout.first[k] - 1)) => coef)
        end
    end
    out
end

"""
    majorana_poly(system, op)

the operator `op` placed on `system` as a `MajoranaPoly`, the ancilla included; an error if it
acts on a site that is not fermionic.
"""
function majorana_poly(system::System, op::IndexedOp{IsOperator})
    layout = MajoranaLayout(system)
    terms = lower(system, op)
    # the parity F_k of each fermionic site as a monomial, and F_<k as their product
    nsites = length(system.sites)
    before = Vector{Tuple{ComplexF64, Vector{Int}}}(undef, nsites)
    acc = (one(ComplexF64), Int[])
    for k in 1:nsites
        before[k] = acc
        if !isempty(layout.majoranas[k])
            for (β, c) in local_expansion(layout, k, matrix(F, system[k]))
                s, α = monomial_product(acc[2], β)
                acc = (acc[1] * s * c, α)
            end
        end
    end
    d = Dict{Vector{Int}, ComplexF64}()
    for t in terms
        add_term!(d, system, layout, before, t)
    end
    with_ancilla(MajoranaPoly(d))
end

# the term t added to d: the factor of site k is its block followed by its strings there, b.m F_k
# (see `string_factors`); a local monomial c(β) of site k, odd, is F_<k c(β) in the global Majoranas
function add_term!(d, system::System, layout::MajoranaLayout, before, t::Term)
    blocks = Dict(b.site => b.m for b in t.blocks)
    factors = Vector{Pair{Vector{Int}, ComplexF64}}[]
    sites = Int[]
    for k in eachindex(system.sites)
        string = on_string(t.kets, k)
        if !haskey(blocks, k) && !string
            continue
        end
        if isempty(layout.majoranas[k])
            # F is the identity there, a block is not
            if haskey(blocks, k)
                throw(ArgumentError("the Gaussian backend does not yet accept operators on sites " *
                                    "$(typeof(system[k]))"))
            end
            continue
        end
        m = haskey(blocks, k) ? blocks[k] : Matrix{ComplexF64}(I, dim(system[k]), dim(system[k]))
        if string
            m = m * matrix(F, system[k])
        end
        push!(factors, local_expansion(layout, k, m))
        push!(sites, k)
    end
    for choice in Iterators.product(factors...)
        coef = ComplexF64(t.coef)
        α = Int[]
        for (k, (β, c)) in zip(sites, choice)
            coef *= c
            if isodd(length(β))
                s, α = monomial_product(α, before[k][2])
                coef *= s * before[k][1]
            end
            s, α = monomial_product(α, β)
            coef *= s
        end
        d[α] = get(d, α, 0) + coef
    end
    d
end

# the system's polynomial on the Majoranas of the Gaussian backend: indices shifted past the ancilla,
# an odd monomial c(α) made i γ₁ c(α)
with_ancilla(p::MajoranaPoly) =
    MajoranaPoly(Dict((isodd(length(α)) ? vcat(1, α .+ 2) : α .+ 2) => (isodd(length(α)) ? im : 1) * c
                      for (α, c) in p.terms))
