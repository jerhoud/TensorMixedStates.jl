# Operators of a system as polynomials in Majoranas, for the Gaussian backend. The modes are those of
# the fermionic sites, in the order of their species (↑ before ↓ on an Electron), and those of the
# sites of dimension 2 (qubits, spins 1/2), each occupied along a Pauli axis. The order of the modes, that
# of Jordan-Wigner, takes the sites block by block, a block per kind in the order the kinds first
# appear, then in the order of the sites; the parities of the sites before a mode in this order
# precede its Majoranas, which for the second block is the Klein factor that makes its odd operators
# anticommute with those of the first. A mode comes first, the ancilla: Majoranas 1 and 2, the
# system's from 3 on, and an odd monomial c(α) of the system becomes i γ₁ c(α), which keeps every
# state even (local/gaussian/experiences.md, « Mode ancillaire »).

import .GaussianKernels: MajoranaPoly, monomial_product

const pauli_x = ComplexF64[0 1; 1 0]
const pauli_y = ComplexF64[0 -im; im 0]
const pauli_z = ComplexF64[1 0; 0 -1]

pauli(v) = v[1] * pauli_x + v[2] * pauli_y + v[3] * pauli_z

"""
    MajoranaLayout(system; axes = Dict())

the Majoranas of `system`: for each site, its local Majoranas as matrices on the site, the parity of
its modes, the index of its first Majorana in the system (the ancilla excluded) and the parities of
the sites before it in the order of the modes; `axes` maps a site of dimension 2 to the Pauli, `:X`,
`:Y` or `:Z`, along which its mode is occupied (the -1 eigenstate), `:Z` by default; an error for a
site whose annihilation operators do not make fermionic modes.
"""
struct MajoranaLayout
    first::Vector{Int}
    majoranas::Vector{Vector{Matrix{ComplexF64}}}
    parity::Vector{Matrix{ComplexF64}}
    before::Vector{Tuple{ComplexF64, Vector{Int}}}
    count::Int
end

# the Majoranas c + c†, i(c† - c) of the species of a fermionic site, which must satisfy
# {γ_a, γ_b} = 2δ_ab to rounding, as TMS checks its species, and span its operators
function fermion_majoranas(site::AbstractSite, species)
    γs = Matrix{ComplexF64}[]
    for c in species
        m = matrix(c, site)
        push!(γs, m + m', im * (m' - m))
    end
    clifford = all(norm(γs[a] * γs[b] + γs[b] * γs[a] - (a == b ? 2I : 0I)) <= rounding_tol * dim(site) for a in eachindex(γs), b in eachindex(γs))
    if !clifford || 2^length(species) != dim(site)
        throw(ArgumentError("the Gaussian backend does not accept sites $(typeof(site)), whose " *
                            "annihilation operators do not make fermionic modes"))
    end
    γs
end

# the axis :X, :Y or :Z as a vector; other axes would leave rounding in the expansions, which a tolerance
# would have to remove (decided 2026-11-26: Pauli axes only until a use asks for others)
function axis_vector(axis::Symbol)
    if !(axis in (:X, :Y, :Z))
        throw(ArgumentError("an axis is :X, :Y or :Z, not :$axis"))
    end
    Float64[axis == :X, axis == :Y, axis == :Z]
end

# σ_b, σ_c for a mode occupied along the Pauli axis a: the next two Paulis in the cyclic order X, Y, Z,
# so that -iσ_bσ_c = σ_a
function qubit_majoranas(a)
    k = findfirst(==(1), a)
    e(j) = Float64[mod1(j, 3) == i for i in 1:3]
    [pauli(e(k + 1)), pauli(e(k + 2))]
end

# :fermion, :qubit or nothing, and the local Majoranas and parity of site k
function site_modes(system::System, k::Int, axes)
    site = system[k]
    species = fermion_species(site)
    if !isnothing(species)
        return :fermion, fermion_majoranas(site, species), Matrix{ComplexF64}(matrix(F, site))
    elseif any(isodd, fermion_parities(site))
        throw(ArgumentError("the Gaussian backend does not accept sites $(typeof(site)), whose " *
                            "annihilation operators do not make fermionic modes"))
    elseif dim(site) == 2
        a = axis_vector(get(axes, k, :Z))
        return :qubit, qubit_majoranas(a), pauli(a)
    end
    nothing, Matrix{ComplexF64}[], Matrix{ComplexF64}(I, dim(site), dim(site))
end

function MajoranaLayout(system::System; axes = Dict{Int, Symbol}())
    n = length(system.sites)
    kinds, majoranas, parity = Any[], Vector{Matrix{ComplexF64}}[], Matrix{ComplexF64}[]
    for k in 1:n
        kind, γs, p = site_modes(system, k, axes)
        push!(kinds, kind)
        push!(majoranas, γs)
        push!(parity, p)
    end
    for k in keys(axes)
        if !(k in 1:n) || kinds[k] != :qubit
            throw(ArgumentError("an axis is for a site of dimension 2 that is not fermionic, not for site $k"))
        end
    end
    blocks = unique(filter(!isnothing, kinds))
    order = [k for b in blocks for k in 1:n if kinds[k] == b]
    first = zeros(Int, n)
    count = 0
    for k in order
        first[k] = count + 1
        count += length(majoranas[k])
    end
    layout = MajoranaLayout(first, majoranas, parity, Vector{Tuple{ComplexF64, Vector{Int}}}(undef, n), count)
    acc = (one(ComplexF64), Int[])
    for k in 1:n
        layout.before[k] = acc
    end
    for k in order
        layout.before[k] = acc
        # the parity of a site is one monomial, the product of its Majoranas up to a phase
        expansion = local_expansion(layout, k, parity[k])
        if length(expansion) != 1
            throw(ArgumentError("the parity of the modes of site $k is not one monomial of its Majoranas"))
        end
        β, c = only(expansion)
        s, α = monomial_product(acc[2], β)
        acc = (acc[1] * s * c, α)
    end
    layout
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

# the matrix m of site k as Σ_β c_β c(β) in its local Majoranas, c_β = tr(c(β)† m)/d, with the global
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
    majorana_poly(system, op; axes = Dict())

the operator `op` placed on `system` as a `MajoranaPoly` on the Majoranas of the Gaussian backend,
the ancilla included, as it acts on its states, which are even; see `MajoranaLayout` for `axes`; an
error if it acts on a site that has no mode.
"""
function majorana_poly(system::System, op::IndexedOp{IsOperator}; axes = Dict{Int, Symbol}())
    layout = MajoranaLayout(system; axes)
    d = Dict{Vector{Int}, ComplexF64}()
    for t in lower(system, op)
        add_term!(d, system, layout, t)
    end
    on_even(with_ancilla(MajoranaPoly(d)), layout.count + 2)
end

# the term t added to d: the factor of site k is its block followed by its strings there, b.m F_k
# (see `string_factors`); a local monomial c(β) of site k, odd, is B_k c(β) in the global
# Majoranas, B_k the parities of the sites before k in the order of the modes
function add_term!(d, system::System, layout::MajoranaLayout, t::Term)
    blocks = Dict(b.site => b.m for b in t.blocks)
    factors = Vector{Pair{Vector{Int}, ComplexF64}}[]
    sites = Int[]
    for k in eachindex(system.sites)
        # the strings of TMS lie on the fermionic sites, F is the identity on the others
        string = on_string(t.kets, k) && matrix(F, system[k]) != I
        if !haskey(blocks, k) && !string
            continue
        end
        if isempty(layout.majoranas[k])
            throw(ArgumentError("the Gaussian backend accepts fermionic sites and sites of " *
                                "dimension 2, not an operator on $(typeof(system[k]))"))
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
                s, α = monomial_product(α, layout.before[k][2])
                coef *= s * layout.before[k][1]
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

# on the even states, the total parity P = (-i)^M c(1:2M) of the M modes is 1, so c(α) = c(α) P: the
# shorter of α and its complement is kept, α when they are as long
function on_even(p::MajoranaPoly, nmajoranas::Int)
    everything = collect(1:nmajoranas)
    phase = (-im)^(nmajoranas ÷ 2)
    d = Dict{Vector{Int}, ComplexF64}()
    for (α, c) in p.terms
        if 2 * length(α) > nmajoranas
            s, β = monomial_product(α, everything)
            d[β] = get(d, β, 0) + c * s * phase
        else
            d[α] = get(d, α, 0) + c
        end
    end
    MajoranaPoly(d)
end
