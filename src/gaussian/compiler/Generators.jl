# Generators of the Gaussian backend: the operators that must be quadratic in Majoranas, refused
# otherwise with the terms that are not and, for qubits, the Pauli axes that would make them so.

import .GaussianKernels: quadratic_form, degree, monomial_product

"""
    gaussian_generator(system, op; axes = Dict(), parity = nothing)

`(c, k, nmajoranas)` with `op = c + (1/4) Σ_ab k_ab γ_a γ_b` on the states of the Gaussian backend,
see `majorana_poly` and `quadratic_form`; with `parity` ±1, on the states of that physical parity when
`op` conserves it (see `in_sector`); an error naming the terms of `op` that are not quadratic, and the
Pauli axes of the qubits that would make them so when it finds some.
"""
function gaussian_generator(system::System, op::IndexedOp{IsOperator}; axes = Dict{Int, Symbol}(),
                            parity::Union{Nothing, Int} = nothing)
    layout = MajoranaLayout(system; axes)
    terms = lower(system, op)
    p = polynomial(system, layout, terms)
    if !isnothing(parity) && conserves_parity(p)
        p = in_sector(p, layout.count, parity)
    end
    if degree(p) <= 2
        c, k = quadratic_form(p, layout.count + 2)
        return c, k, layout.count + 2
    end
    bad = [α for α in keys(p.terms) if length(α) > 2]
    culprits = unique(sort(unique(b.site for b in t.blocks)) for t in terms
                      if (q = polynomial(system, layout, [t]); any(haskey(q.terms, α) for α in bad)))
    largest = maximum(abs(p.terms[α]) for α in bad)
    shown = join(("sites " * join(ks, ", ") for ks in first(culprits, 3)), "; ")
    more = length(culprits) > 3 ? " and $(length(culprits) - 3) others" : ""
    message = "the Gaussian backend needs an operator quadratic in Majoranas, its terms on $shown$more " *
              "are of degree $(maximum(length, bad)) (coefficients up to $(Printf.@sprintf("%.2g", largest)))."
    throw(ArgumentError(message * qubit_advice(system, layout, terms, axes)))
end

# whether the polynomial conserves the physical parity: no monomial has one Majorana of the ancilla
# without the other, as the odd operators of the system give (γ₁γ₂ is the parity of the ancilla, which
# the reduction on the even states brings in)
conserves_parity(p::MajoranaPoly) = !any(α -> (1 in α) ⊻ (2 in α), keys(p.terms))

# p on the states whose physical parity P = (-i)^n c(3:2n+2) is σ, n = count/2 modes, p conserving it:
# there γ₁γ₂ = iσ (the parity of the ancilla is that of the system), and c(α) = σ c(α) P =
# σ (-i)^n c(α) c(3:2n+2) for α in the system, the shorter kept (the term closing a ring of qubits, the
# parity times a quadratic term, becomes quadratic), α when they are as long
function in_sector(p::MajoranaPoly, count::Int, σ::Int)
    system = collect(3:count+2)
    phase = σ * (-im)^(count ÷ 2)
    d = Dict{Vector{Int}, ComplexF64}()
    for (α, c) in p.terms
        if 1 in α
            α, c = α[3:end], c * im * σ
        end
        if 2 * length(α) > count
            s, β = monomial_product(α, system)
            d[β] = get(d, β, 0) + c * s * phase
        else
            d[α] = get(d, α, 0) + c
        end
    end
    MajoranaPoly(d)
end

# the qubits of the system, in the order of the modes
qubit_sites(system::System, layout::MajoranaLayout) =
    [k for k in layout.order if isnothing(fermion_species(system[k]))]

const pauli_matrices = ((:X, pauli_x), (:Y, pauli_y), (:Z, pauli_z))

# the terms split by their Paulis on the qubits: (site => Pauli) sorted by site, and the terms with
# these Paulis, the other blocks unchanged; distinct Paulis give distinct monomials, whatever the axes
function pauli_groups(qubits, terms::Vector{Term})
    groups = Dict{Vector{Pair{Int, Symbol}}, Vector{Term}}()
    for t in terms
        qblocks = [b for b in t.blocks if b.site in qubits]
        others = [b for b in t.blocks if !(b.site in qubits)]
        options = [[(P, tr(σ * b.m) / 2) for (P, σ) in ((:I, I), pauli_matrices...) if tr(σ * b.m) != 0]
                   for b in qblocks]
        for choice in Iterators.product(options...)
            key = [b.site => P for (b, (P, _)) in zip(qblocks, choice) if P != :I]
            blocks = vcat(others, [Block(b.site, pauli(axis_vector(P)), 0)
                                   for (b, (P, _)) in zip(qblocks, choice) if P != :I])
            sort!(blocks; by = b -> b.site)
            coef = t.coef * prod(c for (_, c) in choice; init = one(ComplexF64))
            push!(get!(groups, key, Term[]), Term(coef, blocks, t.kets, t.bras))
        end
    end
    groups
end

# a Pauli other than P, for a qubit whose axis must differ from it: the degree does not depend on which
other_pauli(P::Symbol) = P == :X ? :Y : :X

# the patterns of the group, as sites => whether the axis is off its Pauli, that make it quadratic, by
# compiling it alone; nothing when its support is too large to try them all
function quadratic_patterns(system::System, key, ts::Vector{Term}; largest = 12)
    if length(key) > largest
        return nothing
    end
    patterns = BitVector[]
    for mask in 0:(2^length(key) - 1)
        odd = BitVector([(mask >> (i - 1)) & 1 == 1 for i in eachindex(key)])
        axes = Dict(k => (o ? other_pauli(P) : P) for ((k, P), o) in zip(key, odd))
        if degree(polynomial(system, MajoranaLayout(system; axes), ts)) <= 2
            push!(patterns, odd)
        end
    end
    patterns
end

# the axes of the qubits, as sets, that satisfy one pattern of each group, by backtracking; nothing if
# there are none, :limit if the search stopped
function axes_search(constraints, domains, i::Int, nodes::Ref{Int}; limit = 10^5)
    if i > length(constraints)
        return domains
    end
    key, patterns = constraints[i]
    for odd in patterns
        nodes[] += 1
        if nodes[] > limit
            return :limit
        end
        d = Dict(q => copy(s) for (q, s) in domains)
        for ((k, P), o) in zip(key, odd)
            if o
                setdiff!(d[k], [P])
            else
                intersect!(d[k], [P])
            end
        end
        if all(!isempty, values(d))
            r = axes_search(constraints, d, i + 1, nodes; limit)
            if !isnothing(r)
                return r
            end
        end
    end
    nothing
end

# the Pauli axes of the qubits that make the generator quadratic: each group of terms with the same
# Paulis on the qubits must be quadratic, which depends only on whether the axis of each of its qubits
# is its Pauli or not; the patterns that are are found by compiling, then one per group chosen by
# backtracking; the axes found are checked by compiling again
function qubit_advice(system::System, layout::MajoranaLayout, terms::Vector{Term}, axes)
    qubits = qubit_sites(system, layout)
    if isempty(qubits)
        return ""
    end
    constraints = Tuple{Vector{Pair{Int, Symbol}}, Vector{BitVector}}[]
    complete = true
    gaps = Tuple{Int, Int}[]
    for (key, ts) in pauli_groups(qubits, terms)
        patterns = quadratic_patterns(system, key, ts)
        if isnothing(patterns)
            complete = false
        elseif !isempty(key) || isempty(patterns)
            push!(constraints, (key, patterns))
        end
        sites = sort([k for (k, _) in key]; by = k -> layout.first[k])
        # a gap matters for a group that no axis makes quadratic
        if length(sites) >= 2 && !isnothing(patterns) && isempty(patterns)
            between = [q for q in qubits if layout.first[sites[1]] < layout.first[q] < layout.first[sites[end]]]
            if !all(q in sites for q in between)
                push!(gaps, (sites[1], sites[end]))
            end
        end
    end
    sort!(constraints; by = c -> length(c[2]))
    nodes = Ref(0)
    found = axes_search(constraints, Dict(q => Set([:X, :Y, :Z]) for q in qubits), 1, nodes)
    if found isa Dict
        proposed = Dict(q => (get(axes, q, :Z) in found[q] ? get(axes, q, :Z) : first(sort(collect(found[q]))))
                        for q in qubits)
        if degree(polynomial(system, MajoranaLayout(system; axes = proposed), terms)) <= 2
            return " The axes " * join(("$q => :$(proposed[q])" for q in qubits), ", ") * " make it quadratic."
        end
    end
    advice = complete && found !== :limit ? " No Pauli axis of the qubits makes it quadratic." :
             " No Pauli axis of the qubits found makes it quadratic (search limited)."
    if !isempty(gaps)
        q1, q2 = first(gaps)
        advice *= " Qubits $q1 and $q2 are coupled across others: order the sites so that they are neighbours."
    end
    advice
end
