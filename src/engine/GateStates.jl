# The states built by applying gates to a product state: the dimer states, the Slater
# determinant of given orbitals, and the Fermi sea of a quadratic Hamiltonian, the determinant
# of its lowest orbitals.

export dimer_state, slater_state, fermi_sea

"""
    fill_others!(states, rest, others, words)

the local states `states` with those of the sites `rest` set to `others`, one local state for
all of them or a vector of one for each, which these sites then need; `words` says what they
are, for one site, for several, and after "the sites"
"""
function fill_others!(states::Vector{Any}, rest::Vector{Int}, others, words::NTuple{3, String})
    if isempty(rest)
        return states
    end
    one, many, which = words
    if isnothing(others)
        error(length(rest) == 1 ? "site $(only(rest)) $one: give its state with others" :
                                  "sites $(join(rest, ", ")) $many: give their state with others")
    end
    local_states = others isa AbstractVector && !(eltype(others) <: Number) ? others :
                   fill(others, length(rest))
    if length(local_states) ≠ length(rest)
        error("others gives $(length(local_states)) states for the $(length(rest)) sites $which")
    end
    states[rest] = local_states
    return states
end

"""
    dimer_state(system, pairs, a, b; others, limits = Limits())

the state of a singlet ``(|ab\\rangle - |ba\\rangle)/\\sqrt{2}`` on each pair `(i, j)` of
`pairs`, `a` on `i` in the first term, of the local pure states `a` and `b`, given by their
names or their vectors and orthonormal: `"Up"` and `"Dn"` for spins one half. The sites in no
pair take the local state `others`, one for all of them or a vector of one for each, which
they then need. On fermionic sites, when `a` and `b` are of different parities, ``|ba\\rangle``
is ``|ab\\rangle`` with the fermion moved from one site of the pair to the other by a hopping
``c^\\dagger c``, whose string gives it the sign of the occupied sites between: the singlet of
fermions, which differs from that of the products of local vectors when the pair spans
occupied sites.

The pairs need not be neighbours: pairs of neighbours give the Majumdar-Ghosh state, of bond
dimension 2, and nested pairs `(i, n + 1 - i)` the rainbow state, of bond dimension
``2^{n/2}`` in the middle. `limits` constrains the truncations made as the singlets are
applied.

# Examples

    dimer_state(System(10, Qubit()), [ (i, i + 1) for i in 1:2:9 ], "Up", "Dn")
    dimer_state(System(9, Spin(1/2)), [ (i, i + 1) for i in 1:2:7 ], "1/2", "-1/2";
                others = "1/2")
"""
function dimer_state(system::System, pairs::AbstractVector{Tuple{Int, Int}}, a, b;
                     others = nothing, limits::Limits = Limits())
    n = length(system)
    paired = reduce(vcat, [ [i, j] for (i, j) in pairs ]; init = Int[])
    if !allunique(paired)
        error("a site is in two pairs of $pairs")
    end
    if any(k -> !(1 ≤ k ≤ n), paired)
        error("the pairs $pairs reach beyond the $n sites of the system")
    end
    states = Vector{Any}(undef, n)
    for (i, j) in pairs
        states[i], states[j] = a, b
    end
    fill_others!(states, setdiff(1:n, paired), others,
                 ("is in no pair", "are in no pair", "in no pair"))
    st = State{Pure}(system, states)
    gates = map(pairs) do (i, j)
        va, vb = local_vector(system[i], a), local_vector(system[i], b)
        wa, wb = local_vector(system[j], a), local_vector(system[j], b)
        if norm([va vb]' * [va vb] - I) > 1e-12 || norm([wa wb]' * [wa wb] - I) > 1e-12
            error("a singlet of $a and $b needs them orthonormal")
        end
        # |ab⟩ to the singlet, |ba⟩ to the triplet of no charge, the identity elsewhere
        x, y = kron(va, wb), kron(vb, wa)
        U = I + ((x - y) / sqrt(2) - x) * x' + ((x + y) / sqrt(2) - y) * y'
        Operator{2}("Singlet", U)(i, j)
    end
    return isempty(gates) ? st : apply(prod(gates), st; limits)
end

"""
    species_sites(system)

the species of fermions of `system`, in the order they first appear, as pairs `c => sites` of
the sites holding each
"""
function species_sites(system::System)
    result = Pair{GenericOp{IsOperator, 1}, Vector{Int}}[]
    for k in 1:length(system), c in something(fermion_species(system[k]), ())
        i = findfirst(p -> first(p) == c, result)
        if isnothing(i)
            push!(result, c => [k])
        else
            push!(last(result[i]), k)
        end
    end
    if isempty(result)
        error("the system holds no free fermions, see TensorMixedStates.fermion_species")
    end
    return result
end

"""
    vacuum(site)

the state of `site` that every annihilation operator of its species of fermions sends to zero,
see `fermion_species`, refused when it is not unique
"""
function vacuum(site::AbstractSite)
    v = nullspace(reduce(vcat, [ matrix(c, site) for c in fermion_species(site) ]))
    if size(v, 2) ≠ 1
        error("$site has $(size(v, 2)) states that its fermion species " *
              "$(fermion_species(site)) all annihilate, where a single vacuum is needed")
    end
    return vec(v)
end

"""
    reference_states(system, others)

the local states of the product state of no fermion: the vacuum on the sites holding fermions,
see `vacuum`, and on the others `others`, one local state for all of them or a vector of one
for each
"""
function reference_states(system::System, others)
    n = length(system)
    holds(k) = !isempty(something(fermion_species(system[k]), ()))
    states = Any[ holds(k) ? vacuum(system[k]) : nothing for k in 1:n ]
    return fill_others!(states, [ k for k in 1:n if !holds(k) ], others,
                        ("holds no fermions", "hold no fermions", "holding no fermions"))
end

"""
    givens_circuit(orbitals)

the circuit of rotations of neighbouring modes taking the Slater determinant of the orthonormal
columns of `orbitals` to a product state, after Fishman and White (2015): the rotations
`(j, W)`, `W` a 2×2 unitary on the modes `j` and `j + 1`, in the order they are found, and the
occupation of each mode once they are applied. The modes are taken in turn: on the narrowest
window of the next ones that the precision allows, the eigenvector of the correlations nearest
to an occupation 0 or 1 is brought onto the first mode, which then drops out.
"""
function givens_circuit(orbitals::AbstractMatrix)
    n = size(orbitals, 1)
    # ⟨c†_j c_i⟩, which the rotation d = U c takes to U M U†
    M = Matrix{ComplexF64}(orbitals * orbitals')
    rotations = Tuple{Int, Matrix{ComplexF64}}[]
    occupations = zeros(Int, n)
    for i in 1:n
        # the window grows until an eigenvalue is 0 or 1 to 1e-12, which the whole rest of the
        # modes reaches, the correlations of a determinant being a projector
        λ, v = 0., ComplexF64[]
        for stop in i:n
            e = eigen(Hermitian(M[i:stop, i:stop]))
            k = argmax(abs.(e.values .- 0.5))
            λ, v = e.values[k], ComplexF64.(e.vectors[:, k])
            if min(λ, 1 - λ) ≤ 1e-12
                break
            end
        end
        occupations[i] = λ > 0.5 ? 1 : 0
        # v brought onto the first mode of the window, its components cancelled from the last
        for k in length(v):-1:2
            a, b = v[k-1], v[k]
            if iszero(b)
                continue
            end
            r = hypot(abs(a), abs(b))
            W = [conj(a) conj(b) ; -b a] / r
            v[k-1], v[k] = r, 0
            j = i + k - 2
            M[[j, j + 1], :] = W * M[[j, j + 1], :]
            M[:, [j, j + 1]] = M[:, [j, j + 1]] * W'
            push!(rotations, (j, W))
        end
    end
    return rotations, occupations
end

"""
    mode_gate(c, W, a, b)

the gate on the sites `a < b` rotating their modes of annihilation operator `c` by the 2×2
unitary `W`: ``\\mathcal{G}\\, c_p^\\dagger \\mathcal{G}^\\dagger = \\sum_q W_{qp}
c_q^\\dagger``. `a` and `b` need not be neighbours, the descent laying the strings between them.
"""
function mode_gate(c::GenericOp{IsOperator, 1}, W::AbstractMatrix, a::Int, b::Int)
    K = log(W)
    # c†_b c_a placed on (a, b) is -c_a c†_b, the tensor product being ordered by sites
    return exp(K[1, 1] * ((dag(c) * c) ⊗ Id) + K[2, 2] * (Id ⊗ (dag(c) * c)) +
               K[1, 2] * (dag(c) ⊗ c) - K[2, 1] * (c ⊗ dag(c)))(a, b)
end

"""
    slater_state(system, orbitals...; others, limits = Limits())

the Slater determinant of `orbitals`, one matrix ``\\Phi`` for each species of fermions of the
sites, in the order the species first appear along them (`C` for a `Fermion`, `Cup` then `Cdn`
for an `Electron`). The matrix of a species has a row for each site holding it, in their
order, and a column for each fermion of that species: column `k` holds the amplitudes of
orbital `k`, the columns being orthonormal. The state is
``\\prod_k \\left(\\sum_i \\Phi_{ik} c_i^\\dagger\\right) |0\\rangle``, up to a global
phase, so that ``\\langle c_i^\\dagger c_j \\rangle = \\sum_k \\overline{\\Phi_{ik}}
\\Phi_{jk}``, ``|0\\rangle`` being the vacuum of every site holding fermions.

The sites holding no fermions take the local state `others`, one for all of them or a vector of
one for each, which they then need.

It is built by a circuit of rotations of neighbouring modes (Fishman and White, 2015), exact to
1e-12 on the correlations; `limits` constrains the truncations made as it is applied, the
entanglement of a determinant growing with the size of the system. The sites may conserve the
number of fermions of each species, weakly or strongly.

# Examples

    slater_state(System(10, Fermion()), orbitals)
    slater_state(System(10, Electron(conserve = (Ntot, 2Sz))), up, down;
                 limits = Limits(cutoff = 1e-12))
    slater_state(System([Fermion(), Fermion(), Qubit(), Fermion()]), orbitals; others = "Up")
"""
function slater_state(system::System, orbitals::AbstractMatrix...; others = nothing,
                      limits::Limits = Limits())
    species = species_sites(system)
    if length(orbitals) ≠ length(species)
        error("these sites hold $(length(species)) species of fermions, which take as many " *
              "matrices of orbitals, $(join(first.(species), " and ")) in that order")
    end
    for ((c, sites), o) in zip(species, orbitals)
        if size(o, 1) ≠ length(sites)
            error("the orbitals of $c have $(size(o, 1)) amplitudes, for $(length(sites)) " *
                  "sites holding it")
        end
        if norm(o' * o - I) > 1e-10
            error("the orbitals of $c are not orthonormal")
        end
    end
    circuits = [ givens_circuit(o) for o in orbitals ]
    vectors = map(enumerate(reference_states(system, others))) do (k, v)
        for ((c, sites), (_, occupations)) in zip(species, circuits)
            p = findfirst(==(k), sites)
            if !isnothing(p) && occupations[p] == 1
                v = matrix(dag(c), system[k]) * v
            end
        end
        v
    end
    st = State{Pure}(system, vectors)
    # the determinant is the product state rotated back: the rotations found last act first
    gates = [ mode_gate(c, W', sites[j], sites[j+1])
              for ((c, sites), (rotations, _)) in zip(species, circuits) for (j, W) in rotations ]
    return isempty(gates) ? st : apply(prod(gates), st; limits)
end

"""
    one_body(system, hamiltonian, species, references)

the constant `E0` and, for each species of `species`, the matrix `h` of `hamiltonian` read as
``E_0 + \\sum_{ij} h_{ij} c_i^\\dagger c_j``, from the product state of no fermion, of local
states `references`, and the states of one fermion added to it: ``E_0 = \\langle 0 | H | 0
\\rangle`` and ``h_{ij} = \\langle 1_i | H | 1_j \\rangle - E_0 \\delta_{ij}``. Whether the
Hamiltonian is of that form is left to the caller.
"""
function one_body(system::System, hamiltonian::IndexedOp{IsOperator}, species, references)
    # without charges: the MPO is contracted between states of zero and one fermion, in the
    # basis of the vectors of each site
    plain = weaken(system, ())
    n = length(plain)
    zero = [ local_vector(plain[k], v) for (k, v) in enumerate(references) ]
    w = make_mpo(State{Pure}(plain, zero), hamiltonian)
    s = siteinds(Pure, plain)
    # the matrix of the MPO on site k between the local vectors a and b, a row on the first
    # site and a column on the last
    function transfer(k, a, b)
        t = w[k] * ITensor(conj(a), s[k]') * ITensor(b, s[k])
        if n == 1
            return reshape([scalar(t)], 1, 1)
        elseif k == 1
            return reshape(Array(t, commonind(w[1], w[2])), 1, :)
        elseif k == n
            return reshape(Array(t, commonind(w[n-1], w[n])), :, 1)
        end
        return Array(t, commonind(w[k-1], w[k]), commonind(w[k], w[k+1]))
    end
    empty = [ transfer(k, zero[k], zero[k]) for k in 1:n ]
    # the products of the empty sites on the left of a site and on its right
    left = accumulate(*, empty; init = ones(1, 1))
    right = reverse(accumulate((a, b) -> b * a, reverse(empty); init = ones(1, 1)))
    e0 = only(last(left))
    on_left(k) = k == 1 ? ones(1, 1) : left[k-1]
    on_right(k) = k == n ? ones(1, 1) : right[k+1]
    hs = map(species) do (c, sites)
        one = Dict(k => matrix(dag(c), plain[k]) * zero[k] for k in sites)
        m = length(sites)
        h = zeros(ComplexF64, m, m)
        for (p, a) in enumerate(sites)
            h[p, p] = only(on_left(a) * transfer(a, one[a], one[a]) * on_right(a)) - e0
            # the fermion on the bra side at a and on the ket side at b > a, then the other way
            up = on_left(a) * transfer(a, one[a], zero[a])
            down = on_left(a) * transfer(a, zero[a], one[a])
            for k in a+1:n
                q = findfirst(==(k), sites)
                if !isnothing(q)
                    h[p, q] = only(up * transfer(k, zero[k], one[k]) * on_right(k))
                    h[q, p] = only(down * transfer(k, one[k], zero[k]) * on_right(k))
                end
                up *= empty[k]
                down *= empty[k]
            end
        end
        h
    end
    return e0, hs
end

"""
    fermi_sea(system, hamiltonian, nparticles...; others, limits = Limits())

the ground state of `nparticles` free fermions of each species of the sites, in the order they
first appear, under `hamiltonian`, a quadratic Hamiltonian conserving the number of each,
``E_0 + \\sum_{ij} h_{ij} c_i^\\dagger c_j`` for each species: the Slater determinant of the
`nparticles` lowest orbitals of `h`, see `slater_state`, whose `others` and `limits` it takes.

A Hamiltonian of another form is refused, as one with an interaction, a term mixing species or
a term acting on the sites holding no fermions, and so is a number of fermions filling a
degenerate level in part: the message gives the nearest numbers that fill it, and a small term
added to the Hamiltonian lifts the degeneracy.

# Examples

    sea = fermi_sea(System(20, Fermion()), -sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i)
                                                for i in 1:19), 7)
    ring = -sum(dag(C)(i) * C(j) + dag(C)(j) * C(i) for (i, j) in circle_graph(10))
    fermi_sea(System(10, Fermion()), ring + 1e-3 * N(1), 4)    # 4 alone fills a level in part
    hop(c) = -sum(dag(c)(i) * c(i + 1) + dag(c)(i + 1) * c(i) for i in 1:7)
    fermi_sea(System(8, Electron()), hop(Cup) + hop(Cdn), 3, 3)
"""
function fermi_sea(system::System, hamiltonian::IndexedOp{IsOperator}, nparticles::Int...;
                   others = nothing, limits::Limits = Limits())
    species = species_sites(system)
    if length(nparticles) ≠ length(species)
        error("these sites hold $(length(species)) species of fermions, which take as many " *
              "numbers of fermions, $(join(first.(species), " and ")) in that order")
    end
    references = reference_states(system, others)
    e0, hs = one_body(system, hamiltonian, species, references)
    # rebuilt from what the states of no fermion and of one see of it, so differing by the terms
    # they do not see. Compared by the norms of MPOs: `≈` compares how operators are written,
    # and tells N(1) from dag(C)(1) * C(1)
    rebuilt = e0 * Id(1) + sum(h[p, q] * dag(c)(sites[p]) * c(sites[q])
                               for ((c, sites), h) in zip(species, hs)
                               for p in eachindex(sites), q in eachindex(sites)
                               if !iszero(h[p, q]); init = 0 * Id(1))
    plain = State{Pure}(weaken(system, ()), references)
    if norm(make_mpo(plain, hamiltonian - rebuilt)) >
       1e-10 * norm(make_mpo(plain, hamiltonian))
        error("the fermions of $(join(first.(species), " and ")) are not free under this " *
              "hamiltonian: it is not quadratic, does not conserve the number of fermions of " *
              "each species, or acts on the sites holding none")
    end
    orbitals = map(zip(species, hs, nparticles)) do ((c, sites), h, m)
        if !(0 ≤ m ≤ length(sites))
            error("$m fermions of $c do not fit on the $(length(sites)) sites holding it")
        end
        if norm(h - h') > 1e-10 * max(1, norm(h))
            error("the hamiltonian is not hermitian")
        end
        e = eigen(Hermitian(h))
        ε = e.values
        tol = 1e-10 * max(1, ε[end] - ε[1])
        closed(k) = k == 0 || k == length(sites) || ε[k+1] - ε[k] > tol
        if !closed(m)
            below = findlast(closed, 0:m) - 1
            above = m + findfirst(closed, m:length(sites)) - 1
            error("$m fermions of $c fill a degenerate level in part, which has no single " *
                  "ground state: take $below or $above, or lift the degeneracy by a small term")
        end
        e.vectors[:, 1:m]
    end
    return slater_state(system, orbitals...; others, limits)
end
