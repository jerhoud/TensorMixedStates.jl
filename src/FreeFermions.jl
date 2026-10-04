# The states of free fermions: the Slater determinant of given orbitals, built by a circuit of
# rotations of neighbouring modes applied to a product state, and the Fermi sea of a quadratic
# hamiltonian, the determinant of its lowest orbitals.

export slater_state, fermi_sea

"""
    fermion_species(site)

the annihilation operators of the species of fermions a site holds, a mode of each: `(C,)` for a
`Fermion`, `(Cup, Cdn)` for an `Electron`, and none for a site holding no fermions, the default.
`slater_state` and `fermi_sea` read it off the sites they are given, and a site type of one's
own holding free fermions gets them by a method of it.
"""
fermion_species(::AbstractSite) = ()

"""
    species_sites(system)

the species of fermions of the sites of `system`, see `fermion_species`, in the order they
first appear, each with the sites holding it, its modes, as pairs `c => sites`
"""
function species_sites(system::System)
    result = Pair{SimpleOp, Vector{Int}}[]
    for k in 1:length(system), c in fermion_species(system[k])
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
    reference_states(system, others)

the local states of the product state of no fermion: `"0"` on the sites holding fermions, and
on the others `others`, one local state for all of them or a vector of one for each
"""
function reference_states(system::System, others)
    rest = [ k for k in 1:length(system) if isempty(fermion_species(system[k])) ]
    if isempty(rest)
        return fill!(Vector{Any}(undef, length(system)), "0")
    end
    if isnothing(others)
        which = length(rest) == 1 ? "site $(only(rest)) holds no fermions: give its" :
                                    "sites $(join(rest, ", ")) hold no fermions: give their"
        error("$which state with others")
    end
    local_states = others isa AbstractVector && !(eltype(others) <: Number) ? others :
                   fill(others, length(rest))
    if length(local_states) ≠ length(rest)
        error("others gives $(length(local_states)) states for the $(length(rest)) sites " *
              "holding no fermions")
    end
    result = Vector{Any}(undef, length(system))
    fill!(result, "0")
    result[rest] = local_states
    return result
end

"""
    givens_circuit(orbitals)

the circuit of rotations of neighbouring modes taking the Slater determinant of `orbitals`,
whose columns are orthonormal, to a product state, after Fishman and White (2015): the
rotations `(j, W)`, the unitary `W` of 2×2 acting on the modes `j` and `j + 1`, in the order
they are found, and the occupation of each mode once they are all applied.

The modes are taken one after the other. On a window of the next ones, as narrow as the
precision allows, the eigenvector of the correlations whose eigenvalue is nearest to 0 or 1 is
brought onto the first mode of the window, which is then empty or occupied and drops out.
"""
function givens_circuit(orbitals::AbstractMatrix)
    n = size(orbitals, 1)
    # ⟨c†_j c_i⟩, which the rotation d = U c takes to U M U†
    M = Matrix{ComplexF64}(orbitals * orbitals')
    rotations = Tuple{Int, Matrix{ComplexF64}}[]
    occupations = zeros(Int, n)
    for i in 1:n
        # the window grows until its eigenvalue is a pure occupation to 1e-12, which the full
        # rest of the modes always gives, the correlations of a determinant being a projector
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

the gate on the sites `a < b` that rotates their modes of annihilation operator `c` by `W`, a
unitary of 2×2: ``\\mathcal{G}\\, c_p^\\dagger \\mathcal{G}^\\dagger = \\sum_q W_{qp}
c_q^\\dagger``. It is the exponential of ``\\sum_{pq} K_{pq} c_p^\\dagger c_q`` for
``K = \\log W``, an even operator of two sites, which `apply` takes on sites apart as well.
"""
function mode_gate(c::SimpleOp, W::AbstractMatrix, a::Int, b::Int)
    K = log(W)
    # c†_b c_a placed on (a, b) is -c_a c†_b, the tensor product being ordered by sites
    return exp(K[1, 1] * ((dag(c) * c) ⊗ Id) + K[2, 2] * (Id ⊗ (dag(c) * c)) +
               K[1, 2] * (dag(c) ⊗ c) - K[2, 1] * (c ⊗ dag(c)))(a, b)
end

"""
    slater_state(system, orbitals...; others, limits = Limits())

the Slater determinant of `orbitals`, one matrix for each species of fermions of the sites, as
`fermion_species` gives them, `C` for a `Fermion`, `Cup` and `Cdn` for an `Electron`: the column
`k` of the matrix of a species holds the amplitudes of its orbital `k` on the sites holding that
species, in their order, the columns being orthonormal. The state is
``\\prod_k \\left(\\sum_i \\Phi_{ik} c_i^\\dagger\\right) |0\\rangle``, up to a global
phase, so that ``\\langle c_i^\\dagger c_j \\rangle = \\sum_k \\overline{\\Phi_{ik}}
\\Phi_{jk}``. The species are taken in the order they first appear along the sites.

The sites holding no fermions, qubits, spins or bosons beside the fermions, take the local state
`others`, one for all of them or a vector of one for each, which they then need: an impurity
in a Fermi sea, for instance.

It is built after Fishman and White (2015), as a circuit of rotations of neighbouring modes
applied to a product state, exact to 1e-12 on the correlations: `limits` constrains the
truncations made as the gates are applied, the entanglement of a determinant growing with
the size of the system. The number of fermions of each species is conserved, so the sites may
conserve it, weakly or strongly.

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
                v = matrix(dag(c), system[k]) * (v isa String ? state(system[k], v) : v)
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

the constant `E0` and the matrices `h` of the terms ``\\sum_{ij} h_{ij} c_i^\\dagger c_j`` of
each species `c` of a hamiltonian of free fermions, its modes on the sites `species` gives with
it, read off the product state of no fermion, of local states `references`, see
`reference_states`, and the states of one fermion added to it: ``E_0 = \\langle 0 | H | 0
\\rangle`` and ``h_{ij} = \\langle 1_i | H | 1_j \\rangle - E_0 \\delta_{ij}``. The MPO of the
hamiltonian is contracted between product states by its transfer matrices, without a charge,
so that the basis of a site is that of its vectors. Whether the hamiltonian is of that form is
left to the caller.
"""
function one_body(system::System, hamiltonian::IndexedOp{Pure}, species, references)
    plain = weaken(system, ())
    n = length(plain)
    zero = [ v isa String ? state(plain[k], v) : v for (k, v) in enumerate(references) ]
    w = make_mpo(State{Pure}(plain, zero), hamiltonian)
    s = [ SysIndex{Pure}(plain, k) for k in 1:n ]
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

the ground state of `nparticles` free fermions of each species of the sites, as
`fermion_species` gives them in the order they first appear, under `hamiltonian`, a quadratic
hamiltonian conserving the number of each, ``E_0 + \\sum_{ij} h_{ij} c_i^\\dagger c_j`` for
each species: the Slater determinant of the `nparticles` lowest orbitals of `h`, see
`slater_state`, whose `others` and `limits` it takes.

The hamiltonian is written as for any other function, on any graph, with potentials or
fluxes. A hamiltonian of another form is refused, an interaction, a term mixing species or a
term acting on the sites holding no fermions for instance, and so is a number of fermions
that fills a degenerate level in part, whose ground state is not unique: the message gives the
nearest numbers that fill it, and a small term added to the hamiltonian lifts the degeneracy.

# Examples

    sea = fermi_sea(System(20, Fermion()), -sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i)
                                                for i in 1:19), 7)
    ring = -sum(dag(C)(i) * C(j) + dag(C)(j) * C(i) for (i, j) in circle_graph(10))
    fermi_sea(System(10, Fermion()), ring + 1e-3 * N(1), 4)    # 4 alone fills a level in part
    hop(c) = -sum(dag(c)(i) * c(i + 1) + dag(c)(i + 1) * c(i) for i in 1:7)
    fermi_sea(System(8, Electron()), hop(Cup) + hop(Cdn), 3, 3)
"""
function fermi_sea(system::System, hamiltonian::IndexedOp{Pure}, nparticles::Int...;
                   others = nothing, limits::Limits = Limits())
    species = species_sites(system)
    if length(nparticles) ≠ length(species)
        error("these sites hold $(length(species)) species of fermions, which take as many " *
              "numbers of fermions, $(join(first.(species), " and ")) in that order")
    end
    references = reference_states(system, others)
    e0, hs = one_body(system, hamiltonian, species, references)
    # the hamiltonian rebuilt from what the states of no fermion and of one see of it: the
    # same unless it holds terms they do not see, an interaction, a term mixing species or one
    # acting on the sites without fermions. Compared by the Hilbert-Schmidt norms of MPOs, which
    # `≈` does not replace: it compares how the operators are written, and tells N(1) from
    # dag(C)(1) * C(1)
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
