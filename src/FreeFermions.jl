# The states of free fermions: the Slater determinant of given orbitals, built by a circuit of
# rotations of neighbouring modes applied to a product state, and the Fermi sea of a quadratic
# hamiltonian, the determinant of its lowest orbitals.

export slater_state, fermi_sea

"""
    fermion_species(site)

the annihilation operators of the species of fermions a site holds, a mode of each: `(C,)` for a
`Fermion`, `(Cup, Cdn)` for an `Electron`. `slater_state` and `fermi_sea` read it off the sites
they are given, and a site type of one's own holding free fermions gets them by a method of it.
"""
fermion_species(site::AbstractSite) =
    error("$site holds no free fermions, see TensorMixedStates.fermion_species")

"""
    system_species(system)

the species of fermions of every site of `system`, see `fermion_species`, refused when the
sites do not all hold the same ones
"""
function system_species(system::System)
    species = fermion_species(system[1])
    for k in 2:length(system)
        if fermion_species(system[k]) ≠ species
            error("the sites of the system do not all hold the same species of fermions")
        end
    end
    return species
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
    mode_gate(c, W, j)

the gate on the sites `j` and `j + 1` that rotates the modes of annihilation operator `c` by
`W`, a unitary of 2×2: ``\\mathcal{G}\\, c_p^\\dagger \\mathcal{G}^\\dagger = \\sum_q W_{qp}
c_q^\\dagger``. It is the exponential of ``\\sum_{pq} K_{pq} c_p^\\dagger c_q`` for
``K = \\log W``, an even operator of two sites, which `apply` takes.
"""
function mode_gate(c::SimpleOp, W::AbstractMatrix, j::Int)
    K = log(W)
    # c†_{j+1} c_j placed on (j, j+1) is -c_j c†_{j+1}, the tensor product being ordered by sites
    return exp(K[1, 1] * ((dag(c) * c) ⊗ Id) + K[2, 2] * (Id ⊗ (dag(c) * c)) +
               K[1, 2] * (dag(c) ⊗ c) - K[2, 1] * (c ⊗ dag(c)))(j, j + 1)
end

"""
    slater_state(system, orbitals...; limits = Limits())

the Slater determinant of `orbitals`, one matrix for each species of fermions of the sites, as
`fermion_species` gives them, `C` for a `Fermion`, `Cup` and `Cdn` for an `Electron`: the column
`k` of the matrix of a species holds the amplitudes on the sites of its orbital `k`, the
columns being orthonormal. The state is
``\\prod_k \\left(\\sum_i \\Phi_{ik} c_i^\\dagger\\right) |0\\rangle``, up to a global phase,
so that ``\\langle c_i^\\dagger c_j \\rangle = \\sum_k \\overline{\\Phi_{ik}} \\Phi_{jk}``.

It is built after Fishman and White (2015), as a circuit of rotations of neighbouring modes
applied to a product state, exact to 1e-12 on the correlations: `limits` constrains the
truncations made as the gates are applied, the entanglement of a determinant growing with
the size of the system. The number of fermions of each species is conserved, so the sites may
conserve it, weakly or strongly.

# Examples

    slater_state(System(10, Fermion()), orbitals)
    slater_state(System(10, Electron(conserve = (Ntot, 2Sz))), up, down;
                 limits = Limits(cutoff = 1e-12))
"""
function slater_state(system::System, orbitals::AbstractMatrix...; limits::Limits = Limits())
    species = system_species(system)
    n = length(system)
    if length(orbitals) ≠ length(species)
        error("sites holding $(length(species)) species of fermions take as many matrices of " *
              "orbitals, $(join(species, " and ")) in that order")
    end
    for (c, o) in zip(species, orbitals)
        if size(o, 1) ≠ n
            error("the orbitals of $c have $(size(o, 1)) amplitudes, for $n sites")
        end
        if norm(o' * o - I) > 1e-10
            error("the orbitals of $c are not orthonormal")
        end
    end
    circuits = [ givens_circuit(o) for o in orbitals ]
    vectors = map(1:n) do k
        v = state(system[k], "0")
        for (c, (_, occupations)) in zip(species, circuits)
            if occupations[k] == 1
                v = matrix(dag(c), system[k]) * v
            end
        end
        v
    end
    st = State{Pure}(system, vectors)
    # the determinant is the product state rotated back: the rotations found last act first
    gates = [ mode_gate(c, W', j) for (c, (rotations, _)) in zip(species, circuits)
                                  for (j, W) in rotations ]
    return isempty(gates) ? st : apply(prod(gates), st; limits)
end

"""
    one_body(system, hamiltonian, species)

the constant `E0` and the matrices `h` of the terms ``\\sum_{ij} h_{ij} c_i^\\dagger c_j`` of
each species `c` of a hamiltonian of free fermions, read off the states of no fermion and of
one: ``E_0 = \\langle 0 | H | 0 \\rangle`` and ``h_{ij} = \\langle 1_i | H | 1_j \\rangle -
E_0 \\delta_{ij}``. The MPO of the hamiltonian is contracted between product states by its
transfer matrices, without a charge, so that the basis of a site is that of its vectors.
Whether the hamiltonian is of that form is left to the caller.
"""
function one_body(system::System, hamiltonian::IndexedOp{Pure}, species)
    plain = weaken(system, ())
    n = length(plain)
    vacuum = State{Pure}(plain, "0")
    w = make_mpo(vacuum, hamiltonian)
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
    zero = [ state(plain[k], "0") for k in 1:n ]
    empty = [ transfer(k, zero[k], zero[k]) for k in 1:n ]
    # the products of the empty sites on the left of a site and on its right
    left = accumulate(*, empty; init = ones(1, 1))
    right = reverse(accumulate((a, b) -> b * a, reverse(empty); init = ones(1, 1)))
    e0 = only(last(left))
    on_left(k) = k == 1 ? ones(1, 1) : left[k-1]
    on_right(k) = k == n ? ones(1, 1) : right[k+1]
    hs = map(species) do c
        one = [ matrix(dag(c), plain[k]) * zero[k] for k in 1:n ]
        h = zeros(ComplexF64, n, n)
        for i in 1:n
            h[i, i] = only(on_left(i) * transfer(i, one[i], one[i]) * on_right(i)) - e0
            # the fermion on the bra side at i and on the ket side at j > i, then the other way
            up = on_left(i) * transfer(i, one[i], zero[i])
            down = on_left(i) * transfer(i, zero[i], one[i])
            for j in i+1:n
                h[i, j] = only(up * transfer(j, zero[j], one[j]) * on_right(j))
                h[j, i] = only(down * transfer(j, one[j], zero[j]) * on_right(j))
                up *= empty[j]
                down *= empty[j]
            end
        end
        h
    end
    return e0, hs
end

"""
    fermi_sea(system, hamiltonian, nparticles...; limits = Limits())

the ground state of `nparticles` free fermions of each species of the sites, as
`fermion_species` gives them, under `hamiltonian`, a quadratic hamiltonian conserving the number
of each, ``E_0 + \\sum_{ij} h_{ij} c_i^\\dagger c_j`` for each species: the Slater determinant of
the `nparticles` lowest orbitals of `h`, see `slater_state`, whose `limits` it takes.

The hamiltonian is written as for any other function, on any graph, with potentials or
fluxes. A hamiltonian of another form is refused, an interaction or a term mixing species for
instance, and so is a number of fermions that fills a degenerate level in part, whose ground
state is not unique: the message gives the nearest numbers that fill it, and a small term
added to the hamiltonian lifts the degeneracy.

# Examples

    sea = fermi_sea(System(20, Fermion()), -sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i)
                                                for i in 1:19), 7)
    ring = -sum(dag(C)(i) * C(j) + dag(C)(j) * C(i) for (i, j) in circle_graph(10))
    fermi_sea(System(10, Fermion()), ring + 1e-3 * N(1), 4)     # 4 fills a level of the ring in part
    hop(c) = -sum(dag(c)(i) * c(i + 1) + dag(c)(i + 1) * c(i) for i in 1:7)
    fermi_sea(System(8, Electron()), hop(Cup) + hop(Cdn), 3, 3)
"""
function fermi_sea(system::System, hamiltonian::IndexedOp{Pure}, nparticles::Int...;
                   limits::Limits = Limits())
    species = system_species(system)
    n = length(system)
    if length(nparticles) ≠ length(species)
        error("sites holding $(length(species)) species of fermions take as many numbers of " *
              "fermions, $(join(species, " and ")) in that order")
    end
    e0, hs = one_body(system, hamiltonian, species)
    # the hamiltonian rebuilt from what the states of no fermion and of one see of it: the
    # same unless it holds terms they do not see, an interaction or a term mixing species.
    # Compared by the Hilbert-Schmidt norms of MPOs, which `≈` does not replace: it compares
    # how the operators are written, and tells N(1) from dag(C)(1) * C(1)
    rebuilt = e0 * Id(1) + sum(h[i, j] * dag(c)(i) * c(j) for (c, h) in zip(species, hs)
                               for i in 1:n, j in 1:n if !iszero(h[i, j]); init = 0 * Id(1))
    vacuum = State{Pure}(weaken(system, ()), "0")
    if norm(make_mpo(vacuum, hamiltonian - rebuilt)) >
       1e-10 * norm(make_mpo(vacuum, hamiltonian))
        error("the fermions of $(join(species, " and ")) are not free under this hamiltonian: " *
              "it is not quadratic, or does not conserve the number of fermions of each species")
    end
    orbitals = map(zip(species, hs, nparticles)) do (c, h, m)
        if !(0 ≤ m ≤ n)
            error("$m fermions of $c do not fit on $n sites")
        end
        if norm(h - h') > 1e-10 * max(1, norm(h))
            error("the hamiltonian is not hermitian")
        end
        e = eigen(Hermitian(h))
        ε = e.values
        tol = 1e-10 * max(1, ε[end] - ε[1])
        closed(k) = k == 0 || k == n || ε[k+1] - ε[k] > tol
        if !closed(m)
            below = findlast(closed, 0:m) - 1
            above = m + findfirst(closed, m:n) - 1
            error("$m fermions of $c fill a degenerate level in part, which has no single " *
                  "ground state: take $below or $above, or lift the degeneracy by a small term")
        end
        e.vectors[:, 1:m]
    end
    return slater_state(system, orbitals...; limits)
end
