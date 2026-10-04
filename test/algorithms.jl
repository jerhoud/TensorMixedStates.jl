# Complete physical scenarios.
#
# Goes here: multi phase simulations reproducing a known physical result, ground state
# search, steady state, and the graph helpers used to build them. These are the slowest
# tests of the suite; keep the systems small.

@testset "Graph utilities" begin
    @test line_graph(4) == [(1, 2), (2, 3), (3, 4)]
    @test circle_graph(4) == [(1, 2), (2, 3), (3, 4), (4, 1)]
    @test circle_graph(2) == [(1, 2), (2, 1)]
    @test_throws "a ring has two vertices at least" circle_graph(1)
    @test complete_graph(4) == [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    @test graph_base_size(circle_graph(7)) == 7
    # the snake runs down the columns: 1 4 5 on the first row, 2 3 6 on the second
    @test square_lattice(3, 2) == [(1, 2), (3, 4), (5, 6), (1, 4), (2, 3), (4, 5), (3, 6)]
    @test square_lattice(3) == square_lattice(3, 3)
    for (nx, ny) in [(2, 2), (3, 2), (4, 3), (10, 3), (1, 5), (5, 1)]
        g = square_lattice(nx, ny)
        @test graph_base_size(g) == nx * ny
        @test length(g) == ny * (nx - 1) + nx * (ny - 1)
        @test allunique(g)
        # this is what the snake buys: no bond spans more than 2ny - 1 sites, whatever
        # the length of the lattice along x
        @test maximum(abs(b - a) for (a, b) in g) ≤ max(1, 2ny - 1)
    end
end

@testset "A failed Check fails test_phases" begin
    # what every test_phases of a check relies on: a Check that fails stops the run with its
    # error, rather than being logged. The tolerance bounds the norm of the difference
    @test_throws "failed with values" redirect_stdout(devnull) do
        test_phases([CreateState{Pure}(2, Qubit(), "Up"),
                     Evolve(duration = 0.2, time_step = 0.1, algo = Tdvp(), evolver = -im * Z(1),
                            final_measurements = check(Z(1), -1))])
    end
end

@testset "Complete graphs" begin
    @test_ok test_phases(create_graph_state(complete_graph(4);
        final_measurements = check([X, Y, Z, X(1)Z(2)Z(3)Z(4), (Y, Y)],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], 1, [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]])))
    @test_ok test_phases([ create_graph_state(complete_graph(4)), ToMixed(;
        final_measurements = check([X, Y, Z, X(1)Z(2)Z(3)Z(4), (Y, Y)],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], 1, [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]]))])

end

@testset "Slater determinants and Fermi seas" begin
    LA = TensorMixedStates.LinearAlgebra
    orth(n, m) = Matrix(LA.qr(randn(ComplexF64, n, m)).Q)[:, 1:m]

    # the amplitude of each configuration is the determinant of the orbitals on its sites, in
    # the order of the Jordan-Wigner convention, up to a global phase
    n, m = 6, 3
    Φ = orth(n, m)
    d = zeros(ComplexF64, 2^n)
    for b in 0:2^n-1
        occ = [ (b >> (n - k)) & 1 for k in 1:n ]
        if sum(occ) == m
            d[b + 1] = LA.det(Φ[findall(==(1), occ), :])
        end
    end
    @test abs(LA.dot(d / LA.norm(d), dense_vector(slater_state(System(n, Fermion()), Φ)))) ≈ 1

    # the correlations, on sites conserving the number of fermions, and Wick's theorem
    n, m = 8, 3
    Φ = orth(n, m)
    ψ = slater_state(System(n, Fermion(conserve = N)), Φ)
    Λ = conj(Φ) * transpose(Φ)
    @test expect2(ψ, (dag(C), C)) ≈ Λ atol = 1e-10
    @test [ expect(ψ, N(i) * N(j)) for i in 1:n for j in 1:n if i ≠ j ] ≈
          [ Λ[i, i] * Λ[j, j] - Λ[i, j] * Λ[j, i] for i in 1:n for j in 1:n if i ≠ j ] atol = 1e-10

    # electrons, a matrix of orbitals per spin, under a strong conservation
    up, dn = orth(5, 2), orth(5, 3)
    ψ = slater_state(System(5, Electron(conserve = (TensorMixedStates.strong(Ntot), 2Sz))),
                     up, dn)
    @test expect2(ψ, (dag(Cup), Cup)) ≈ conj(up) * transpose(up) atol = 1e-10
    @test expect2(ψ, (dag(Cdn), Cdn)) ≈ conj(dn) * transpose(dn) atol = 1e-10

    # the Fermi sea of a hamiltonian written in several ways, with a potential, a constant and a
    # flux through a ring: its energy is the sum of the lowest levels
    n, m = 12, 5
    flux = exp(0.4im)
    H = -sum((dag(C) ⊗ C)(i, i+1) for i in 1:n-1) + sum(C(i) * dag(C)(i+1) for i in 1:n-1) +
        0.3 * N(4) + 2 * Id(1) - flux * dag(C)(n) * C(1) - conj(flux) * dag(C)(1) * C(n)
    h = zeros(ComplexF64, n, n)
    for i in 1:n-1
        h[i, i+1] = h[i+1, i] = -1
    end
    h[4, 4] = 0.3
    h[n, 1], h[1, n] = -flux, -conj(flux)
    ψ = fermi_sea(System(n, Fermion()), H, m)
    @test real(expect(ψ, H)) ≈ 2 + sum(LA.eigvals(LA.Hermitian(h))[1:m]) atol = 1e-9
    @test sum(real(expect1(ψ, N))) ≈ m
    hop(c, k) = -sum(dag(c)(i) * c(i + 1) + dag(c)(i + 1) * c(i) for i in 1:k-1)
    he = hop(Cup, 6) + 0.5 * hop(Cdn, 6)
    ψ = fermi_sea(System(6, Electron(conserve = (Ntot, 2Sz))), he, 3, 1)
    ε = LA.eigvals([ abs(i - j) == 1 ? -1. : 0. for i in 1:6, j in 1:6 ])
    @test real(expect(ψ, he)) ≈ sum(ε[1:3]) + 0.5 * ε[1] atol = 1e-9

    # a level filled in part has no single ground state, which a small term settles
    ring = -sum(dag(C)(i) * C(j) + dag(C)(j) * C(i) for (i, j) in circle_graph(10))
    @test_throws "take 3 or 5" fermi_sea(System(10, Fermion()), ring, 4)
    @test sum(real(expect1(fermi_sea(System(10, Fermion()), ring + 1e-3 * N(1), 4), N))) ≈ 4

    # refused: what is not a hamiltonian of free fermions, and arguments that do not fit
    free = hop(C, 4)
    @test_throws "not free" fermi_sea(System(4, Fermion()), free + N(1) * N(2), 2)
    @test_throws "not free" fermi_sea(System(4, Fermion()),
                                      free + C(1) * C(2) + dag(C)(2) * dag(C)(1), 2)
    @test_throws "not free" fermi_sea(System(4, Electron()),
                                      hop(Cup, 4) + dag(Cup)(1) * Cdn(1) + dag(Cdn)(1) * Cup(1), 1, 1)
    @test_throws "take as many numbers" fermi_sea(System(4, Electron()), hop(Cup, 4), 1)
    @test_throws "not orthonormal" slater_state(System(4, Fermion()), ones(4, 1))
    @test_throws "holds no free fermions" slater_state(System(4, Qubit()), orth(4, 1))

    # sites holding no fermions take the local state given, and the strings run across them:
    # the determinant on the fermions, the qubit down
    sys = System([Fermion(), Fermion(), Qubit(), Fermion(), Fermion()])
    fs = [1, 2, 4, 5]
    Φ = orth(4, 2)
    ψ = slater_state(sys, Φ; others = "Dn")
    d = zeros(ComplexF64, 2^5)
    for b in 0:2^5-1
        occ = [ (b >> (5 - k)) & 1 for k in 1:5 ]
        if occ[3] == 1 && sum(occ[fs]) == 2
            d[b + 1] = LA.det(Φ[findall(==(1), occ[fs]), :])
        end
    end
    @test abs(LA.dot(d / LA.norm(d), dense_vector(ψ))) ≈ 1
    Λ = conj(Φ) * transpose(Φ)
    @test [ expect(ψ, dag(C)(fs[p]) * C(fs[q])) for p in 1:4, q in 1:4 ] ≈ Λ atol = 1e-10

    # an impurity in a Fermi sea, which a term acting on it or coupling it would make not free
    hf = -sum(dag(C)(i) * C(j) + dag(C)(j) * C(i) for (i, j) in [(1, 2), (2, 4), (4, 5)])
    ψ = fermi_sea(sys, hf, 2; others = "Up")
    @test real(expect(ψ, hf)) ≈ sum(LA.eigvals([ abs(i - j) == 1 ? -1. : 0. for i in 1:4, j in 1:4 ])[1:2])
    @test expect(ψ, Z(3)) ≈ 1
    @test_throws "acts on the sites holding none" fermi_sea(sys, hf + Z(3), 2; others = "Up")
    @test_throws "not free" fermi_sea(sys, hf + X(3) * N(1), 2; others = "Up")
    @test_throws "site 3 holds no fermions" fermi_sea(sys, hf, 2)

    # each species on the sites holding it, in the order they first appear
    ψ = slater_state(System([Electron(), Fermion(), Electron(), Fermion()]),
                     orth(2, 1), orth(2, 1), orth(2, 1))
    @test real(expect(ψ, Nup(1) + Nup(3))) ≈ 1
    @test real(expect(ψ, Ndn(1) + Ndn(3))) ≈ 1
    @test real(expect(ψ, N(2) + N(4))) ≈ 1
end

@testset "GHZ states" begin
    g = ghz_state(System(6, Qubit()), "Up", "Dn")
    @test maxlinkdim(g) == 2
    @test dense_vector(g) ≈ [ k == 1 || k == 64 ? 1 / sqrt(2) : 0 for k in 1:64 ]
    @test expect(g, prod(X(i) for i in 1:6)) ≈ 1
    @test real(expect(ghz_state(System(1, Qubit()), "Up", "Dn"), X(1))) ≈ 1
    @test norm(ghz_state(System(5, Qudit(3)), "0", "1", "2")) ≈ 1
    # electrons of either spin, which conserving their number allows, both product states having
    # the same, but not up and down qubits counted by N
    @test expect(ghz_state(System(4, Electron(conserve = Ntot)), "Up", "Dn"), Sz(1) * Sz(4)) ≈ 0.25
    @test_throws "no definite charge" ghz_state(System(4, Qubit(conserve = N)), "Up", "Dn")
    @test_throws "one local state at least" ghz_state(System(4, Qubit()))
end

@testset "States from the tensors of their MPS and from their dense vector" begin
    LA = TensorMixedStates.LinearAlgebra
    # the cluster state, against the graph state of a chain
    A = zeros(2, 2, 2)
    for l in 1:2, s in 1:2
        A[l, s, s] = l == s == 2 ? -1 : 1
    end
    c = mps_state(System(6, Qubit()), [ A[1:1, :, :], fill(A, 4)..., sum(A; dims = 3) ])
    @test maxlinkdim(c) == 2
    @test abs(LA.dot(dense_vector(c), dense_vector(graph_state(line_graph(6))))) ≈ 1
    @test_throws "3 tensors cannot" mps_state(System(4, Qubit()), [A[1:1, :, :], A, sum(A; dims = 3)])
    @test_throws "left dimension of 2 where 1" mps_state(System(2, Qubit()), [A, sum(A; dims = 3)])
    @test_throws "right dimension of 2" mps_state(System(2, Qubit()), [A[1:1, :, :], A])
    @test_throws "3 states for a site of dimension 2" mps_state(System(2, Qubit()),
                                                                [zeros(1, 3, 1), zeros(1, 2, 1)])
    # a pure state, and one of definite charge on fermions conserving it
    ψ = [ sin(3k) + im * cos(5k) for k in 1:32 ]
    @test dense_vector(dense_state(System(5, Qubit()), ψ)) ≈ ψ / LA.norm(ψ)
    φ = [ count_ones(b) == 2 ? sin(3b) : 0. for b in 0:15 ]
    dc = dense_state(System(4, Fermion(conserve = N)), φ)
    @test sum(expect1(dc, N)) ≈ 2
    @test dense_vector(dc) ≈ φ / LA.norm(φ)
    @test maxlinkdim(dense_state(System(8, Qubit()), kron(fill([1., 1.], 8)...))) == 1
    @test_throws "no definite charge" dense_state(System(4, Fermion(conserve = N)), ones(16))
    @test_throws "cannot be a state of 16" dense_state(System(4, Qubit()), ones(8))
    # a density matrix, against the expectations it gives
    M = [ sin(k + 2l) + im * cos(3k - l) for k in 1:8, l in 1:8 ]
    ρ = M * M' / LA.tr(M * M')
    m = dense_state(System(3, Qubit()), ρ)
    z, x, id = [1 0 ; 0 -1], [0 1 ; 1 0], [1 0 ; 0 1]
    @test expect(m, Z(1) * Z(3)) ≈ real(LA.tr(ρ * kron(z, id, z)))
    @test expect(m, X(2) * Z(3)) ≈ real(LA.tr(ρ * kron(id, x, z)))
    @test trace2(m) ≈ real(LA.tr(ρ * ρ))
    # on fermions, a coherence within a sector, under weak and strong conservations
    ρf = [ 0 0 0 0 ; 0 0.5 0.3im 0 ; 0 -0.3im 0.5 0 ; 0 0 0 0 ]
    hop = jw_matrix(Fermion(), 2, [dag(C) => 1, C => 2])
    for site in (Fermion(conserve = N), Fermion(conserve = TensorMixedStates.strong(N)))
        mf = dense_state(System(2, site), ρf)
        @test expect(mf, dag(C)(1) * C(2)) ≈ LA.tr(ρf * hop)
    end
    @test_throws "no definite charge" dense_state(System(2, Fermion(conserve = N)), ones(4, 4))
    @test_throws "size (4, 4) cannot" dense_state(System(3, Qubit()), ones(4, 4))
end

@testset "Superpositions and mixtures of product states" begin
    up, dn = [1., 0.], [0., 1.]
    s = superposition(System(4, Qubit()), [1 => ["Up", "Dn", "Up", "Dn"],
                                           -1 => ["Dn", "Up", "Dn", "Up"], 0.5im => "+"])
    ref = kron(up, dn, up, dn) - kron(dn, up, dn, up) + 0.5im * kron(fill([1., 1.] / sqrt(2), 4)...)
    @test abs(TensorMixedStates.LinearAlgebra.dot(dense_vector(s), ref / norm(ref))) ≈ 1
    @test maxlinkdim(s) == 3
    # on charged sites, the terms of the same charge
    sc = superposition(System(4, Fermion(conserve = N)),
                       [1 => ["Occ", "Emp", "Occ", "Emp"], 2 => ["Emp", "Occ", "Emp", "Occ"]])
    @test expect(sc, N(1)) ≈ 0.2
    @test_throws "no definite charge" superposition(System(4, Fermion(conserve = N)),
                                                    [1 => "Occ", 1 => "Emp"])
    # a term of zero coefficient has no channel and no say in the charge
    @test maxlinkdim(superposition(System(4, Qubit()), [1 => "Up", 0 => "Dn"])) == 1
    @test expect(superposition(System(4, Fermion(conserve = N)), [1 => "Occ", 0 => "Emp"]), N(1)) ≈ 1
    @test_throws "a nonzero term" superposition(System(4, Qubit()), [0 => "Up", 0 => "Dn"])
    @test_throws "a nonzero term" mixture(System(4, Qubit()), [0 => "Up"])
    @test_throws "cannot be one of 4 sites" superposition(System(4, Qubit()), [1 => "Up", 0 => ["Dn"]])
    @test maxlinkdim(superposition(System(8, Qubit()), [ k => isodd(k) ? "Up" : "Dn" for k in 1:10 ];
                                   limits = Limits(maxdim = 1))) == 1
    # a classical ensemble of two Néel states, and one of mixed local states
    m = mixture(System(4, Qubit()), [0.5 => ["Up", "Dn", "Up", "Dn"], 0.5 => ["Dn", "Up", "Dn", "Up"]])
    @test trace2(m) ≈ 0.5
    @test expect(m, Z(1) * Z(2)) ≈ -1
    mm = mixture(System(2, Qubit()), [1 => ["FullyMixed", "Up"], 3 => [[0.5 0.5 ; 0.5 0.5], "Dn"]])
    @test expect(mm, Z(2)) ≈ -0.5
    @test expect(mm, X(1)) ≈ 0.75
    @test_throws "real and not negative" mixture(System(2, Qubit()), [-1 => "Up"])
end

@testset "Dicke and W states" begin
    d = dicke_state(System(5, Qubit()), 2, "Up", "Dn")
    @test maxlinkdim(d) == 3
    @test dense_vector(d) ≈ [ count_ones(b) == 2 ? 1 / sqrt(10) : 0. for b in 0:31 ]
    # every product state of a Dicke state has the same charge
    dc = dicke_state(System(8, Qubit(conserve = N)), 4, "Up", "Dn")
    @test sum(expect1(dc, N)) ≈ 4
    @test expect1(w_state(System(6, Qubit()), "Up", "Dn"), Z) ≈ fill(2 / 3, 6)
    @test expect1(dicke_state(System(3, Qubit()), 0, "Up", "Dn"), Z) ≈ [1, 1, 1]
    @test expect1(dicke_state(System(3, Qubit()), 3, "Up", "Dn"), Z) ≈ [-1, -1, -1]
    @test_throws "do not fit on 3 sites" dicke_state(System(3, Qubit()), 4, "Up", "Dn")
    LA = TensorMixedStates.LinearAlgebra
    # a wave packet of one excitation, a site of zero amplitude included
    up, dn = [1., 0.], [0., 1.]
    φ = [1, 2im, 0, -1]
    w = w_state(System(4, Qubit()), "Up", "Dn", φ)
    ref = sum(φ[j] * kron([ k == j ? dn : up for k in 1:4 ]...) for j in 1:4)
    @test dense_vector(w) ≈ ref / LA.norm(ref)
    @test maxlinkdim(w) == 2
    @test sum(expect1(w_state(System(6, Qubit(conserve = N)), "Up", "Dn", 1:6), Z)) ≈ 4
    @test_throws "a nonzero amplitude" w_state(System(3, Qubit()), "Up", "Dn", zeros(3))
    @test_throws "2 amplitudes for 3 sites" w_state(System(3, Qubit()), "Up", "Dn", [1, 1])
end

@testset "Dimer states" begin
    # pairs that cross, against the product of the two singlets, |ab> - |ba> with a first
    up, dn = [1., 0.], [0., 1.]
    ref = sum(c1 * c2 * kron(s1, s2, s3, s4) for (s1, s3, c1) in ((up, dn, 1), (dn, up, -1))
                                              for (s2, s4, c2) in ((up, dn, 1), (dn, up, -1))) / 2
    @test dense_vector(dimer_state(System(4, Qubit()), [(1, 3), (2, 4)], "Up", "Dn")) ≈ ref
    # Majumdar-Ghosh, on charged sites
    mg = dimer_state(System(10, Qubit(conserve = N)), [ (i, i + 1) for i in 1:2:9 ], "Up", "Dn")
    @test maxlinkdim(mg) == 2
    @test expect(mg, Z(1) * Z(2)) ≈ -1
    @test expect(mg, Z(2) * Z(3)) ≈ 0 atol = 1e-12
    # the rainbow state, nested pairs
    rb = dimer_state(System(8, Qubit()), [ (i, 9 - i) for i in 1:4 ], "Up", "Dn")
    @test maxlinkdim(rb) == 16
    @test expect(rb, X(1) * X(8)) ≈ -1
    # a site in no pair, and electrons, whose spins make the singlet
    sp = dimer_state(System(5, Spin(1/2)), [(1, 2), (4, 5)], "1/2", "-1/2"; others = "1/2")
    @test expect(sp, Sx(1) * Sx(2) + Sy(1) * Sy(2) + Sz(1) * Sz(2)) ≈ -0.75
    @test expect(sp, Sz(3)) ≈ 0.5
    el = dimer_state(System(4, Electron(conserve = (Ntot, 2Sz))), [(1, 4), (2, 3)], "Up", "Dn")
    @test expect(el, Sz(1) * Sz(4) + (Sp(1) * Sm(4) + Sm(1) * Sp(4)) / 2) ≈ -0.75
    @test_throws "a site is in two pairs" dimer_state(System(4, Qubit()), [(1, 2), (2, 3)],
                                                      "Up", "Dn")
    @test_throws "reach beyond the 4 sites" dimer_state(System(4, Qubit()), [(1, 5)], "Up", "Dn";
                                                        others = "Up")
    @test_throws "site 3 is in no pair" dimer_state(System(3, Qubit()), [(1, 2)], "Up", "Dn")
    @test_throws "needs them orthonormal" dimer_state(System(2, Qubit()), [(1, 2)], "Up", "+")
end

@testset "AKLT states" begin
    # an exact ground state whatever the virtual spins at the ends, which set the total Sz
    n = 10
    ss(i) = Sz(i) * Sz(i+1) + (Sp(i) * Sm(i+1) + Sm(i) * Sp(i+1)) / 2
    h = sum(ss(i) + ss(i) * ss(i) / 3 for i in 1:n-1)
    for site in (Spin(1), Spin(1, conserve = Sz)), (l, r, sz) in (("Up", "Up", 0), ("Up", "Dn", 1))
        a = aklt_state(System(n, site); left = l, right = r)
        @test maxlinkdim(a) == 2
        @test expect(a, h) ≈ -2 * (n - 1) / 3
        @test variance(a, h) ≈ 0 atol = 1e-10
        @test sum(expect1(a, Sz)) ≈ sz atol = 1e-12
    end
    @test_throws "a chain of spins one" aklt_state(System(3, Qubit()))
end

@testset "Fully mixed states of a sector" begin
    LA = TensorMixedStates.LinearAlgebra
    strong = TensorMixedStates.strong
    # the projector on the sector over its dimension, whatever the sites conserve
    n, m = 6, 2
    for site in (Fermion(), Fermion(conserve = N), Fermion(conserve = strong(N)))
        ρ = fully_mixed(System(n, site), N => m)
        @test trace(ρ) ≈ 1
        @test trace2(ρ) ≈ 1 / binomial(n, m)
        @test expect(ρ, N(2)) ≈ m / n
        @test expect(ρ, N(1) * N(4)) ≈ m * (m - 1) / (n * (n - 1))
        @test maxlinkdim(ρ) == m + 1
    end
    # two quantities, the second a half integer, under a strong conservation: 36 states of four
    # electrons on four sites with no magnetization
    ρ = fully_mixed(System(4, Electron(conserve = (strong(Ntot), 2Sz))), Ntot => 4, Sz => 0)
    @test trace2(ρ) ≈ 1 / 36
    @test sum(expect1(ρ, Sz)) ≈ 0 atol = 1e-12

    # the canonical thermal state of fermions hopping, which Thermalize reaches from it
    n, m, β = 4, 2, 0.9
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:n-1) + 0.4 * N(1)
    c(j) = foldl(kron, [ k < j ? [1. 0. ; 0. -1.] : k == j ? [0. 1. ; 0. 0.] : [1. 0. ; 0. 1.]
                         for k in 1:n ])
    dense = -sum(c(i)' * c(i+1) + c(i+1)' * c(i) for i in 1:n-1) + 0.4 * c(1)' * c(1)
    sector = LA.Diagonal(Float64.(round.(LA.diag(sum(c(i)' * c(i) for i in 1:n))) .== m))
    g = exp(-β * dense) * sector
    for site in (Fermion(conserve = N), Fermion(conserve = strong(N)))
        _, ρ = thermal_state(h, β, fully_mixed(System(n, site), N => m); nsteps = 18,
                             limits = Limits(cutoff = 1e-14, maxdim = 64))
        @test expect(ρ, h) ≈ LA.tr(g * dense) / LA.tr(g) atol = 1e-5
    end

    @test_throws "no basis state of the system has N = 5" fully_mixed(System(4, Fermion()), N => 5)
    @test_throws "not diagonal" fully_mixed(System(4, Qubit()), X => 2)
    @test_throws "does not fix the charges" fully_mixed(
        System(4, Electron(conserve = (strong(Ntot), 2Sz))), Sz => 0)
end

@testset "Dmrg" begin
    @test_ok test_phases([
        CreateState(
            type = Pure(),
            system = System(5, Qubit()),
            randomize = 10,
            ),
            GroundState(
                hamiltonian = sum(-Z(i) for i in 1:5),
                limits = Limits(maxdim = 10),
                nsweeps = 2,
                final_measurements = check([X, Y, Z, Norm], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1], 1], 1e-7)
                )
                ])
            end
            @testset "GHZ" begin
                @test_ok begin
                    sys = System(6, Qubit())
                    # adding mixed states adds density matrices, so summing the two mixed
                    # product states gives their classical mixture and not a GHZ state.
                    # The superposition has to be made in pure representation, where the
                    # addition is on amplitudes.
                    ghz = mix((State{Pure}(sys, "Up") + State{Pure}(sys, "Dn")) / sqrt(2))
                    test_phases([
        CreateState(type = Mixed(), state = ghz,
            final_measurements = check([Purity, prod(X(i) for i in 1:6)], [1, 1], 1e-10)),
        PartialTrace(
            keep_positions = [2, 3, 5],
            final_measurements = check([X, Y, Z, (Z, Z)], [[0, 0, 0], [0, 0, 0], [0, 0, 0], [1 1 1 ; 1 1 1 ; 1 1 1]])
            )
            ])
        end
    end

@testset "Ising chain" begin
    # The reference values are exact, computed by diagonalizing the 64 dimensional Hilbert
    # space in test/reference/ising_ed.jl, which shares no code with what is tested here.
    # The tolerance is therefore the error of the evolution below, about 2e-8, and not the
    # precision of the references. The last one is analytic rather than computed: the
    # product of all X commutes with the hamiltonian and starts at 1, so it stays at 1.
    @test_ok test_phases([
        CreateState{Pure}(6, Qubit(), "X+"),
        Evolve(
    algo =  ApproxW(order = 4, w = 2),
    limits = Limits(maxdim = 8),
    duration = 1.0,
    time_step = 0.02,
    evolver =
        -im*(sum(Z(i)*Z(i+1) for i in 1:5)+Z(6)*Z(1)-sum(X(i) for i in 1:6)),
    final_measurements = [
        check([X,Y,Z],[[0.48881258418,0.48881258418,0.48881258418,0.48881258418,0.48881258418,0.48881258418],[0.0,0,0,0,0,0],[0.0,0,0,0,0,0]],1e-7),
        check([Z(1)Z(2),Z(2)Z(3),Z(1)Z(6)],[-0.51118741582,-0.51118741582,-0.51118741582],1e-7),
        check([Y(1)Y(2),Y(2)Y(3),Y(1)Y(6)],[-0.2518341076,-0.2518341076,-0.2518341076],1e-7),
        check([X(1)X(2),X(2)X(3),X(1)X(6)],[0.1123126174,0.1123126174,0.1123126174],1e-7),
        check([X(1)X(2)X(3)X(4),X(2)X(3)X(4)X(5),X(4)X(5)X(6)X(1)],[0.1123126174,0.1123126174,0.1123126174],1e-7),
        check(EntanglementEntropy(3), 1.15220908566, 1e-7),
        check(X(1)X(2)X(3)X(4)X(5)X(6),1.0,1e-8)
        ])
        ])
    end

@testset "An expansion before the first step" begin
    # the ring above from a product state, with tdvp: the first step, of bond dimension one,
    # left the tangent space through Z(6)Z(1), and the expansion, coming after it, left an
    # error of order the time step, 0.016 on X(1), whatever expand_period
    for expand_period in (1, 2)
        @test_ok test_phases([
            CreateState{Pure}(6, Qubit(), "X+"),
            Evolve(algo = Tdvp(; expand_period), limits = Limits(maxdim = 8, cutoff = 1e-14),
                   duration = 1.0, time_step = 0.1,
                   evolver = -im * (sum(Z(i) * Z(i + 1) for i in 1:5) + Z(6) * Z(1) -
                                    sum(X(i) for i in 1:6)),
                   final_measurements = [check(X(1), 0.48881258418, 1e-8),
                                     check(Z(1)Z(6), -0.51118741582, 1e-8)])])
    end
end

@testset "Free fermions with source" begin
    # The reference values are exact, from the dense Lindblad evolution of the 32
    # dimensional Fock space in test/reference/fermion_lindblad.jl, which shares no code
    # with what is tested here. The tolerance is the error of the evolution below, about
    # 8e-8.
    @test_ok test_phases([
        CreateState{Mixed}(5, Fermion(), "0"),
        Evolve(
            algo=Tdvp(),
            limits = Limits(maxdim = 16),
            duration = 1.0,
            time_step = 0.05,
            evolver =
                -im * sum(dag(C)(i)*C(i+1)+dag(C)(i+1)*C(i) for i in 1:4) + Dissipator(sqrt(2*0.2)*dag(C))(3),
            final_measurements = [
                check(N, [0.0125582080327,0.063008590052,0.1950187333854,0.063008590052,0.0125582080327], 1e-6),
                check([dag(C)(3)*C(i) for i in 1:5],
                [-0.023529279887, -0.0783682258151im,0.1950187333854,-0.0783682258151im,-0.023529279887],1e-6),
                check(Purity,0.5256648672193,1e-6)
            ])
    ])
end

@testset "Free bosons with source" begin
    # The reference values are exact, from the Gaussian moments of this quadratic Lindbladian
    # in test/reference/boson_gaussian.jl, which shares no code with what is tested here: with
    # four sites of dimension 7 a dense reference, on a vectorized Liouvillian of 2401^2, is
    # out of reach. The evolution below is 6e-7 to 8e-7 away from them.
    #
    # The bond dimension is 16 rather than the 10 first used, and that is not a detail. At
    # 10 the truncation itself is unstable: which singular values survive depends on the
    # rounding, so the answer follows the order the sums happen to be taken in. Changing
    # nothing but the BLAS thread count on one machine moved this correlation by 4.4e-4
    # relative, an overall deviation swinging between 8.8e-7 and 8.0e-6, and the Windows
    # job, landing on a different thread count from the Linux and macOS ones, missed the
    # 1e-5 tolerance by one percent. At 16 the dependence is gone, the deviation being
    # 8.6e-7 whatever the thread count, so the tolerance measures the accuracy of the
    # method again rather than the arithmetic of the runner. It is left at 1e-5 and not
    # tightened, since 8.6e-7 leaves too little margin under 1e-6.
    @test_ok test_phases([
        CreateState{Mixed}(4, Boson(7), "0"),
        Evolve(
            algo = ApproxW(order = 4, w = 2),
            limits = Limits(maxdim = 16),
            duration = 0.3,
            time_step = 0.1,
            evolver =
                -im*sum(A(i)*dag(A)(i+1)+dag(A)(i)*A(i+1) for i in 1:3) + Dissipator(2*sqrt(0.1)*dag(A))(2),
            final_measurements = [
                    check(N, [0.0036328913502, 0.1200236303285, 0.0035680174576, 4.86613106e-5], 1e-5),
                    check([dag(A)(2)*A(i) for i in 1:4],
                          [-0.0180013634952im, 0.1200236303285, -0.0178675798942im, -0.0017865872242], 1e-5),
                    check(Purity, 0.7965328475916, 1e-5)
                ])
    ])
end

@testset "Steady state" begin
    @test_ok test_phases([
        CreateState{Mixed}(2, Qubit(), "+"),
        SteadyState(
            lindbladian = Dissipator(Sp)(1) + Dissipator(Sm)(2),
            nsweeps = 20,
            limits = Limits(cutoff = 1e-10, maxdim = 10),
            final_measurements = check([X, Y, Z], [[0, 0], [0, 0], [1, -1]], 1e-10)
        )
    ])
    @test_ok test_phases([
        CreateState(
            type = Mixed(),
            system = System(5, Qubit()),
            randomize = 10,
        ),
        # the tolerances of these three searches hid an error of the default Krylov search of
        # ITensorMPS, 1e-5 on the first, which the default of steady_state now resolves
        SteadyState(
            lindbladian = -im * (-sum(Z(i)Z(i+1) for i in 1:4)) + sum(Dissipator(Sp)(i) for i in 1:5),
            limits = Limits(maxdim = 10, cutoff = 1e-10),
            nsweeps = 40,
            final_measurements = check([X, Y, Z], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1]], 1e-8)
        )
    ])
    @test_ok test_phases([
        CreateState{Mixed}(4, Qubit(), "FullyMixed"),
        SteadyState(
            lindbladian =
                -im * (sum(X(i)X(i+1)+Y(i)Y(i+1) for i in 1:3))
                + Dissipator(Sp)(1) + Dissipator(Sm)(4),
            nsweeps = 200,
            limits = Limits(cutoff = 1e-10, maxdim = 10),
            final_measurements = [
                check(Z, [0.05882352941176472, 0.0, 0.0,-0.05882352941176472], 1e-7),
                check([2(X(i)Y(i+1)-Y(i)X(i+1)) for i in 1:3], fill(0.9411764705882353, 3), 1e-7),
            ]
        )
    ])
    # the state found is a density matrix, of trace one, where the eigenvector of (L+)L has
    # norm one and a sign of its own
    ρ0 = mix(State{Pure}(System(2, Qubit()), "+"))
    _, ρ = steady_state(Dissipator(Sp)(1) + Dissipator(Sm)(2), ρ0;
                        nsweeps = 10, limits = Limits(cutoff = 1e-10, maxdim = 10))
    @test trace(ρ) ≈ 1
end

@testset "The options of steady_state" begin
    # it resumes at `first_sweep` as dmrg does, which is how a SteadyState phase continues
    # after a checkpoint, and takes a noise as GroundState does
    L = Dissipator(Sp)(1) + Dissipator(Sm)(2)
    ρ0 = State{Mixed}(System(2, Qubit()), "FullyMixed")
    lim = Limits(maxdim = 16)
    for first_sweep in 1:3
        obs = TensorMixedStates.ITensorMPS.DMRGObserver()
        steady_state(L, ρ0; nsweeps = 3, first_sweep, limits = lim, observer! = obs)
        @test length(obs.energies) == 4 - first_sweep
    end
    _, ρ = steady_state(L, ρ0; nsweeps = 2, limits = lim, noise = 1e-6)
    @test expect1(ρ, Z) ≈ [1, -1]
end

@testset "dmrg with no sweep left" begin
    # a search resumed after its last sweep is asked for none, on which ITensorMPS gave an
    # energy of 0: it is that of the state given
    h = -sum(Z(i) for i in 1:2) - 0.5 * X(1)
    st = State{Pure}(System(2, Qubit()), "Up")
    e, s = dmrg(h, st; nsweeps = 2, first_sweep = 3)
    @test e ≈ real(expect(st, h))
    @test s === st
end

@testset "The Krylov parameters reach dmrg" begin
    # a Krylov space of a single vector holds nothing but the state it starts from, so that
    # the search does not leave it: dmrg ends where it started
    n = 4
    h = sum(X(i) * X(i + 1) for i in 1:n - 1)
    up = State{Pure}(System(n, Qubit()), "Up")
    lim = Limits(maxdim = 16)
    @test first(dmrg(h, up; nsweeps = 3, limits = lim)) ≈ -3
    e, st = dmrg(h, up; nsweeps = 3, limits = lim, krylov = Krylov(dim = 1))
    @test e ≈ 0 atol = 1e-12
    @test expect1(st, Z) ≈ ones(n)
    # and steady_state, from the fully mixed state
    L = Dissipator(Sp)(1) + Dissipator(Sm)(2)
    ρ0 = State{Mixed}(System(2, Qubit()), "FullyMixed")
    _, ρ = steady_state(L, ρ0; nsweeps = 2, limits = lim, krylov = Krylov(dim = 1))
    @test expect1(ρ, Z) ≈ [0, 0] atol = 1e-12
    # the phases hand them on: the state stays up through both searches, where the ground
    # state of `h` has no magnetization and the steady state of `L` a spin down on site 2
    sim = runTMS(SimData(phases = [
            CreateState(type = Pure(), state = up),
            GroundState(hamiltonian = h, limits = lim, nsweeps = 3, krylov = Krylov(dim = 1)),
            ToMixed(),
            SteadyState(lindbladian = L, limits = lim, nsweeps = 2, krylov = Krylov(dim = 1))]);
        output = devnull)
    @test expect1(sim.state, Z) ≈ ones(n)
end

@testset "dmrg and steady_state refuse the options of ITensorMPS" begin
    # they no longer pass on what they do not know, `outputlevel` included
    up = State{Pure}(System(2, Qubit()), "Up")
    ρ0 = State{Mixed}(System(2, Qubit()), "FullyMixed")
    @test_throws MethodError dmrg(X(1) * X(2), up; outputlevel = 1)
    @test_throws MethodError steady_state(Dissipator(Sm)(1), ρ0; outputlevel = 1)
end

@testset "Dmrg of a hamiltonian on a mixed state" begin
    # it would minimise ρ ↦ Hρ + ρH, whose lowest eigenvector is neither the ground state nor a
    # density matrix, and is refused. A superoperator given as such is left to the caller
    h = -Z(1) * Z(2) - 0.5 * (X(1) + X(2))
    ρ = mix(RandomState{Pure}(System(2, Qubit()), 2))
    lim = Limits(maxdim = 4)
    @test_throws "ground state of a pure state" dmrg(h, ρ; nsweeps = 2, limits = lim)
    @test_throws "ground state of a pure state" runTMS(SimData(phases = [
        CreateState{Mixed}(2, Qubit(), "+"),
        GroundState(hamiltonian = h, nsweeps = 2, limits = lim)]); output = devnull)
    @test_ok dmrg(sum(Left(Z)(i) + Right(Z)(i) for i in 1:2), ρ; nsweeps = 2, limits = lim)
end

@testset "Ground and steady states on a charged system" begin
    strong = TensorMixedStates.strong
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:3)
    lind = -im * h + sum(Dissipator(sqrt(0.3) * N)(i) for i in 1:4)
    lim = Limits(cutoff = 1e-12, maxdim = 32)
    start(site) = begin
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        return (p(["Occ", "Emp", "Occ", "Emp"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
    end

    # dmrg searches inside the sector its starting state lives in, and steady_state builds
    # the MPO of `(L+)L`, whose flux is zero whenever that of `L` is. Both must land where
    # the dense computation does
    e = first(dmrg(h, start(Fermion()); nsweeps = 3, limits = lim))
    # the two modes of negative energy filled, -2cos(π/5) - 2cos(2π/5)
    @test e ≈ -sqrt(5)
    z = real(first(steady_state(lind, mix(start(Fermion())); nsweeps = 2, limits = lim)))
    # an operator prepared once, as the solvers take it, for the system it was prepared on
    s0 = start(Fermion())
    @test first(dmrg(PreMPO(s0, h), s0; nsweeps = 3, limits = lim)) ≈ e
    ρ0 = mix(s0)
    @test real(first(steady_state(PreMPO(ρ0, lind), ρ0; nsweeps = 2, limits = lim))) ≈ z atol = 1e-10
    @test_throws "prepared on another System" dmrg(PreMPO(start(Fermion()), h), s0; nsweeps = 1)
    for site in (Fermion(conserve = N), Fermion(conserve = strong(N)))
        @test first(dmrg(h, start(site); nsweeps = 3, limits = lim)) ≈ e
        @test real(first(steady_state(lind, mix(start(site)); nsweeps = 2, limits = lim))) ≈ z atol = 1e-10
    end
end

@testset "Weakening between two phases" begin
    strong = TensorMixedStates.strong
    # dephasing commutes with the number of particles and allows a strong symmetry, loss does
    # not. Weakening in between must give what a weak symmetry gives from the start
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:3)
    dephasing = -im * h + sum(Dissipator(sqrt(0.3) * N)(i) for i in 1:4)
    loss = -im * h + Dissipator(sqrt(0.2) * C)(2)
    evolve(ev) = Evolve(algo = Tdvp(), limits = Limits(cutoff = 1e-12, maxdim = 32),
                        duration = 0.5, time_step = 0.05, evolver = ev)
    run(site, middle) = runTMS(SimData(phases = [
            CreateState{Mixed}(4, site, ["Occ", "Emp", "Occ", "Emp"]),
            evolve(dephasing), middle..., evolve(loss) ]); output = devnull).state
    numbers(s) = [ real(trace(s)); real.(expect1(s, N)); expect(s, dag(C)(1) * C(3)) ]
    reference = numbers(run(Fermion(conserve = N), []))

    s = run(Fermion(conserve = strong(N)), [Weaken()])
    @test repr(symmetries(s.system)) == "N"
    @test numbers(s) ≈ reference
    s = run(Fermion(conserve = strong(N)), [Weaken(target = ())])
    @test !TensorMixedStates.is_charged(s.system)
    @test numbers(s) ≈ reference

    @test_throws "drop `strong`" run(Fermion(conserve = strong(N)), [])
end
