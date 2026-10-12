# States prepared by applying gates to a product state, graph states, dimer states, Slater
# determinants and Fermi seas, and the graphs they are built on.
#
# The testsets of 1.x on the prepared states of layer 1 checked with expect are still in
# v1/algorithms.jl.

# a site holding a fermion beside a level of its own, which its annihilation operator leaves
# alone: two states have no fermion, and no single vacuum
struct Doubled <: AbstractSite end

TensorMixedStates.dim(::Doubled) = 4

TensorMixedStates.fermion_parities(::Doubled) = (0, 1, 0, 1)

@def_operators Doubled [ Cd = [0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 1. ; 0. 0. 0. 0.] ]

TensorMixedStates.fermion_species(::Doubled) = (Cd,)

@testset "Graph utilities" begin
    @test line_graph(4) == [(1, 2), (2, 3), (3, 4)]
    @test circle_graph(4) == [(1, 2), (2, 3), (3, 4), (4, 1)]
    @test circle_graph(2) == [(1, 2), (2, 1)]
    @test_throws "a ring has two vertices at least" circle_graph(1)
    @test complete_graph(4) == [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    @test graph_base_size(circle_graph(7)) == 7
    # a graph without edges, which a reduction over nothing refused
    @test graph_base_size(line_graph(1)) == graph_base_size(complete_graph(1)) == 0
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

    # periodic conditions: the bonds closing a direction come after the others, twice for a
    # direction of length 2, as in circle_graph. A cylinder keeps the bound on the bonds, and on
    # a torus every site has four neighbours
    @test square_lattice(3, 2; yperiodic = true) == [square_lattice(3, 2) ; [(1, 2), (3, 4), (5, 6)]]
    @test square_lattice(2, 3; xperiodic = true) == [square_lattice(2, 3) ; [(1, 6), (2, 5), (3, 4)]]
    @test square_lattice(3; yperiodic = true) == square_lattice(3, 3; yperiodic = true)
    for (nx, ny) in [(3, 3), (4, 3), (5, 4)]
        cylinder = square_lattice(nx, ny; yperiodic = true)
        @test length(cylinder) == length(square_lattice(nx, ny)) + nx && allunique(cylinder)
        @test maximum(abs(b - a) for (a, b) in cylinder) ≤ 2ny - 1
        torus = square_lattice(nx, ny; xperiodic = true, yperiodic = true)
        @test length(torus) == 2nx * ny && allunique(torus)
        @test all(k -> count(e -> k in e, torus) == 4, 1:nx * ny)
    end
    @test_throws "a ring has two vertices at least" square_lattice(1, 3; xperiodic = true)
    @test_throws "a ring has two vertices at least" square_lattice(3, 1; yperiodic = true)
end

@testset "Graph states" begin
    # the graph state of the complete graph of 4 vertices, pure and mixed: no Pauli of one site
    # has a value, its stabilizers have 1, and so does every Y Y. 1.x checked it through
    # create_graph_state, whose phases come with layer 5
    ψ = graph_state(complete_graph(4))
    for s in (ψ, mix(ψ))
        for a in (X, Y, Z)
            @test expect1(s, a) ≈ zeros(4) atol = 1e-12
        end
        @test expect(s, X(1) * Z(2) * Z(3) * Z(4)) ≈ 1
        @test expect2(s, (Y, Y)) ≈ ones(4, 4)
    end
    @test_throws "a graph with edges" graph_state(Tuple{Int, Int}[])
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
    ψ = slater_state(System(5, Electron(conserve = (TMS.strong(Ntot), 2Sz))), up, dn)
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
    # the reference state is the vacuum of each site, the state its annihilation operators all
    # send to zero, which 1.x took as the state named "0": a site where it is not unique
    @test_throws "2 states that its fermion species" slater_state(System(2, Doubled()), orth(2, 1))

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
    @test_throws "needs them orthonormal" dimer_state(System(2, Qubit()), [(1, 2)], [2., 0.],
                                                      [0., 1.])
end
