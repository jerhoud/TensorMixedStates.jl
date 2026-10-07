# Functions that extract numbers from a state, checked against exact values.
#
# Goes here: sampling, one and two point correlations, entanglement and entropy
# measures. These are the functions with an analytic answer on simple states, so tests
# here should compare against that answer rather than against a recorded output.

@testset "Qubit sampling" begin
    # `sample` takes the generator to draw from, and these tests pass one of their own
    # rather than leaning on the global one: the frequencies below then come out the same
    # whether the whole suite runs or only this group, so the tolerances are met or missed
    # once and for all instead of now and then
    rng = Xoshiro(20260920)

    sys = System(3, Qubit())
    st = State{Pure}(sys, ["Up", "Dn", "+"])  # site 3 is a 50/50 superposition

    # deterministic sites always give the same outcome
    for _ in 1:20
        @test sample(st, 1; rng) == 0
        @test sample(st, 2; rng) == 1
    end

    # superposed site: frequency should be close to 1/2
    n = 2000
    s3 = [sample(st, 3; rng) for _ in 1:n]
    @test all(x -> x in (0, 1), s3)
    @test isapprox(sum(s3) / n, 0.5; atol = 0.05)

    # sampling the whole state is consistent with per-site sampling
    for _ in 1:20
        r = sample(st; rng)
        @test r[1] == 0 && r[2] == 1 && r[3] in (0, 1)
    end

    # classical diagonal mixture with known probabilities
    p0 = 0.2
    stm = State{Mixed}(System(1, Qubit()), [p0 0. ; 0. 1 - p0])

    nm = 4000
    samples = [sample(stm, 1; rng) for _ in 1:nm]
    @test all(x -> x in (0, 1), samples)
    @test isapprox(count(==(0), samples) / nm, p0; atol = 0.03)

    samples_full = [sample(stm; rng)[1] for _ in 1:nm]
    @test isapprox(count(==(0), samples_full) / nm, p0; atol = 0.03)

    # the vector carried along a mixed chain held the probability of the outcomes drawn so far,
    # which fell below the smallest float after some 1074 qubits: every later site came out in
    # its last state. A fully mixed chain draws ones and zeros to its end
    far = sample(State{Mixed}(System(1200, Qubit()), "FullyMixed"); rng)
    @test 0 < sum(far[1101:1200]) < 100
end

@testset "Sampling a site of a charged state" begin
    # the site legs of a charged state have a direction, which the projector on a basis state
    # has to take
    sq = State{Pure}(System(3, Fermion(conserve = N)), ["Occ", "Emp", "Occ"])
    @test sample(sq, 1) == 1
    @test sample(sq, 2) == 0
end

@testset "Fermionic correlations" begin
    # dag(C) must be recognized as fermionic, and a product of two
    # fermionic operators must not be
    @test isfermionic(C)
    @test isfermionic(dag(C))
    @test !isfermionic(dag(C) * C)
    @test !isfermionic(N)
    # the function of ITensors, extended rather than shadowed: after `using ITensors` the
    # two exported names clashed and isfermionic was not defined at all
    @test isfermionic === TensorMixedStates.ITensors.isfermionic
    # expect2 inserts the Jordan-Wigner strings itself: check the whole
    # correlation matrix, at every distance, against an exact computation
    # a deterministic state with non zero correlations at every distance:
    # a product state let to spread by a hopping hamiltonian
    n = 6
    sys = System(n, Fermion())
    st = tdvp(-im * sum(dag(C)(i)C(i + 1) + dag(C)(i + 1)C(i) for i in 1:n - 1), 0.7,
              State{Pure}(sys, ["1", "0", "1", "0", "1", "0"]); limits = Limits(maxdim = 32))
    ref = exact_fermionic_correlations(st, Fermion())
    @test ref ≈ ref'                                    # sanity check on the reference
    @test real([ref[i, i] for i in 1:n]) ≈ real(expect1(st, N))
    @test expect2(st, (dag(C), C)) ≈ ref atol=1e-10
    @test expect2(st, (C, dag(C))) ≈ [ i == j ? 1 - ref[i, i] : -ref[j, i] for i in 1:n, j in 1:n ] atol=1e-10
    # a single fermionic operator is odd, so its expectation value is 0 on any physical
    # state; expect1 says so rather than computing a column of zeros
    @test_throws "vanishes on any state of definite fermion parity" expect1(st, C)
    @test expect(st, C(2)) ≈ 0 atol=1e-10          # the route the message points to
    # a pair mixing the two parities is not a correlation anyone can ask for: its product
    # is odd, and the branch that inserts the Jordan-Wigner string is chosen on the first
    # operator alone, so it has to be refused by name rather than by accident
    @test_throws "one is fermionic and the other is not" expect2(st, (dag(C), N))
    @test_throws "one is fermionic and the other is not" expect2(st, (N, dag(C)))
    @test_throws "one is fermionic and the other is not" expect2(st, [(N, N), (C, Id)])
    # but mixing fermionic and non fermionic pairs in one call must not disturb either
    m = expect2(st, [(N, N), (dag(C), C)])
    @test m[1] ≈ expect2(st, (N, N)) atol=1e-10
    @test m[2] ≈ ref atol=1e-10
    # the mixed representation must agree with the pure one, computed by another algorithm
    @test expect2(mix(st), (dag(C), C)) ≈ ref atol=1e-8
    b = RandomState{Pure}(System(4, Boson(3)), 4)
    # the entries (j, i) of equal operators, and of adjoint ones on a pure state, follow from (i, j)
    pairs = [(A, dag(A)), (N, A), (dag(A), N), (A, A), (N, N)]
    @test all(expect2(b, pairs) .≈ expect2(mix(b), pairs))
    ψp = RandomState(State{Pure}(System(n, Fermion(conserve = parity(N))), "Emp"), 6)
    cc = [ expect(ψp, C(i) * C(j)) for i in 1:n, j in 1:n ]
    @test maximum(abs, cc) > 0.01
    @test expect2(ψp, (C, C)) ≈ cc atol=1e-10
    # expect normalises its argument itself, so an operator already put through simplify
    # and the same product written directly must agree, and both must match the reference.
    # What simplify inserts here is the Jordan-Wigner strings (Multi_F for two sites or
    # more), which is exactly what a bare product would otherwise be missing.
    stm = mix(st)
    for d in 1:n - 1
        op = simplify(dag(C)(1) * C(1 + d))
        @test expect(st, op) ≈ ref[1, 1 + d] atol=1e-10
        @test expect(stm, op) ≈ ref[1, 1 + d] atol=1e-8
        @test expect(st, dag(C)(1) * C(1 + d)) ≈ ref[1, 1 + d] atol=1e-10
        @test expect(stm, dag(C)(1) * C(1 + d)) ≈ ref[1, 1 + d] atol=1e-8
        # factors in descending order used to recontract a tensor already consumed, and
        # swapping two fermionic operators costs the anticommutation sign
        @test expect(st, C(1 + d) * dag(C)(1)) ≈ -ref[1, 1 + d] atol=1e-10
        # a sum on the first site is gathered into one factor, which the string of the
        # later operator crosses. <c_j c_1> vanishes, the number of particles being fixed
        @test expect(st, C(1 + d) * (C + dag(C))(1)) ≈ -ref[1, 1 + d] atol=1e-10
        @test expect(stm, C(1 + d) * (C + dag(C))(1)) ≈ -ref[1, 1 + d] atol=1e-8
        # a renamed fermionic operator keeps its string, on either side of the product
        @test expect(st, C(1 + d) * named(dag(C), "Cd")(1)) ≈ -ref[1, 1 + d] atol=1e-10
        @test expect(st, dag(C)(1) * named(2C, "C2")(1 + d)) ≈ 2ref[1, 1 + d] atol=1e-10
    end
end

@testset "Electron and Tj against Jordan-Wigner strings written by hand" begin
    # two species on a site, whose order on the site the matrices of Cup and Cdn hold, and
    # strings through the whole sites on the left
    n = 4
    for site in (Electron(), Tj())
        st = RandomState{Pure}(System(n, site), 4)
        psi = dense_vector(st)
        exact(factors) = psi' * jw_matrix(site, n, factors) * psi
        for a in (Cup, Cdn), b in (Cup, Cdn)
            ref = [ exact([dag(a) => i, b => j]) for i in 1:n, j in 1:n ]
            @test [ expect(st, dag(a)(i) * b(j)) for i in 1:n, j in 1:n ] ≈ ref atol = 1e-10
            @test expect2(st, (dag(a), b)) ≈ ref atol = 1e-10
        end
        # a pair hopping, four fermionic factors on two sites apart
        @test expect(st, dag(Cup)(1) * dag(Cdn)(1) * Cdn(3) * Cup(3)) ≈
              exact([dag(Cup) => 1, dag(Cdn) => 1, Cdn => 3, Cup => 3]) atol = 1e-10
    end
end

@testset "Superoperators of placed operators, unnormalized expectations and matrix elements" begin
    strong = TensorMixedStates.strong
    sys = System(4, Fermion())
    ψ = normalize(RandomState{Pure}(ComplexF64, sys, 4))
    ρ = mix(ψ)
    a = dag(C)(1) * C(3) + 0.5im * N(2)
    adjoint_a = dag(C)(3) * C(1) - 0.5im * N(2)
    b = dag(C)(3) * C(1) + N(4)
    # Left and Right of a placed operator, a sum with its strings and its coefficients, against
    # the pure state: tr(b a ρ) = <ψ|b a|ψ> and tr(b ρ a†) = <ψ|a† b|ψ>
    @test expect(apply(make_mpo(ρ, Left(a)), ρ), b; normalize = false) ≈ inner(ψ, b * a, ψ)
    @test expect(apply(make_mpo(ρ, Right(a)), ρ), b; normalize = false) ≈
          inner(ψ, adjoint_a * b, ψ)
    @test norm(apply(Left(dag(C)(1) * C(3)), ρ) - apply(Left(dag(C))(1) * Left(C)(3), ρ)) < 1e-12

    # not divided by the trace: of the identity, of a pure state, and of c_1 ρ, of trace zero on
    # a state of a definite number of fermions, which the normalized one cannot divide by
    @test expect(2 * ρ, 3 * Id(1); normalize = false) ≈ 6
    @test expect(2 * ψ, N(1); normalize = false) ≈ 4 * expect(ψ, N(1))
    for site in (Fermion(), Fermion(conserve = N), Fermion(conserve = strong(N)))
        s = System(4, site)
        φ = normalize(State{Pure}(s, ["Occ", "Emp", "Occ", "Emp"]) +
                      0.5 * State{Pure}(s, ["Emp", "Occ", "Occ", "Emp"]))
        σ = apply(Left(C(1)), mix(φ))
        @test expect(σ, dag(C)(1); normalize = false) ≈ expect(mix(φ), N(1))
        @test expect(σ, [dag(C)(1), dag(C)(1) * N(3)]; normalize = false) ≈ [0.8, 0.8]
    end

    # the matrix element of an operator between two states, and of a superoperator between two
    # mixed ones: <<mix(u)| Gate(X) |mix(v)>> = |<u|X|v>|^2
    q = System(4, Qubit())
    u, v = RandomState{Pure}(ComplexF64, q, 3), RandomState{Pure}(ComplexF64, q, 3)
    m = X(1) * Z(3) + 0.3 * Y(2)
    @test inner(u, m, v) ≈ inner(u, apply(make_mpo(v, m), v))
    @test dot(u, m, v) ≈ inner(u, m, v)
    @test inner(mix(u), Gate(X)(1), mix(v)) ≈ abs2(inner(u, X(1), v))
end

@testset "The adjoint of a tensor product of fermions" begin
    # (A ⊗ B)(i, j) is A(i) * B(j), whose adjoint reverses the factors: two fermions
    # anticommute, so dag(C ⊗ C) is the opposite of dag(C) ⊗ dag(C)
    st = RandomState{Pure}(System(3, Fermion()), 4)
    for (i, j) in [(1, 2), (1, 3), (3, 1)]
        @test expect(st, dag(C ⊗ C)(i, j)) ≈ conj(expect(st, C(i) * C(j)))
        @test expect(st, dag(dag(C) ⊗ C)(i, j)) ≈ conj(expect(st, dag(C)(i) * C(j)))
    end
    # a hopping term written with its adjoint is self adjoint, and the dissipator of a
    # fermionic jump preserves the trace
    h = dag(C) ⊗ C
    @test abs(imag(expect(st, sum((h + dag(h))(i, i + 1) for i in 1:2)))) < 1e-12
    ρ = mix(st)
    for l in (C ⊗ C, dag(C) ⊗ C, C ⊗ dag(C))
        @test abs(trace(apply(make_mpo(ρ, Dissipator(l)(1, 2)), ρ))) < 1e-12
    end
end

@testset "Functions of a fermionic operator" begin
    # exp and mod of an operator that is not even, placed after a fermionic site, take the
    # string on their odd part, which is checked against the dense matrices of the explicit
    # c₁ = c ⊗ 1 ⊗ 1, c₂ = F ⊗ c ⊗ 1 and c₃ = F ⊗ F ⊗ c
    fe = Fermion()
    st = RandomState{Pure}(System(3, fe), 4)
    idx = [ SysIndex{Pure}(st.system, k) for k in 1:3 ]
    psi = reshape(Array(reduce(*, [st.state[k] for k in 1:3]), reverse(idx)...), 8)
    ev(m) = psi' * m * psi / (psi' * psi)
    c, f, id = matrix(C, fe), matrix(F, fe), matrix(Id, fe)
    cs = [kron(c, id, id), kron(f, c, id), kron(f, f, c)]
    g = 0.7 * (C + dag(C))
    for i in 1:3
        b = 0.7 * (cs[i] + cs[i]')
        @test expect(st, exp(g)(i)) ≈ ev(exp(b))
        @test expect(st, mod(g, 3)(i)) ≈ ev(exp(2im * π * b / 3))
        @test expect(st, exp(C)(i)) ≈ ev(kron(id, id, id) + cs[i])
        @test expect(mix(st), exp(g)(i)) ≈ ev(exp(b))
        # a non integer power is a function of the operator too
        @test expect(st, (g^0.5)(i)) ≈ ev(b^0.5)
        # simplifying the result again leaves it as it is
        @test simplify(simplify(exp(g)(i))) == simplify(exp(g)(i))
    end
end

@testset "A projector on a fermionic site has a definite parity" begin
    # a projector is taken as even, the F of a string commuting across it: one on a state of
    # mixed parity, which is no observable of the mode, is refused
    fe = Fermion()
    st = RandomState{Pure}(System(3, fe), 4)
    idx = [ SysIndex{Pure}(st.system, k) for k in 1:3 ]
    psi = reshape(Array(reduce(*, [st.state[k] for k in 1:3]), reverse(idx)...), 8)
    ev(m) = psi' * m * psi / (psi' * psi)
    c, f, id = matrix(C, fe), matrix(F, fe), matrix(Id, fe)
    @test_throws "no definite fermionic parity" expect(st, C(3) * Proj([1., 1.] / sqrt(2))(1))
    p = Proj([0., 1.])
    @test matrix(simplify(F * p), fe) ≈ matrix(F * p, fe)
    pm = kron([0. 0. ; 0. 1.], id, id)
    c3 = kron(f, f, c)
    @test expect(st, C(3) * p(1)) ≈ ev(c3 * pm)
    @test expect(st, p(1) * C(3)) ≈ ev(pm * c3)
end

@testset "Entanglement and entropies" begin
    L2 = log(2)
    s2 = System(2, Qubit())
    s4 = System(4, Qubit())
    # a mindim above the Schmidt rank keeps singular values of exactly zero, whose 0 * log(0)
    # made the entropy NaN
    _, g = dmrg(sum(Z(i) for i in 1:4), State{Pure}(s4, "Up"); nsweeps = 2,
                limits = Limits(mindim = 3, maxdim = 3))
    @test entanglement_entropy(g, 2)[1] ≈ 0 atol=1e-12
    bell = (State{Pure}(s2, "Up") + State{Pure}(s2, "Dn")) / sqrt(2)
    ghz = (State{Pure}(s4, "Up") + State{Pure}(s4, "Dn")) / sqrt(2)
    prod = State{Pure}(s4, "Up")
    fullymixed = State{Mixed}(s4, "FullyMixed")

    # a maximally entangled pair carries one bit, with a flat spectrum
    e, spectrum = entanglement_entropy(bell, 1)
    @test e ≈ L2
    @test spectrum ≈ [0.5, 0.5]
    # a GHZ state carries that same one bit wherever it is cut
    for cut in 1:3
        @test entanglement_entropy(ghz, cut)[1] ≈ L2
    end
    # the cut is on the right of the site, so the last site separates the whole state from
    # nothing, and a position outside the chain is no cut at all
    @test entanglement_entropy(ghz, 4)[1] ≈ 0 atol=1e-12
    @test_throws "entanglement entropy at site 0 of a 4 site state" entanglement_entropy(ghz, 0)
    @test_throws "entanglement entropy at site 5 of a 4 site state" entanglement_entropy(ghz, 5)
    # a product state carries none
    e, spectrum = entanglement_entropy(prod, 2)
    @test e ≈ 0 atol=1e-12
    @test spectrum ≈ [1.]

    # a pure state stays pure, whichever representation holds it
    for st in [bell, ghz, mix(ghz)]
        @test trace2(st) ≈ 1
        @test renyi2(st) ≈ 0 atol=1e-12
    end
    # but one of its sites, seen alone, is maximally mixed
    @test renyi2(mix(ghz), [1]) ≈ L2
    @test renyi2(mix(ghz), [1, 2]) ≈ L2
    # and the pure representation must answer the same: a subsystem of a pure state is not
    # pure, so this is an entanglement measure and not the 0. the whole state gives
    @test renyi2(ghz, [1]) ≈ L2
    @test renyi2(ghz, [1, 2]) ≈ L2
    @test renyi2(ghz, [1, 3]) ≈ L2                  # the sites need not be contiguous
    @test renyi2(prod, [1, 2]) ≈ 0 atol=1e-12
    # across a cut the mutual information is twice it, which is how it is computed there
    @test mutual_info_renyi2(bell, 1) ≈ 2 * renyi2(bell, [1])

    # the infinite temperature state of n qubits has purity 2^-n and entropy n log 2
    @test trace2(fullymixed) ≈ 1 / 16
    @test renyi2(fullymixed) ≈ 4L2
    @test renyi2(fullymixed, [1]) ≈ L2
    @test hermiticity(fullymixed) ≈ 1

    # mutual information: 2 log 2 for GHZ, none between independent halves
    @test mutual_info_renyi2(ghz, 2) ≈ 2L2
    @test mutual_info_renyi2(mix(ghz), 2) ≈ 2L2
    @test mutual_info_renyi2(fullymixed, 2) ≈ 0 atol=1e-12

    # on a pure state the mutual information across a cut is read from the entanglement
    # spectrum rather than from partial traces: both routes must agree, on a state whose
    # spectrum is not degenerate
    k = 5
    spread = tdvp(-im * (sum(Z(i)Z(i + 1) for i in 1:k - 1) + sum(0.7X(i) for i in 1:k)),
                  1.3, State{Pure}(System(k, Qubit()), "Z+"); limits = Limits(maxdim = 32))
    for cut in 1:k - 1
        @test mutual_info_renyi2(spread, cut) ≈ mutual_info_renyi2(mix(spread), cut)
        # the cut form must agree with the list form on the same bipartition
        @test mutual_info_renyi2(spread, cut) ≈ mutual_info_renyi2(spread, collect(1:cut))
    end
    # a non contiguous part of a pure state is twice its entropy as well, one partial trace
    # standing for the two parts and the whole of the mixed state
    @test mutual_info_renyi2(spread, [1, 3]) ≈ mutual_info_renyi2(mix(spread), [1, 3])
    # and so on fermions, for a state of definite parity
    ψf = RandomState(State{Pure}(System(6, Fermion(conserve = parity(N))), "Emp"), 8)
    @test mutual_info_renyi2(ψf, [2, 5]) ≈ mutual_info_renyi2(mix(ψf), [2, 5])

    # the same quantities as measurements must agree with the functions
    m(state, f) = only(last.(measure(state, Measure(f))))
    @test m(ghz, Trace) ≈ 1
    @test m(ghz, Norm) ≈ 1
    @test m(ghz, Trace2) ≈ 1
    @test m(ghz, Purity) ≈ 1
    @test m(ghz, Renyi2) ≈ 0 atol=1e-12
    @test m(ghz, Hermiticity) ≈ 1
    @test m(ghz, TraceError) ≈ 0 atol=1e-12
    @test m(ghz, HermiticityError) ≈ 0 atol=1e-12
    @test m(ghz, MaxLinkdim) == 2
    @test m(ghz, EntanglementEntropy(2)) ≈ L2
    @test m(fullymixed, Purity) ≈ 1 / 16
    @test m(fullymixed, Renyi2) ≈ 4L2
    @test m(mix(ghz), SubRenyi2([1])) ≈ L2
    @test m(ghz, MutualInfoRenyi2(2)) ≈ 2L2
    @test m(mix(ghz), MutualInfoRenyi2(2)) ≈ 2L2
    @test MutualInfoRenyi2(2).name == "MutualInfoRenyi2(1,2)"
    # a link is named after the sites on its left, which it stands for: MutualInfoRenyi2(3)
    # was named as the one site part [3], so that the two could not be measured together
    st4 = RandomState{Pure}(System(4, Qubit()), 3)
    r = Dict(measure(st4, [MutualInfoRenyi2(3), MutualInfoRenyi2([3])]))
    @test r["MutualInfoRenyi2(1:3)"] ≈ mutual_info_renyi2(st4, [1, 2, 3])
    @test r["MutualInfoRenyi2(3)"] ≈ mutual_info_renyi2(st4, [3])
    # the spectrum asked for is written whole, the eigenvalues beyond the bond dimension at the
    # cut being zeros, so that the rows of a file keep their width
    v = last(only(measure(State{Pure}(System(3, Qubit()), "Up"), EntanglementEntropy(1, 4))))
    @test length(v) == 5
    @test v[3:5] == zeros(3)
end

@testset "Entanglement resolved by sector" begin
    IT = TensorMixedStates.ITensors
    # the sectors add up to the entanglement entropy, the charge fluctuating between the two
    # sides bringing the number entropy on top of what each sector holds
    total(d) = sum(v.weight * v.entropy for v in values(d)) -
               sum(v.weight * log(v.weight) for v in values(d))
    f = System(4, Fermion(conserve = N))
    p(v) = State{Pure}(f, v)
    st = (p(["Occ", "Emp", "Emp", "Occ"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
    for pos in 1:4
        @test total(entanglement_by_sector(st, pos)) ≈ first(entanglement_entropy(st, pos))
    end

    # cut after one site, the charge on the left is either one particle or none
    d = entanglement_by_sector(st, 1)
    @test sort(collect(keys(d)); by = string) == [IT.QN("N", 0), IT.QN("N", 1)]
    @test d[IT.QN("N", 1)].weight ≈ 0.8
    @test d[IT.QN("N", 0)].weight ≈ 0.2
    # cut after two, it is one particle in both terms, and the entanglement lives inside it
    d = entanglement_by_sector(st, 2)
    @test collect(keys(d)) == [IT.QN("N", 1)]
    @test d[IT.QN("N", 1)].weight ≈ 1
    @test d[IT.QN("N", 1)].spectrum ≈ [0.8, 0.2]
    @test d[IT.QN("N", 1)].entropy ≈ first(entanglement_entropy(st, 2))

    # several quantities at once, and a charge modulo 2 over a basis whose charges repeat
    e = System(3, Electron(conserve = (Ntot, 2Sz)))
    q(v) = State{Pure}(e, v)
    d = entanglement_by_sector((q(["Up", "Dn", "Emp"]) + 0.5 * q(["UpDn", "Emp", "Emp"])) / sqrt(1.25), 1)
    @test d[IT.QN(("Ntot", 1), ("2Sz", 1))].weight ≈ 0.8
    @test d[IT.QN(("Ntot", 2), ("2Sz", 0))].weight ≈ 0.2
    b = System(3, Boson(4, conserve = parity(N)))
    r(v) = State{Pure}(b, v)
    d = entanglement_by_sector((r(["1", "2", "0"]) + 0.5 * r(["0", "3", "0"])) / sqrt(1.25), 1)
    @test d[IT.QN("parity(N)", 1, 2)].weight ≈ 0.8
    @test d[IT.QN("parity(N)", 0, 2)].weight ≈ 0.2

    # without charges there is a single sector holding the whole spectrum
    fd = System(4, Fermion())
    pd(v) = State{Pure}(fd, v)
    sd = (pd(["Occ", "Emp", "Emp", "Occ"]) + 0.5 * pd(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
    d = entanglement_by_sector(sd, 2)
    @test collect(keys(d)) == [IT.QN()]
    @test d[IT.QN()].spectrum ≈ last(entanglement_entropy(sd, 2))

    # a mixed representation is refused, its links pairing the charge of a ket with a bra
    @test_throws "takes a pure state" entanglement_by_sector(mix(st), 1)
    @test_throws "between 1 and 4" entanglement_by_sector(st, 5)
end

@testset "Composed operators as measurements" begin
    # `SimpleOp` is abstract, so a measurement may be any pure one site operator and not
    # just a named one. Only `Operator` has a name of its own; every other form is a
    # composition and is named by how it prints
    sys = System(3, Qubit())
    st = State{Pure}(sys, ["X+", "Y+", "Z+"])
    colnames(o) = first.(measure(st, Measure([o]), 0.))
    # a superoperator is not an observable, and expect says so
    @test_throws "is a superoperator" expect(mix(st), Left(X)(1))
    @test colnames(X) == ["X"]
    @test colnames(X * Y) == ["X*Y"]
    @test colnames(2X) == ["2X"]
    @test colnames(X + Y) == ["X+Y"]
    @test colnames(dag(S)) == ["dag(S)"]
    # these two used to need a method of their own, for want of a name field
    @test colnames(Id) == ["Id"]
    @test colnames((Id, X)) == ["IdX"]
    @test colnames((X, Y)) == ["XY"]
    # a coefficient rides outside the AtIndex that `(c * A)(i)` builds, so it has to come
    # off where the one site tensor is made and not only where the name is
    @test expect1(st, 2X) ≈ 2 * expect1(st, X)
    @test expect1(mix(st), 2X) ≈ 2 * expect1(mix(st), X)
    @test expect1(st, 0.5 * (X + Y)) ≈ 0.5 * (expect1(st, X) + expect1(st, Y))
    @test expect2(st, (2X, 3Y)) ≈ 6 * expect2(st, (X, Y))
end

@testset "Expectation values on an unnormalised state" begin
    # The left environments carry the 1/trace normalisation of `expect`, injected once in
    # l[1] and riding along the recursion. The shortcut taken when the mps is already left
    # orthogonal writes the environment from scratch and has to carry it too. Two non
    # unitary gates produce both conditions at once: a norm that is not one, and an
    # orthogonality centre pushed to the right, Sp(1) acting first and Sp(3) after it.
    sys = System(6, Qubit())
    st = apply(Sp(3) * Sp(1), State{Pure}(sys, "+"))
    # the test only means anything while the shortcut is actually taken, and that depends
    # on where ITensor leaves the orthogonality centre, so it is asserted rather than hoped
    @test TensorMixedStates.ITensorMPS.leftlim(st.state) ≥ 1
    @test trace(st) ≈ 0.25
    # Sp|+> = |0>/sqrt(2), so sites 1 and 3 are up and the other four are still |+>
    @test expect1(st, Z) ≈ [1., 0., 1., 0., 0., 0.] atol = 1e-12
    @test expect1(st, X) ≈ [0., 1., 0., 1., 1., 1.] atol = 1e-12
    # the answer must not depend on whether the caller normalised first. Only a term whose
    # leftmost factor sits on site 1 reads l[1], which is why the bug showed on expect1 and
    # expect2 while a product starting at site 1 stayed right
    @test expect1(st, X) ≈ expect1(normalize(st), X)
    @test expect(st, Z(1) * Z(3)) ≈ expect(normalize(st), Z(1) * Z(3))
    @test expect2(st, (X, X)) ≈ expect2(normalize(st), (X, X))

    # on a charged state the shortcut has arrows to respect, and it had them the other way
    # round: an ordinary gate pushing the orthogonality centre to the right was enough to make
    # every measurement starting past site 1 fail
    a = apply(Swap(2, 3), State{Pure}(System(4, Qubit(conserve = N)), ["Up", "Dn", "Up", "Dn"]))
    @test TensorMixedStates.ITensorMPS.leftlim(a.state) ≥ 1
    @test expect1(a, N) ≈ [0., 0., 1., 1.] atol = 1e-12
    @test expect(a, N(3) * N(4)) ≈ 1
end

@testset "Positions given as a range" begin
    # the position arguments took a `Vector{Int}` and nothing else, so a range was refused
    # by a `MethodError`. The name of the measurement was built happily by
    # `compact_positions`, so the failure came only when the measurement was taken
    sys = System(6, Qubit())
    stp = RandomState{Pure}(sys, 4)
    stm = mix(stp)
    @test renyi2(stm, 1:3) ≈ renyi2(stm, [1, 2, 3])
    @test renyi2(stp, 1:3) ≈ renyi2(stp, [1, 2, 3])
    @test mutual_info_renyi2(stm, 1:3) ≈ mutual_info_renyi2(stm, [1, 2, 3])
    @test length(partial_trace(stm, 1:3; keep = true)) == 3
    # `measure` hands back a vector of name => value pairs, one per measurement
    @test last(only(measure(stm, SubRenyi2(1:3)))) ≈ renyi2(stm, [1, 2, 3])
    @test last(only(measure(stm, MutualInfoRenyi2(1:3)))) ≈ mutual_info_renyi2(stm, [1, 2, 3])
    # a partial trace needs a density matrix, and says so rather than raising a MethodError
    @test_throws "mixed representation" partial_trace(stp, [1, 2])
end

@testset "A partial trace keeps the trace" begin
    # the first site was taken from the left environment of expect, normalised by the trace:
    # the result had trace 1 whatever the state, and a traceless state gave NaN
    stm = mix(RandomState{Pure}(System(4, Qubit()), 4))
    for k in ([1], [2], [4], [1, 3])
        @test trace(partial_trace(2 * stm, k)) ≈ 2 * trace(stm)
    end
    @test abs(trace(partial_trace(stm - stm, [1]))) < 1e-12
end

@testset "A check on a symbol not given" begin
    # an empty value, which a symbol not given has, was subtracted from the other: the check
    # failed on a MethodError, and with it a whole simulation measuring it at its end
    st = State{Pure}(System(2, Qubit()), "Up")
    @test last(only(measure(st, Check("e", :energy, -1.0)))) == Any[[], -1.0, []]
    @test last(only(measure(st, Check("e", :energy, -1.0); energy = -1.0))) == Any[-1.0, -1.0, 0.0]
    @test_throws "Check e has nothing to compare" measure(st, Check("e", :energy, -1.0, 1e-8))
end

@testset "Measuring a multiple of the identity" begin
    # a one site operator simplified to c Id went to the tensor of an identity placed on no
    # site, which has none: expect1, expect2 and measure failed on it, Id alone measured right
    st = RandomState{Pure}(System(3, Qubit()), 4)
    for s in (st, mix(st))
        @test expect1(s, 2Id) ≈ fill(2, 3)
        @test expect1(s, (2X)^2) ≈ fill(4, 3)
        @test expect2(s, (2Id, X)) ≈ 2 * expect2(s, (Id, X))
        @test last(only(measure(s, (0.5Z)^2))) ≈ fill(0.25, 3)
    end
end

@testset "A partial trace keeps the fermionic signs" begin
    # a fermion traced out has to be moved past the fermions kept on its right, which the
    # trace of the spins forgot: on (|110⟩ + |011⟩)/√2 the hopping from 3 to 1 crosses the
    # fermion of site 2, and gave +1/2 on the reduced state instead of -1/2
    fe = Fermion()
    s3 = System(3, fe)
    ρ = mix(normalize(State{Pure}(s3, ["Occ", "Occ", "Emp"]) + State{Pure}(s3, ["Emp", "Occ", "Occ"])))
    @test expect(partial_trace(ρ, [2]), dag(C)(1) * C(2)) ≈ expect(ρ, dag(C)(1) * C(3)) ≈ -0.5
    # a state of no definite parity takes the sign on a single fermionic operator as well
    ρ = mix(normalize(State{Pure}(s3, ["Emp", "Occ", "Emp"]) + State{Pure}(s3, ["Emp", "Occ", "Occ"])))
    @test expect(partial_trace(ρ, [1, 2]), C(1) + dag(C)(1)) ≈ expect(ρ, C(3) + dag(C)(3)) ≈ -1
    # the reduced state is the right one: every operator of the sites kept has on it the value
    # it has on the state. With sites traced out on both sides, the entropies of the trace of
    # the spins were wrong too, for a state of definite parity
    ρ = mix(RandomState(State{Pure}(System(5, Fermion(conserve = N)),
                                    ["Occ", "Occ", "Emp", "Emp", "Occ"]), 4))
    red = partial_trace(ρ, [1, 4]; keep = true)
    odd = (C, dag(C))
    for a in (Id, N, C, dag(C)), b in (Id, N, C, dag(C))
        if (a in odd) == (b in odd)
            @test expect(red, a(1) * b(2)) ≈ expect(ρ, a(1) * b(4)) atol = 1e-12
        end
    end
    es = System(3, Electron())
    ρ = mix(RandomState{Pure}(es, 4))
    @test expect(partial_trace(ρ, [2]), dag(Cdn)(1) * Cup(2)) ≈ expect(ρ, dag(Cdn)(1) * Cup(3))
    # no fermion kept on the right of one traced out: no sign, and the bond dimension stays
    ρ = mix(RandomState{Pure}(System(4, fe), 4))
    @test maxlinkdim(partial_trace(ρ, [3, 4])) ≤ maxlinkdim(ρ)
end

@testset "Reduced density matrix and von Neumann entropy" begin
    # tr(ρ_A M) is the expectation value of the operator whose matrix on the kept sites is M,
    # in the order of the sites and the Jordan-Wigner basis of the kept sites alone, which
    # also checks the order of kets and bras with operators that are not hermitian
    LA = TensorMixedStates.LinearAlgebra
    q = RandomState{Pure}(System(4, Qubit()), 4)
    ρ = reduced_density_matrix(q, [1, 3])
    @test LA.tr(ρ) ≈ 1
    @test LA.ishermitian(round.(ρ; digits = 12))
    for (a, b) in [(X, Z), (Sp, Sm), (Y, Id), (Sm, Sp)]
        m = matrix(a ⊗ b, Qubit(), Qubit())
        @test LA.tr(ρ * m) ≈ expect(q, a(1) * b(3)) atol = 1e-12
    end
    @test reduced_density_matrix(mix(q), [1, 3]) ≈ ρ
    @test reduced_density_matrix(q, Any[3, 1]) ≈ ρ
    s = RandomState{Pure}(System([Qubit(), Spin(1), Qubit()]), 4)
    ρs = reduced_density_matrix(s, [2, 3])
    @test size(ρs) == (6, 6)
    m = matrix(Sp ⊗ X, Spin(1), Qubit())
    @test LA.tr(ρs * m) ≈ expect(s, Sp(2) * X(3)) atol = 1e-12
    f = RandomState{Pure}(System(4, Fermion()), 4)
    ρf = reduced_density_matrix(f, [1, 3])
    for (a, b) in [(dag(C), C), (C, C), (N, dag(C))]
        m = matrix(a ⊗ b, Fermion(), Fermion())
        @test LA.tr(ρf * m) ≈ expect(f, a(1) * b(3)) atol = 1e-12
    end
    strong_ρ = State{Mixed}(System(3, Fermion(conserve = strong(N))), ["Occ", "Emp", "Occ"])
    @test real(LA.diag(reduced_density_matrix(strong_ρ, [1, 2]))) ≈ [0, 0, 1, 0]
    @test reduced_density_matrix(q, []) == ones(1, 1)
    @test_throws "given site 5, which the state does not have" reduced_density_matrix(q, [5])
    # the entropy of half a Bell pair is log 2, that of the pair 0, and on a cut it is the
    # entanglement entropy
    bell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(3, Qubit()), "Up"))
    @test vonneumann_entropy(bell, [1]) ≈ log(2)
    @test abs(vonneumann_entropy(bell, [1, 2])) < 1e-12
    @test vonneumann_entropy(q, 1:2) ≈ first(entanglement_entropy(q, 2))
    @test vonneumann_entropy(q, []) === 0.0
    values = Dict(measure(q, [ReducedDensityMatrix([1, 3]), VonNeumannEntropy([1])]))
    @test values["ReducedDensityMatrix(1,3)"] ≈ ρ
    @test values["VonNeumannEntropy(1)"] ≈ vonneumann_entropy(q, [1])
end

@testset "Reduced density matrix of a pure state" begin
    # contracted on the pure state, it is the one of its mixed representation, fermionic signs
    # included: sites kept on both sides of sites traced out, and a site traced out first
    for (st, positions) in [
            (RandomState{Pure}(System(6, Qubit()), 8), [[6], [5, 2], [1, 3, 6]]),
            (RandomState(State{Pure}(System(6, Fermion(conserve = N)),
                                     ["Occ", "Emp", "Occ", "Emp", "Occ", "Emp"]), 8),
             [[2, 4], [2, 3, 5]]),
            (RandomState(State{Pure}(System(4, Electron(conserve = (Ntot, 2Sz))),
                                     ["Up", "Dn", "Emp", "UpDn"]), 6),
             [[1, 3], [2, 4]]),
            (RandomState{Pure}(System([Qubit(), Fermion(), Boson(3), Fermion(), Qubit()]), 6),
             [[2, 4], [1, 3, 5]])]
        for pos in positions
            @test reduced_density_matrix(st, pos) ≈ reduced_density_matrix(mix(st), pos)
        end
    end
end

@testset "Measuring a site" begin
    # the result is drawn with its probability, and the state left is the state projected onto
    # it and normalized, which the outcome on one half of a Bell pair fixes on the other
    rng = Xoshiro(20261007)
    plus = State{Pure}(System(2, Qubit()), "+")
    ps = probabilities(plus, 1)
    @test first.(ps) == [0, 1]
    @test last.(ps) ≈ [0.5, 0.5]
    ps = probabilities(plus, X(1))
    @test first.(ps) ≈ [-1, 1]
    @test last.(ps) ≈ [0, 1] atol = 1e-12
    @test sample(plus, X(1); rng) ≈ 1
    @test first(collapse(plus, 1; rng)) isa Int
    bell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(2, Qubit()), "Up"))
    for st in (bell, mix(bell)), _ in 1:10
        x, s = collapse(st, 1; rng)
        @test real(expect(s, Z(2))) ≈ (x == 0 ? 1 : -1)
        @test real(trace(s)) ≈ 1
    end
    r = RandomState{Pure}(System(4, Qubit()), 4)
    for st in (r, mix(r)), _ in 1:4
        x, s = collapse(st, 2; rng)
        @test norm(s - normalize(apply(Proj(x)(2), st))) < 1e-12
        λ, s = collapse(st, X(3); rng)
        @test norm(s - normalize(apply(Proj(X => λ)(3), st))) < 1e-12
    end
    up = State{Mixed}(System(1, Qubit()), "Up")
    @test count(≈(1), [ sample(up, X(1); rng) for _ in 1:2000 ]) / 2000 ≈ 0.5 atol = 0.05
    # Ntot does not tell Up from Dn, whose coherence the measurement keeps
    e = State{Pure}(System(1, Electron()), [[0., 1., 1., 0.] / sqrt(2)])
    for st in (e, mix(e))
        λ, s = collapse(st, Ntot(1); rng)
        @test λ ≈ 1
        @test norm(s - st) < 1e-12
    end
    # on a strong symmetry, and on a simulation, which keeps its time
    fs = mix(apply(sqrt(Swap)(1, 2), State{Pure}(System(2, Qubit(conserve = strong(N))), ["Dn", "Up"])))
    x, s = collapse(fs, 1; rng)
    @test real(expect(s, N(1))) ≈ x
    @test real(expect(s, N(1) + N(2))) ≈ 1
    x, sim = collapse(Simulation(bell; time = 2.), 1; rng)
    @test sim.time == 2.
    @test real(expect(sim.state, Z(2))) ≈ (x == 0 ? 1 : -1)
end

@testset "Measuring an operator of several sites" begin
    rng = Xoshiro(20261008)
    # a parity measured on |++> leaves a Bell pair, and Swap tells the singlet from the triplets
    plus = State{Pure}(System(2, Qubit()), "+")
    ps = probabilities(plus, Z(1) * Z(2))
    @test first.(ps) ≈ [-1, 1]
    @test last.(ps) ≈ [0.5, 0.5]
    for st in (plus, mix(plus))
        s, after = collapse(st, Z(1) * Z(2); rng)
        @test real(expect(after, Z(1) * Z(2))) ≈ s
        @test real(expect(after, X(1) * X(2))) ≈ 1
    end
    @test last.(probabilities(State{Pure}(System(2, Qubit()), ["Up", "Dn"]), Swap(1, 2))) ≈ [0.5, 0.5]
    # the number of excitations of three sites, from the diagonal of their density matrix
    r = RandomState{Pure}(System(5, Qubit()), 4)
    ρ = reduced_density_matrix(r, [1, 3, 5])
    counts = zeros(4)
    for (i, b) in enumerate(Iterators.product(0:1, 0:1, 0:1))
        counts[sum(b) + 1] += real(ρ[i, i])
    end
    ps = probabilities(r, N(1) + N(3) + N(5))
    @test first.(ps) ≈ 0:3
    @test last.(ps) ≈ counts
    # a hop over site 2, of spectrum -1, 0, 1: its projectors are polynomials of it, which
    # make_mpo places with their Jordan-Wigner strings, and so must the projection
    f = RandomState(State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Occ", "Emp", "Emp"]), 4)
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    projs = Dict(-1 => (h * h - h) / 2, 0 => Id(1) - h * h, 1 => (h * h + h) / 2)
    for (λ, p) in probabilities(f, h)
        @test p ≈ real(expect(f, projs[round(Int, λ)])) atol = 1e-12
    end
    for st in (f, mix(f)), _ in 1:6
        λ, s = collapse(st, h; rng)
        ref = normalize(apply(make_mpo(f, projs[round(Int, λ)]), f))
        @test norm(s - (st isa State{Mixed} ? mix(ref) : ref)) < 1e-10
    end
    # an operator moving a charge has probabilities, but no projection on a conserving system
    q = State{Pure}(System(3, Qubit(conserve = N)), ["Up", "Dn", "Up"])
    @test last.(probabilities(q, X(1))) ≈ [0.5, 0.5]
    @test first(collapse(q, Z(1) * Z(2); rng)) ≈ -1
    @test_throws "no definite flux" collapse(q, X(1); rng)
    @test_throws "commuting with the fermionic parity" probabilities(f, (C + dag(C))(1))
    @test_throws "needs a Hermitian operator" probabilities(plus, Sp(1) * Z(2))
    @test_throws "needs an operator placed on sites" probabilities(plus, Id(1))
    @test_throws "which the system does not have" sample(plus, X(7))
    @test_throws "which the state does not have" collapse(plus, 3)
end

@testset "Collapse on a given result" begin
    # the probability of the result and the state projected onto it, on each way of measuring:
    # the result post-selected
    plus = State{Pure}(System(2, Qubit()), "+")
    for st in (plus, mix(plus))
        p, after = collapse(st, Z(1) * Z(2) => -1)
        @test p ≈ 0.5
        @test real(expect(after, Z(1) * Z(2))) ≈ -1
        @test real(expect(after, X(1) * X(2))) ≈ 1
    end
    ψ = RandomState{Pure}(System(6, Qubit()), 6)
    P = Z(1) * X(3) * Y(5)
    op = N(1) + N(3) + N(5)
    r = TensorMixedStates.measured_spectrum(ψ.system, op, "dense")
    for st in (ψ, mix(ψ))
        for (s, q) in probabilities(st, P)
            p, after = collapse(st, P => s)
            ref = normalize(apply(make_mpo(ψ, (Id(1) + s * P) / 2), ψ))
            @test p ≈ q
            @test norm(after - (st isa State{Mixed} ? mix(ref) : ref)) < 1e-12
        end
        for (v, q) in probabilities(st, op)
            p, after = collapse(st, op => v)
            k = findfirst(≈(v), r[2])
            @test p ≈ q
            @test norm(after - normalize(apply(Operator{3}("P", r[3][k], selfadjoint_op)(1, 3, 5), st))) < 1e-12
        end
    end
    p, singlet = collapse(State{Pure}(System(2, Qubit()), ["Up", "Dn"]), Swap(1, 2) => -1)
    @test p ≈ 0.5
    @test real(expect(singlet, Swap(1, 2))) ≈ -1
    bell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(2, Qubit()), "Up"))
    for st in (bell, mix(bell))
        p, after = collapse(st, 1 => 1)
        @test p ≈ 0.5
        @test real(expect(after, Z(2))) ≈ -1
    end
    p, sim = collapse(Simulation(bell; time = 3.), 2 => 0)
    @test p ≈ 0.5
    @test sim.time == 3.
    up = State{Pure}(System(2, Qubit()), "Up")
    @test_throws "0.5 is not a result of measuring X(1)" collapse(plus, X(1) => 0.5)
    @test_throws "of probability" collapse(up, Z(1) => -1)
    @test_throws "is not a result" collapse(up, 1 => 2)
    @test_throws "which the state does not have" collapse(up, 3 => 0)
end

@testset "Projector of several sites" begin
    # applied as a gate, with the Jordan-Wigner strings of the sites in between, the order of
    # its sites included; written as a sum of products by named, to be measured
    ud = State{Pure}(System(2, Qubit()), ["Up", "Dn"])
    @test real(expect(ud, named(Proj(Swap => -1), "Singlet", Qubit())(1, 2))) ≈ 0.5
    for st in (ud, mix(ud))
        @test real(expect(normalize(apply(Proj(Swap => -1)(1, 2), st)), Swap(1, 2))) ≈ -1
    end
    @test_throws "acts on several sites at once" expect(ud, Proj(Swap => -1)(1, 2))
    f = RandomState(State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Occ", "Emp", "Emp"]), 4)
    projected(a, λ) = normalize(apply(make_mpo(f, (a * a + λ * a) / 2), f))
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    j = im * (dag(C)(1) * C(3) - dag(C)(3) * C(1))
    hop = dag(C) ⊗ C - C ⊗ dag(C)
    current = im * (dag(C) ⊗ C + C ⊗ dag(C))
    for (g, ref) in [(Proj(hop => 1)(1, 3), projected(h, 1)), (Proj(hop => 1)(3, 1), projected(h, 1)),
                     (Proj(current => 1)(1, 3), projected(j, 1)),
                     (Proj(current => 1)(3, 1), projected(j, -1))]
        @test norm(normalize(apply(g, f)) - ref) < 1e-12
        @test norm(normalize(apply(g, mix(f))) - mix(ref)) < 1e-12
    end
end

@testset "Measuring a product of involutions" begin
    # measured by its mean value on any number of sites, it agrees with the matrix of the
    # operator on a few, and its projection, the state plus its image, with (1 + sP)/2
    rng = Xoshiro(20261009)
    TMS = TensorMixedStates
    dense(st, op) = (r = TMS.measured_spectrum(st.system, op, "dense");
                     collect(zip(r[2], TMS.outcome_probabilities(st, r[1], r[3]))))
    ψ = RandomState{Pure}(System(6, Qubit()), 6)
    for op in (Z(1) * X(3) * Y(5), -2 * Z(1) * Z(2), X(2) * X(3) * X(4) * X(5)), st in (ψ, mix(ψ))
        @test !isnothing(TMS.involution(st.system, op))
        ps, ds = probabilities(st, op), dense(st, op)
        @test first.(ps) ≈ first.(ds)
        @test last.(ps) ≈ last.(ds)
    end
    P = Z(1) * X(3) * Y(5)
    projected(λ) = normalize(apply(make_mpo(ψ, (Id(1) + λ * P) / 2), ψ))
    for _ in 1:4
        λ, after = collapse(ψ, P; rng)
        @test norm(after - projected(λ)) < 1e-12
        λ, after = collapse(mix(ψ), P; rng)
        @test norm(after - mix(projected(λ))) < 1e-12
    end
    p = last(last(probabilities(ψ, P)))
    @test count(_ -> first(collapse(ψ, P; rng)) > 0, 1:1000) / 1000 ≈ p atol = 0.05
    # on twenty sites, the bond dimension at most doubles
    big = RandomState{Pure}(System(20, Qubit()), 8)
    string20 = prod(X(i) for i in 1:2:20) * prod(Z(i) for i in 2:2:20)
    s, after = collapse(big, string20; rng)
    @test real(expect(after, string20)) ≈ s
    @test maxlinkdim(after) ≤ 2 * maxlinkdim(big)
    # a diagonal one on a conserving system, a parity of fermions, and a strong symmetry
    q = RandomState(State{Pure}(System(6, Qubit(conserve = N)),
                                ["Up", "Dn", "Up", "Dn", "Up", "Dn"]), 4)
    Pz = Z(1) * Z(4) * Z(6)
    @test last.(probabilities(q, Pz)) ≈ last.(dense(q, Pz))
    @test isnothing(TMS.involution(q.system, X(1) * X(2)))
    s, after = collapse(q, Pz; rng)
    @test real(expect(after, Pz)) ≈ s
    f = RandomState(State{Pure}(System(5, Fermion(conserve = N)),
                                ["Occ", "Emp", "Occ", "Emp", "Occ"]), 4)
    Pf = F(1) * F(3) * F(4)
    @test last.(probabilities(f, Pf)) ≈ last.(dense(f, Pf))
    for st in (f, mix(f))
        s, after = collapse(st, Pf; rng)
        @test real(expect(after, Pf)) ≈ s
        @test real(trace(after)) ≈ 1
    end
    fs = mix(apply(sqrt(Swap)(1, 2),
                   State{Pure}(System(3, Qubit(conserve = strong(N))), ["Dn", "Up", "Up"])))
    s, after = collapse(fs, Z(1) * Z(2); rng)
    @test real(trace(after)) ≈ 1
    @test real(expect(after, N(1) + N(2) + N(3))) ≈ 1
end

@testset "Measuring a sum of operators of one site" begin
    # its totals, integers or half integers, are counted on any number of sites, in agreement
    # with the matrix of the operator on a few, and its projection by an MPO carrying the
    # partial sums with the projector on its eigenspace
    rng = Xoshiro(20261010)
    TMS = TensorMixedStates
    dense(st, op) = (r = TMS.measured_spectrum(st.system, op, "dense");
                     collect(zip(r[2], TMS.outcome_probabilities(st, r[1], r[3]))))
    ψ = RandomState{Pure}(System(6, Qubit()), 6)
    s = RandomState{Pure}(System([Qubit(), Spin(1), Qubit(), Spin(3/2)]), 6)
    for (st, op) in [(ψ, N(1) + N(3) + N(5)), (ψ, Sz(1) + Sz(2) + Sz(4)), (ψ, X(1) + X(2) + X(3)),
                     (ψ, 2 * (N(2) - N(6))), (s, Sz(1) + Sz(2) + Sz(3) + Sz(4))], x in (st, mix(st))
        @test !isnothing(TMS.counting(x.system, op, "test"))
        ps, ds = probabilities(x, op), dense(x, op)
        @test first.(ps) ≈ first.(ds)
        @test last.(ps) ≈ last.(ds) atol = 1e-12
    end
    op = N(1) + N(3) + N(5)
    r = TMS.measured_spectrum(ψ.system, op, "dense")
    projected(st, v) = normalize(apply(Operator{3}("P", r[3][findfirst(≈(v), r[2])],
                                                   selfadjoint_op)(1, 3, 5), st))
    for st in (ψ, mix(ψ)), _ in 1:3
        v, after = collapse(st, op; rng)
        @test norm(after - projected(st, v)) < 1e-12
    end
    # on twenty sites, the number of particles of the region becomes sharp
    big = RandomState{Pure}(System(30, Qubit()), 8)
    NA = sum(N(i) for i in 6:25)
    ps = probabilities(big, NA)
    @test sum(last.(ps)) ≈ 1
    @test sum(first.(ps) .* last.(ps)) ≈ real(expect(big, NA))
    v, after = collapse(big, NA; rng)
    @test real(expect(after, NA)) ≈ v
    @test real(expect(after, NA * NA)) ≈ v^2
    # on charged fermions, the total stays, and on a strong symmetry
    f = RandomState(State{Pure}(System(8, Fermion(conserve = N)),
                                ["Occ", "Emp", "Occ", "Emp", "Occ", "Emp", "Occ", "Emp"]), 6)
    Nf = sum(N(i) for i in 2:6)
    @test last(last(probabilities(f, Nf))) ≈ 0 atol = 1e-12
    for st in (f, mix(f))
        v, after = collapse(st, Nf; rng)
        @test real(expect(after, Nf)) ≈ v
        @test real(expect(after, sum(N(i) for i in 1:8))) ≈ 4
        @test real(trace(after)) ≈ 1
    end
    fs = mix(apply(sqrt(Swap)(1, 2),
                   State{Pure}(System(3, Qubit(conserve = strong(N))), ["Dn", "Up", "Dn"])))
    v, after = collapse(fs, N(1) + N(3); rng)
    @test real(expect(after, N(1) + N(3))) ≈ v
    @test real(trace(after)) ≈ 1
end

@testset "Logarithmic negativity" begin
    # log 2 for a Bell pair, and for the Werner state p Bell + (1 - p) I/4, the depolarization
    # of the pair, log((1 + 3p)/2) above p = 1/3, separable below; the partial transpose of
    # either part gives the same, which lets a part without fermionic site take it
    bell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(3, Qubit()), "Up"))
    @test log_negativity(bell, [1], [2]) ≈ log(2)
    @test abs(log_negativity(bell, [1], [3])) < 1e-12
    for p in [0.6, 0.2]
        werner = apply(depolarizing_gate(1 - p, 2)(1, 2), mix(bell))
        @test log_negativity(werner, [1], [2]) ≈ max(0, log((1 + 3p) / 2)) atol = 1e-12
    end
    q = RandomState{Pure}(System(4, Qubit()), 4)
    @test log_negativity(q, [1], [3, 4]) ≈ log_negativity(q, [3, 4], [1])
    s = RandomState{Pure}(System([Fermion(), Qubit(), Fermion()]), 4)
    @test log_negativity(s, [1], [2]) ≈ log_negativity(s, [2], [1])
    @test_throws "a part without fermionic sites" log_negativity(s, [1], [3])
    @test_throws "sites in both parts" log_negativity(q, [1, 2], [2])
    @test log_negativity(q, [], [2]) == 0
    name, value = only(measure(bell, LogNegativity([1], [2])))
    @test name == "LogNegativity(1;2)"
    @test value ≈ log(2)
end

@testset "A random fermionic oracle" begin
    # every implementation of the Jordan-Wigner strings checked on random products against
    # the dense matrices of the strings written out: simplify, through expect on a pure and on
    # a mixed state, the MPO and the gates, the gates of functions of several sites, expect2,
    # the matrix of a tensor product and partial_trace. In the dense basis C(j) is F ⊗ … ⊗ F ⊗ c ⊗ 1 ⊗ … ⊗ 1, and any other
    # operator is written from those, so that the reference knows nothing of the package
    fe = Fermion()
    mc, mf, mi = matrix(C, fe), matrix(F, fe), matrix(Id, fe)
    cd(j, n) = foldl(kron, [ k < j ? mf : k == j ? mc : mi for k in 1:n ])
    singles = [ (C, (j, n) -> cd(j, n)), (dag(C), (j, n) -> cd(j, n)'),
                (N, (j, n) -> cd(j, n)' * cd(j, n)),
                (C + 2dag(C), (j, n) -> cd(j, n) + 2cd(j, n)') ]
    count = Ref(0)
    # a factor on the positions 1:m, as an operator and as the dense matrix on the n sites of
    # the state, `site` taking a position to its site: an operator of one site, or a tensor
    # product of two, renamed or not, on two positions in either order
    function factor(m, site, n)
        a, da = rand(singles)
        i = rand(1:m)
        if m == 1 || rand(Bool)
            return a(i), da(site(i), n)
        end
        b, db = rand(singles)
        j = rand(setdiff(1:m, i))
        t = rand(Bool) ? a ⊗ b : named(a ⊗ b, "T$(count[] += 1)")
        return t(i, j), da(site(i), n) * db(site(j), n)
    end
    function product(m, site, n)
        fs = [ factor(m, site, n) for _ in 1:rand(1:3) ]
        return prod(first.(fs)), prod(last.(fs))
    end
    n = 4
    sys = System(n, fe)
    idx = [ SysIndex{Pure}(sys, k) for k in 1:n ]
    dense_vec(st) = reshape(Array(reduce(*, [ st.state[k] for k in 1:n ]), reverse(idx)...), 2^n)
    ψ = RandomState{Pure}(sys, 4)
    ρ = mix(ψ)
    v = dense_vec(ψ)
    value(d) = v' * d * v / (v' * v)
    for _ in 1:30
        o, d = product(n, identity, n)
        @test expect(ψ, o) ≈ value(d) atol = 1e-10
        @test expect(ρ, o) ≈ value(d) atol = 1e-10
        @test dense_vec(apply(make_mpo(ψ, o), ψ)) ≈ d * v atol = 1e-10
        @test dense_vec(apply(o, ψ)) ≈ d * v atol = 1e-10
        @test norm(apply(o, ρ) - mix(apply(o, ψ))) < 1e-10
        keep = sort(randperm(n)[1:rand(1:3)])
        o, d = product(length(keep), p -> keep[p], n)
        @test expect(partial_trace(ρ, keep; keep = true), o) ≈ value(d) atol = 1e-10
    end
    # a function of an even product applied as a gate, on sites in any order and apart, whose
    # strings come from a conjugation by diagonal gates
    odd = singles[[1, 2, 4]]
    for _ in 1:10
        fs = Any[ rand(odd), rand(odd) ]
        if rand(Bool)
            insert!(fs, rand(1:3), singles[3])
        end
        sites = randperm(n)[1:length(fs)]
        a = reduce(⊗, first.(fs))
        da = prod(f[2](s, n) for (f, s) in zip(fs, sites))
        g = exp(0.3 * a - 0.3 * dag(a))(sites...)
        @test dense_vec(apply(g, ψ)) ≈ exp(0.3 * da - 0.3 * da') * v atol = 1e-10
        @test norm(apply(g, ρ) - mix(apply(g, ψ))) < 1e-10
    end
    # and on electrons, whose parity counts both spins
    el = Electron()
    mu, md, mfe, mie = matrix(Cup, el), matrix(Cdn, el), matrix(F, el), matrix(Id, el)
    ce(j, m) = foldl(kron, [ k < j ? mfe : k == j ? m : mie for k in 1:3 ])
    es = System(3, el)
    eidx = [ SysIndex{Pure}(es, k) for k in 1:3 ]
    dense_e(st) = reshape(Array(reduce(*, [ st.state[k] for k in 1:3 ]), reverse(eidx)...), 64)
    ψe = RandomState{Pure}(es, 4)
    ge = exp(-0.4im * (dag(Cup) ⊗ Cdn + dag(dag(Cup) ⊗ Cdn)))
    for (i, j) in [(1, 3), (3, 1), (1, 2)]
        h = ce(i, mu)' * ce(j, md) + ce(j, md)' * ce(i, mu)
        @test dense_e(apply(ge(i, j), ψe)) ≈ exp(-0.4im * h) * dense_e(ψe) atol = 1e-10
    end
    for (a, da) in singles, (b, db) in singles
        if isfermionic(a) == isfermionic(b)
            ref = [ value(da(i, n) * db(j, n)) for i in 1:n, j in 1:n ]
            @test expect2(ψ, (a, b)) ≈ ref atol = 1e-10
        end
    end
    # a tensor product on consecutive sites, of operators of one site and of renamed tensor
    # products of two, each with its width and its dense matrix from its first site
    for _ in 1:20
        fs = map(1:rand(2:3)) do _
            a, da = rand(singles)
            if rand(Bool)
                return (a, 1, da)
            end
            b, db = rand(singles)
            return (named(a ⊗ b, "T$(count[] += 1)"), 2, (j, n) -> da(j, n) * db(j + 1, n))
        end
        starts = cumsum([ 1; [ w for (_, w, _) in fs ] ])
        m = starts[end] - 1
        ref = prod(d(s, m) for ((_, _, d), s) in zip(fs, starts))
        @test matrix(reduce(⊗, first.(fs)), fe) ≈ ref atol = 1e-12
    end
end

@testset "Measurements that were refused" begin
    # a part that is all the system or nothing shares nothing with the rest
    ρ = mix(RandomState{Pure}(System(3, Qubit()), 2))
    @test mutual_info_renyi2(ρ, [1, 2, 3]) == 0
    @test mutual_info_renyi2(ρ, Int[]) == 0
    @test_throws "does not have" mutual_info_renyi2(ρ, [4])
    # a range stands for the equal vector
    st = State{Pure}(System(2, Qubit()), "Up")
    @test measure(st, [1:2]) == measure(st, [[1, 2]])
    # a reference conserving less than the measured state is refused by the message that says
    # what to do, weakening the state at every measurement costing a whole conversion of it
    sq = State{Pure}(System(2, Qubit(conserve = N)), "Up")
    @test_throws "conserves less than the measured state" measure(sq, Fidelity(State{Pure}(System(2, Qubit()), "Up")))
end

@testset "A cut at either end" begin
    # a cut at 0 gave 0 on a mixed state and an error about the entanglement entropy on a
    # pure one, and a negative cut was taken for an empty part on a mixed state
    ψ = RandomState{Pure}(System(3, Qubit()), 2)
    for st in (ψ, mix(ψ))
        @test mutual_info_renyi2(st, 0) == 0
        @test mutual_info_renyi2(st, 3) == 0
        @test_throws "the cut -1" mutual_info_renyi2(st, -1)
        @test_throws "the cut 4" mutual_info_renyi2(st, 4)
    end

    # SubRenyi2(3) was accepted and failed when measured: a link stands for the sites on its
    # left, as for MutualInfoRenyi2, and on a pure state is read off the spectrum
    for cut in 1:2
        @test renyi2(ψ, cut) ≈ renyi2(mix(ψ), collect(1:cut))
        @test renyi2(mix(ψ), cut) ≈ renyi2(mix(ψ), collect(1:cut))
    end
    @test renyi2(ψ, 0) == 0
    # no site has no entropy, as the cut 0, where partial_trace refused to keep nothing
    @test renyi2(ψ, Int[]) == 0
    @test renyi2(mix(ψ), Int[]) == 0
    @test last(only(measure(ψ, SubRenyi2(Int[])))) == 0
    # positions in a vector of another element type failed at the first measurement
    @test last(only(measure(ψ, SubRenyi2([])))) == 0
    @test renyi2(ψ, Any[1, 2]) ≈ renyi2(ψ, [1, 2])
    @test mutual_info_renyi2(ψ, Any[1]) ≈ mutual_info_renyi2(ψ, [1])
    @test mutual_info_renyi2(ψ, []) == 0
    @test renyi2(mix(ψ), 3) ≈ renyi2(mix(ψ))
    @test_throws "renyi2 was given the cut 4" renyi2(ψ, 4)
    @test SubRenyi2(2).name == "SubRenyi2(1,2)"
    @test last(only(measure(ψ, SubRenyi2(2)))) ≈ renyi2(mix(ψ), [1, 2])
end

@testset "Inner products and fidelities" begin
    sys = System(3, Qubit())
    up    = State{Pure}(sys, "Up")
    dn    = State{Pure}(sys, "Dn")
    plus  = State{Pure}(sys, "+")
    mixup = State{Pure}(sys, ["+", "+", "Up"])
    icplx = State{Pure}(sys, ["i", "+", "+"])

    # <+++|++Up> = 1 * 1 * <+|Up> = 1/sqrt(2)
    @test inner(plus, mixup) ≈ 1/√2
    @test dot(plus, mixup) == inner(plus, mixup)          # dot is an alias
    @test inner(plus, plus) ≈ 1
    @test inner(up, dn) ≈ 0 atol = 1e-14
    # the first argument is the one conjugated
    @test inner(icplx, plus) ≈ conj(inner(plus, icplx))
    @test imag(inner(icplx, plus)) ≉ 0                    # the test would be empty otherwise

    # fidelity is normalised, so neither the norm nor the trace of its arguments matters
    @test fidelity(plus, mixup) ≈ 0.5
    @test fidelity(plus, plus) ≈ 1
    @test fidelity(up, dn) ≈ 0 atol = 1e-14
    @test fidelity(3 * plus, 2 * mixup) ≈ fidelity(plus, mixup)

    # against a mixed representation, in either order. Mixing a pure state must not change
    # the answer, which is what makes the two methods one quantity
    fm = State{Mixed}(sys, "FullyMixed")
    @test fidelity(plus, mix(mixup)) ≈ fidelity(plus, mixup)
    @test fidelity(mix(mixup), plus) ≈ fidelity(plus, mixup)
    @test fidelity(up, fm) ≈ 1/8                          # <psi| I/2^3 |psi>
    @test fidelity(2 * up, fm) ≈ 1/8

    # the normalised Hilbert-Schmidt overlap of two mixed states. On two mixed pure states
    # it coincides with the fidelity, the traces of the squares being one
    @test hs_fidelity(mix(plus), mix(plus)) ≈ 1
    @test hs_fidelity(mix(plus), mix(mixup)) ≈ fidelity(plus, mixup)
    @test hs_fidelity(mix(up), fm) ≈ √2/4
    @test hs_fidelity(2 * mix(up), 3 * fm) ≈ hs_fidelity(mix(up), fm)

    # a pure and a mixed representation are vectors of different spaces
    @test_throws "pure and a mixed representation" inner(plus, fm)
    @test_throws "pure and a mixed representation" inner(fm, plus)

    # and two mixed ones have no fidelity to speak of, the Uhlmann one needing a spectrum.
    # The refusal names what to reach for instead
    @test_throws "no fidelity between two mixed" fidelity(fm, mix(plus))
    @test_throws "hs_fidelity" fidelity(fm, mix(plus))
end

@testset "Fidelity and Overlap as measurements" begin
    # the reference is written when the simulation is described, so it lives on a system
    # of its own; the state functions are the ones putting it where it can be contracted
    ref = State{Pure}(System(3, Qubit()), "+")
    st = State{Pure}(System(3, Qubit()), ["+", "+", "Up"])
    @test first(only(measure(st, Fidelity(ref)))) == "Fidelity"
    @test last(only(measure(st, Fidelity(ref)))) ≈ 0.5
    @test last(only(measure(st, Overlap(ref)))) ≈ 1/√2
    # and through a run, where the system does not exist until the phase creates it
    sim = runTMS(SimData(phases = [
            CreateState{Pure}(3, Qubit(), ["+", "+", "Up"];
                final_measurements = Data("d") => [Fidelity(ref), Overlap(ref)])]);
        output = devnull)
    @test only(sim.data["d"]["Fidelity"]["data"]) ≈ 0.5
    @test only(sim.data["d"]["Overlap"]["data"]) ≈ 1/√2

    # a reference built before the state was weakened is weakened in its turn to what the
    # state conserves, on a system mixing sites that conserve and sites that do not as well
    for sites in ([Fermion(conserve = strong(N)) for _ in 1:3],
                  [Fermion(conserve = strong(N)), Qubit(), Fermion(conserve = strong(N))])
        p = State{Pure}(System(sites), ["Occ", sites[2] isa Qubit ? "Up" : "Emp", "Occ"])
        @test last(only(measure(weaken(mix(p)), Fidelity(p)))) ≈ 1
        @test last(only(measure(weaken(mix(p), ()), Fidelity(p)))) ≈ 1
        @test last(only(measure(weaken(p), Overlap(p)))) ≈ 1
    end
    # but a reference conserving less than the state is refused, weakening the state at every
    # measurement costing a whole conversion of it
    strong_state = State{Pure}(System(2, Fermion(conserve = strong(N))), "Occ")
    weak_ref = State{Pure}(System(2, Fermion(conserve = N)), "Occ")
    @test_throws "conserves less than the measured state" measure(mix(strong_state), Fidelity(weak_ref))
end

@testset "Real, imaginary and complex values" begin
    # the kind of a measurement is decided on the measurement, never on a value: an operator
    # is real when `dag` gives it back, imaginary when it gives its opposite, complex when
    # the symbolic test proves neither
    rk, ik, ck = TensorMixedStates.real_kind, TensorMixedStates.imaginary_kind,
                 TensorMixedStates.complex_kind
    kind(o) = TensorMixedStates.make_obs(o).kind
    hop = dag(C)(1) * C(2)
    @test kind(X(1)) == rk
    @test kind(Z(1)Z(2) + X(1)) == rk
    @test kind(Proj("Up")(1)) == rk
    @test kind(hop + dag(C)(2) * C(1)) == rk
    @test kind(im * X(1) * Y(2)) == ik
    @test kind(hop - dag(C)(2) * C(1)) == ik
    @test kind(Sp(1)) == ck
    @test kind(hop) == ck
    # a miss of the test comes out complex, never real: `Sm` is a name whose relation to
    # `Sp` it does not know
    @test kind(Sp(1)Sm(2) + Sm(1)Sp(2)) == ck
    @test kind(X) == rk
    @test kind(Sp) == ck
    # a correlation matrix takes one kind, complex when its entries differ: the diagonal of
    # (X, Y) is <XY> = i<Z> and the rest is real
    @test kind((X, X)) == rk
    @test kind((X, Y)) == ck
    @test kind((Sp, Sm)) == ck
    @test kind((dag(C), C)) == ck
    @test kind((N, N)) == rk
    # a function is real unless declared otherwise, a number takes the kind of its type
    @test kind(Purity) == rk
    @test kind(t -> exp(im * t)) == rk
    @test kind(0.5) == rk
    @test kind(0.5im) == ck
    @test kind(Overlap(State{Pure}(System(2, Qubit()), "+"))) == ck

    # whatever the test says real or imaginary is: on states with complex amplitudes the
    # part it drops is rounding
    for (st, ops) in ((RandomState{Pure}(System(3, Qubit()), 4),
                       [X(1), Z(1)Z(2) + X(3), Proj("Dn")(2), exp(0.3X)(1), im * X(1) * Y(2),
                        (im * Z)(3), im * (X(1)Y(2) - Y(1)X(2))]),
                      (RandomState{Pure}(System(3, Fermion()), 4),
                       [hop + dag(C)(2) * C(1), hop - dag(C)(2) * C(1),
                        im * (dag(C)(1) * C(3) - dag(C)(3) * C(1)),
                        (C + dag(C))(1) * (C + dag(C))(2)]))
        for o in ops
            k = kind(o)
            v = expect(st, o)
            @test k ≠ ck
            @test abs(k == rk ? imag(v) : real(v)) < 1e-12
        end
    end

    # measure gives each value the kind of its measurement whatever the element type of the
    # state: a real value is a Float64 on a complex state, a complex one a ComplexF64 on a
    # real state
    value(st, m) = last(only(measure(st, m)))
    up = State{Pure}(System(2, Qubit()), "Up")
    ph = apply(Phase(0.7)(1), State{Pure}(System(2, Qubit()), "+"))
    plus_i = State{Pure}(System(2, Qubit()), ["+", "i"])
    @test value(ph, X(1)) isa Float64
    @test value(ph, X(1)) ≈ cos(0.7)
    @test value(up, Sp(1)) isa ComplexF64
    @test value(ph, Sp(1)) ≈ exp(0.7im) / 2
    @test value(up, Overlap(up)) isa ComplexF64
    # an imaginary value is its imaginary part, under a name saying so
    iv = only(measure(plus_i, im * X(1) * Y(2)))
    @test first(iv) == "Im(im*X(1)*Y(2))"
    @test last(iv) ≈ 1
    # the diagonal of (X, Y), which used to be written as zero
    @test value(up, (X, Y)) isa Matrix{ComplexF64}
    @test value(up, (X, Y)) ≈ [im 0; 0 im]

    # a declaration overrides the kind found by the test, holds for each part of a vector,
    # and the innermost one holds
    xy = Sp(1)Sm(2) + Sm(1)Sp(2)
    @test value(ph, xy) isa ComplexF64
    @test value(ph, RealValue(xy)) isa Float64
    @test value(ph, ComplexValue(X(1))) isa ComplexF64
    @test value(ph, RealValue(ComplexValue(X(1)))) isa ComplexF64
    vs = @test_logs (:warn, r"large imaginary part: time 0.0, Sp\(1\)") measure(
        ph, RealValue([Sp(1), ComplexValue(Sp(2))]))
    @test first.(vs) == ["Sp(1)", "Sp(2)"]
    @test last(vs[1]) ≈ cos(0.7) / 2
    @test last(vs[2]) isa ComplexF64
    iv = @test_logs (:warn, r"large real part: time 0.0, X\(1\)") only(
        measure(ph, ImaginaryValue(X(1))))
    @test first(iv) == "Im(X(1))"
    # a declaration holds for each element of any array, as a measurement without one is
    # taken element by element: a range or a view was made a single measurement, which no
    # method could compute
    @test measure(ph, RealValue(1:3)) == measure(ph, 1:3)
    vs = measure(ph, ComplexValue(view([X(1), X(2), Z(1)], 1:2)))
    @test first.(vs) == ["X(1)", "X(2)"]
    @test all(v -> last(v) isa ComplexF64, vs)
    @test value(ph, Check("c", RealValue(1:3), [1, 2, 3], 1e-10))[3] == 0

    # the part dropped is reported relative to the modulus, which a zero real part leaves
    # finite, and a function is told how to keep it; a complex value drops nothing
    @test_logs (:warn, r"rel  1.0e\+00") measure(plus_i, RealValue(im * X(1) * Y(2)))
    @test_logs (:warn, r"sp .*ComplexValue keeps it") measure(ph,
        StateFunc("sp", st -> expect(st, Sp(1))))
    @test_logs measure(ph, Sp(1))

    # a Check compares the values as computed: a complex reference given as a function of
    # time, real unless declared otherwise, is written as its real part but compared whole
    c = @test_logs (:warn, r"func .*ComplexValue keeps it") value(ph,
        Check("c", Sp(1), _ -> exp(0.7im) / 2, 1e-10))
    @test c[1] ≈ exp(0.7im) / 2
    @test c[2] ≈ cos(0.7) / 2
    @test c[3] < 1e-10
    # each part keeps its own type, rather than all three becoming complex
    @test c[2] isa Float64
    @test c[3] isa Float64
    # a constant reference keeps its kind, and an imaginary part stays a complex number on
    # the line of a check, where no name could say what it is
    c = @test_logs value(ph, Check("c", Sp(1), exp(0.7im) / 2, 1e-10))
    @test c[2] isa ComplexF64
    c = value(plus_i, Check("c", im * X(1) * Y(2), 1im, 1e-10))
    @test c[1] isa ComplexF64
    @test c[1] ≈ 1im
    c = value(State{Pure}(System(2, Qubit()), "+"), Check("c", Sp, [0.5, 0.5], 1e-10))
    @test c[1] isa Vector{ComplexF64}
    @test c[2] == [0.5, 0.5]

    # a vector in a set stands for its measurements, each under its own name, and the names
    # are then checked for duplicates
    @test first.(measure(ph, [[X(1), Z(2)]])) == ["X(1)", "Z(2)"]
    @test_throws "named \"X(1)\"" Measure([X(1), [X(1)]])
end

@testset "Energy variance" begin
    sys = System(6, Qubit())
    ising = -sum(Z(i) * Z(i+1) for i in 1:5)
    field = -sum(X(i) for i in 1:6)
    h = ising + field

    # zero exactly on an eigenstate, whichever one
    @test variance(State{Pure}(sys, "Up"), ising) ≈ 0 atol = 1e-12
    @test variance(State{Pure}(sys, ["Up", "Dn", "Up", "Dn", "Up", "Dn"]), ising) ≈ 0 atol = 1e-12
    @test variance(State{Pure}(sys, "+"), field) ≈ 0 atol = 1e-12
    # |+...+> is an eigenstate of the field but not of the couplings. <ZZ> is zero there and
    # so are the cross terms, so the variance is the number of couplings
    @test variance(State{Pure}(sys, "+"), h) ≈ 5

    # against the naive route, the one that forms H^2 and that the implementation avoids
    r = RandomState{Pure}(sys, 8)
    @test variance(r, h) ≈ real(expect(r, h * h)) - real(expect(r, h))^2
    # neither the norm nor a phase of the state changes it
    @test variance(3r, h) ≈ variance(r, h)
    @test variance(im * r, h) ≈ variance(r, h)
    # an mpo already built gives the same answer
    @test variance(r, make_mpo(r, h)) ≈ variance(r, h)

    # what it is for: a converged ground state has a variance the energy alone cannot show
    e, gs = dmrg(h, RandomState{Pure}(sys, 16); nsweeps = 12,
                 limits = Limits(cutoff = 1e-14, maxdim = 32))
    # and its energy is the lowest eigenvalue of the dense hamiltonian
    LA = TensorMixedStates.LinearAlgebra
    on(m, i) = kron([ k == i ? m : [1. 0. ; 0. 1.] for k in 1:6 ]...)
    dense = -sum(on([1. 0. ; 0. -1.], i) * on([1. 0. ; 0. -1.], i + 1) for i in 1:5) -
            sum(on([0. 1. ; 1. 0.], i) for i in 1:6)
    @test e ≈ minimum(LA.eigvals(LA.Symmetric(dense))) atol = 1e-10
    @test variance(gs, h) < 1e-8
    @test variance(gs, h) < variance(r, h)

    # as a measurement
    @test first(only(measure(gs, Variance(h)))) == "Variance"
    @test last(only(measure(gs, Variance(h)))) ≈ variance(gs, h)

    # a hamiltonian on a density matrix has no variance in this sense
    @test_throws "needs a pure representation" variance(mix(r), h)
end

@testset "Observables do not depend on the charges" begin
    # the same physics twice, on sites that declare a conserved quantity and on sites that
    # do not. Declaring one splits the tensors into blocks, it does not change what is
    # measured, so every number below has to come out the same both times. This says
    # nothing about which values are right, which the rest of this file covers, and
    # everything about the charges leaving them alone
    function chain(site)
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        # three configurations of the same particle number, so their sum has a charge
        s = p(["Occ", "Emp", "Occ", "Emp"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"]) -
            0.3 * p(["Occ", "Occ", "Emp", "Emp"])
        return s / norm(s)
    end
    d, q = chain(Fermion()), chain(Fermion(conserve = N))
    # compared with an absolute tolerance, because several of these are exactly zero: the
    # Renyi 2 entropy of a pure state seen as mixed, and the trace of a dissipator applied
    # to a state, a dissipator being traceless. `≈` alone has no tolerance against zero, so
    # the two paths rounding to 0.0 and 4e-16 would count as a difference, which they are not
    both(f) = @test isapprox(f(d), f(q); atol = 1e-12)

    both(s -> expect(s, N(2)))
    both(s -> expect2(s, (N, N)))
    both(s -> [ expect(s, dag(C)(i) * C(j)) for i in 1:4, j in 1:4 ])
    both(s -> norm(s))
    both(s -> entanglement_entropy(s, 2)[1])
    both(s -> collect(entanglement_entropy(s, 2)[2]))

    both(s -> trace(mix(s)))
    both(s -> trace2(mix(s)))
    both(s -> renyi2(mix(s)))
    both(s -> hermiticity(mix(s)))
    both(s -> expect(mix(s), N(2)))
    both(s -> trace(partial_trace(mix(s), [1, 3])))
    both(s -> mutual_info_renyi2(mix(s), 2))
    both(s -> trace(apply(Gate(F)(1), mix(s))))
    both(s -> trace(apply(Dissipator(N)(2), mix(s))))
    both(s -> expect(apply(SetState("Occ")(2), mix(s)), N(2)))
end

@testset "A mixed state built without going through mix" begin
    # its local matrix used to be written straight onto the mixed index. That index is a
    # combination, and combining charged indices merges and sorts their sectors, so the
    # flat order of its basis is not the order of the matrix: the diagonal landed off the
    # diagonal, where the charges annihilate it, and the state came out with a trace of
    # zero and no complaint
    direct(site) = State{Mixed}(System(3, site), ["2", "0", "1"])
    bd, bq = direct(Boson(4)), direct(Boson(4, conserve = N))
    @test trace(bd) ≈ 1
    @test trace(bq) ≈ 1
    @test [ expect(bq, N(k)) for k in 1:3 ] ≈ [2, 0, 1]
    @test [ expect(bq, N(k)) for k in 1:3 ] ≈ [ expect(bd, N(k)) for k in 1:3 ]

    # a density matrix given as a matrix goes the same way, and a diagonal one is a
    # mixture over sectors, which the charges allow
    rho = [0.1 0. 0. 0.; 0. 0.2 0. 0.; 0. 0. 0.3 0.; 0. 0. 0. 0.4]
    md = State{Mixed}(System(2, Boson(4)), rho)
    mq = State{Mixed}(System(2, Boson(4, conserve = N)), rho)
    @test trace(mq) ≈ 1
    @test expect(mq, N(1)) ≈ 2
    @test expect(mq, N(1)) ≈ expect(md, N(1))
end

@testset "Observables do not depend on the kind of symmetry" begin
    strong = TensorMixedStates.strong
    # the same physics conserved weakly and strongly. The two label the tensors differently,
    # one keeping the difference of the ket and bra charges and the other keeping them
    # apart, and neither changes a measured number. Under a strong symmetry the trace stops
    # being a product of one vector per site, so a strong state is measured through its weak
    # form, and this covers that detour
    function chain(site)
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        s = p(["Occ", "Emp", "Occ", "Emp"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"])
        return mix(s / norm(s))
    end
    w, s = chain(Fermion(conserve = N)), chain(Fermion(conserve = strong(N)))
    both(f) = @test isapprox(f(w), f(s); atol = 1e-12)

    @test trace(s) ≈ 1
    both(trace)
    both(trace2)
    both(x -> expect(x, N(2)))
    both(x -> expect1(x, N))
    both(hermiticity)
    # a correlation whose two ends do not conserve the charge is measured all the same, the
    # two shifts it brings cancelling
    both(x -> [ expect(x, dag(C)(i) * C(j)) for i in 1:4, j in 1:4 ])
    both(x -> trace(apply(Dissipator(N)(2), x)))

    # the adjoint exchanges ket and bra, which swaps the two charges a strong symmetry keeps
    # apart. It relabels the charges so that the exchange becomes a permutation of zero flux
    @test hermiticity(s) ≈ 1
    @test trace(hermitianize(s)) ≈ 1
    @test flux(dag(s).state) == flux(s.state)

    # on a Hermitian state the adjoint is invisible, so it is checked on one that is not,
    # holding a coherence between sites 1 and 3 that an off diagonal observable reads
    function skew(site)
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        c1, c2 = ["Occ", "Emp", "Emp", "Occ"], ["Emp", "Emp", "Occ", "Occ"]
        a, b = p(c1) + 0.7im * p(c2), p(c2) - 0.4 * p(c1)
        return State(mix(a), mix(a).state + 0.6im * mix(b).state)
    end
    ws, ss = skew(Fermion(conserve = N)), skew(Fermion(conserve = strong(N)))
    adjoint_numbers(x) = [ expect(dag(x), dag(C)(1) * C(3)), expect(dag(x), N(3)),
                           inner(dag(x), x) ]
    @test isapprox(adjoint_numbers(ws), adjoint_numbers(ss); atol = 1e-12)
    @test abs(expect(dag(ss), dag(C)(1) * C(3))) > 0.1
    @test expect(dag(ss), dag(C)(1) * C(3)) ≈ conj(expect(ss, dag(C)(3) * C(1)))
    @test hermiticity(ss) < 0.99

    # tracing part of the sites out is the one thing that cannot be done: what is left is a
    # mixture over several sectors, and a state keeping the two charges apart has only one.
    # It is refused rather than weakened behind the user's back, who is to see the strong
    # symmetry go
    @test_ok partial_trace(w, [1, 3])
    @test_throws "weaken it first" partial_trace(s, [1, 3])
    @test_ok partial_trace(weaken(s), [1, 3])

    # a number read off a partial trace is measured all the same, through the weak form
    both(x -> renyi2(x, [1, 2]))
    both(x -> mutual_info_renyi2(x, [1, 2]))
    both(x -> mutual_info_renyi2(x, 2))
end

@testset "A placed operator where a generic one is expected" begin
    # each raised a MethodError naming no way out, expect1 one on `length`
    ψ = State{Pure}(System(2, Qubit()), "Up")
    @test_throws "measure X(1) with expect" expect1(ψ, X(1))
    @test_throws "measure X(1) with expect" expect1(ψ, [Z, X(1)])
    @test_throws "(X, Y) rather than (X(1), Y(2))" expect2(ψ, (X(1), X(2)))
    @test_throws "(X, Y) rather than (X(1), Y(2))" expect2(ψ, [(X, X), (X(1), X(2))])
    @test_throws "rather than Dissipator(Sm(1))" Dissipator(Sm(1))
end
