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
end

@testset "Fermionic correlations" begin
    # dag(C) must be recognized as fermionic, and a product of two
    # fermionic operators must not be
    @test isfermionic(C)
    @test isfermionic(dag(C))
    @test !isfermionic(dag(C) * C)
    @test !isfermionic(N)
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
    # the mixed representation must agree with the pure one
    @test expect2(mix(st), (dag(C), C)) ≈ ref atol=1e-8
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

@testset "Entanglement and entropies" begin
    L2 = log(2)
    s2 = System(2, Qubit())
    s4 = System(4, Qubit())
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
    # a non contiguous part has no such shortcut, but must still be accepted on a pure state
    @test mutual_info_renyi2(spread, [1, 3]) ≈ mutual_info_renyi2(mix(spread), [1, 3])

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
    # the deprecated spelling forwards to the new one, label included
    @test MutualInfoRenyi2(2).name == "MutualInfoRenyi2(2)"
    @test Mutual_Info_Renyi2(2).name == MutualInfoRenyi2(2).name
end

@testset "Entanglement resolved by sector" begin
    Q = TensorMixedStates.ITensors
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
    @test sort(collect(keys(d)); by = string) == [Q.QN("N", 0), Q.QN("N", 1)]
    @test d[Q.QN("N", 1)].weight ≈ 0.8
    @test d[Q.QN("N", 0)].weight ≈ 0.2
    # cut after two, it is one particle in both terms, and the entanglement lives inside it
    d = entanglement_by_sector(st, 2)
    @test collect(keys(d)) == [Q.QN("N", 1)]
    @test d[Q.QN("N", 1)].weight ≈ 1
    @test d[Q.QN("N", 1)].spectrum ≈ [0.8, 0.2]
    @test d[Q.QN("N", 1)].entropy ≈ first(entanglement_entropy(st, 2))

    # several quantities at once, and a charge modulo 2 over a basis whose charges repeat
    e = System(3, Electron(conserve = (Ntot, 2Sz)))
    q(v) = State{Pure}(e, v)
    d = entanglement_by_sector((q(["Up", "Dn", "Emp"]) + 0.5 * q(["UpDn", "Emp", "Emp"])) / sqrt(1.25), 1)
    @test d[Q.QN(("Ntot", 1), ("2Sz", 1))].weight ≈ 0.8
    @test d[Q.QN(("Ntot", 2), ("2Sz", 0))].weight ≈ 0.2
    b = System(3, Boson(4, conserve = parity(N)))
    r(v) = State{Pure}(b, v)
    d = entanglement_by_sector((r(["1", "2", "0"]) + 0.5 * r(["0", "3", "0"])) / sqrt(1.25), 1)
    @test d[Q.QN("parity(N)", 1, 2)].weight ≈ 0.8
    @test d[Q.QN("parity(N)", 0, 2)].weight ≈ 0.2

    # without charges there is a single sector holding the whole spectrum
    fd = System(4, Fermion())
    pd(v) = State{Pure}(fd, v)
    sd = (pd(["Occ", "Emp", "Emp", "Occ"]) + 0.5 * pd(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
    d = entanglement_by_sector(sd, 2)
    @test collect(keys(d)) == [Q.QN()]
    @test d[Q.QN()].spectrum ≈ last(entanglement_entropy(sd, 2))

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
    # orthogonality centre pushed to the right.
    sys = System(6, Qubit())
    st = apply(Sp(1) * Sp(3), State{Pure}(sys, "+"))
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
    @test length(partial_trace(stm, 1:3; keepers = true)) == 3
    # `measure` hands back a vector of name => value pairs, one per measurement
    @test last(only(measure(stm, SubRenyi2(1:3)))) ≈ renyi2(stm, [1, 2, 3])
    @test last(only(measure(stm, MutualInfoRenyi2(1:3)))) ≈ mutual_info_renyi2(stm, [1, 2, 3])
    # a partial trace needs a density matrix, and says so rather than raising a MethodError
    @test_throws "mixed representation" partial_trace(stp, [1, 2])
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
                final_measures = Data("d") => [Fidelity(ref), Overlap(ref)])]);
        output = devnull)
    @test only(sim.data["d"]["Fidelity"]["data"]) ≈ 0.5
    @test only(sim.data["d"]["Overlap"]["data"]) ≈ 1/√2
end

@testset "Energy variance" begin
    sys = System(6, Qubit())
    ising = -sum(Z(i) * Z(i+1) for i in 1:5)
    field = -sum(X(i) for i in 1:6)
    h = ising + field

    # zero exactly on an eigenstate, whichever one
    @test variance(ising, State{Pure}(sys, "Up")) ≈ 0 atol = 1e-12
    @test variance(ising, State{Pure}(sys, ["Up", "Dn", "Up", "Dn", "Up", "Dn"])) ≈ 0 atol = 1e-12
    @test variance(field, State{Pure}(sys, "+")) ≈ 0 atol = 1e-12
    # |+...+> is an eigenstate of the field but not of the couplings. <ZZ> is zero there and
    # so are the cross terms, so the variance is the number of couplings
    @test variance(h, State{Pure}(sys, "+")) ≈ 5

    # against the naive route, the one that forms H^2 and that the implementation avoids
    r = RandomState{Pure}(sys, 8)
    @test variance(h, r) ≈ real(expect(r, h * h)) - real(expect(r, h))^2
    # neither the norm nor a phase of the state changes it
    @test variance(h, 3r) ≈ variance(h, r)
    @test variance(h, im * r) ≈ variance(h, r)
    # an mpo already built gives the same answer
    @test variance(make_mpo(r, h), r) ≈ variance(h, r)

    # what it is for: a converged ground state has a variance the energy alone cannot show
    _, gs = dmrg(h, RandomState{Pure}(sys, 16); nsweeps = 12,
                 limits = Limits(cutoff = 1e-14, maxdim = 32))
    @test variance(h, gs) < 1e-8
    @test variance(h, gs) < variance(h, r)

    # as a measurement
    @test first(only(measure(gs, Variance(h)))) == "Variance"
    @test last(only(measure(gs, Variance(h)))) ≈ variance(h, gs)

    # a hamiltonian on a density matrix has no variance in this sense
    @test_throws "needs a pure representation" variance(h, mix(r))
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
    both(s -> trace(apply(Dissipator(C)(2), mix(s))))
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
