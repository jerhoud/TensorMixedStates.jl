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
    @test m(ghz, Linkdim) == 2
    @test m(ghz, EE(2)) ≈ L2
    @test m(fullymixed, Purity) ≈ 1 / 16
    @test m(fullymixed, Renyi2) ≈ 4L2
    @test m(mix(ghz), SubRenyi2([1])) ≈ L2
    @test m(ghz, Mutual_Info_Renyi2(2)) ≈ 2L2
    @test m(mix(ghz), Mutual_Info_Renyi2(2)) ≈ 2L2
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
