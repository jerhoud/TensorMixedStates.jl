# Functions that extract numbers from a state, checked against exact values.
#
# Goes here: sampling, one and two point correlations, entanglement and entropy
# measures. These are the functions with an analytic answer on simple states, so tests
# here should compare against that answer rather than against a recorded output.

@testset "Qubit sampling" begin
    sys = System(3, Qubit())
    st = State{Pure}(sys, ["Up", "Dn", "+"])  # site 3 is a 50/50 superposition

    # deterministic sites always give the same outcome
    for _ in 1:20
        @test sample(st, 1) == 0
        @test sample(st, 2) == 1
    end

    # superposed site: frequency should be close to 1/2
    n = 2000
    s3 = [sample(st, 3) for _ in 1:n]
    @test all(x -> x in (0, 1), s3)
    @test isapprox(sum(s3) / n, 0.5; atol = 0.05)

    # sampling the whole state is consistent with per-site sampling
    for _ in 1:20
        r = sample(st)
        @test r[1] == 0 && r[2] == 1 && r[3] in (0, 1)
    end

    # classical diagonal mixture with known probabilities
    p0 = 0.2
    stm = State{Mixed}(System(1, Qubit()), [p0 0. ; 0. 1 - p0])

    nm = 4000
    samples = [sample(stm, 1) for _ in 1:nm]
    @test all(x -> x in (0, 1), samples)
    @test isapprox(count(==(0), samples) / nm, p0; atol = 0.03)

    samples_full = [sample(stm)[1] for _ in 1:nm]
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
    # mixing fermionic and non fermionic pairs in one call must not disturb either
    m = expect2(st, [(N, N), (dag(C), C)])
    @test m[1] ≈ expect2(st, (N, N)) atol=1e-10
    @test m[2] ≈ ref atol=1e-10
    # the mixed representation must agree with the pure one
    @test expect2(mix(st), (dag(C), C)) ≈ ref atol=1e-8
    # expect works on indexed operators, but only once they went through simplify,
    # which is what inserts the Jordan-Wigner strings (Multi_F for two sites or more)
    stm = mix(st)
    for d in 1:n - 1
        op = simplify(dag(C)(1) * C(1 + d))
        @test expect(st, op) ≈ ref[1, 1 + d] atol=1e-10
        @test expect(stm, op) ≈ ref[1, 1 + d] atol=1e-8
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
