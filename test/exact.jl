# The exact representation of small systems: an ExactState created from the descriptions of State
# and from a State, its measurements against those of the State, its conversions and its file.
# The terms acting on it are tested in operators.jl and fermions.jl.

@testset "Creation from descriptions" begin
    q = System(3, Qubit())
    for st in ["Up", ["+", "0", "i"]]
        @test TMS.ExactState{Pure}(q, st).state ≈ Vector(State{Pure}(q, st))
    end
    ρ = State{Mixed}(q, ["Up", "+", "Dn"])
    @test TMS.ExactState{Mixed}(q, ["Up", "+", "Dn"]).state ≈ ket_bra_matrix(ρ)
    mix = TMS.ExactState{Mixed}(q, Mixture([0.25 => "Up", 0.75 => "Dn"]))
    @test mix.state ≈ LA.Diagonal([0.25; zeros(6); 0.75])
    # the projector on the 6 states of two fermions out of four, divided by 6
    f = System(4, Fermion())
    @test TMS.ExactState{Mixed}(f, FullyMixed(N => 2)).state ≈
          LA.Diagonal([ count_ones(a - 1) == 2 ? 1 / 6 : 0 for a in 1:16 ])
    @test_throws "not pure" TMS.ExactState{Pure}(f, FullyMixed(N => 2))
end

# `proj`: Hermitian operators measured projectively, a product of involutions, a sum of operators
# of one site, and one diagonalised on its sites, past a site carrying a string
@testset "Measurements on $name" for (name, sys, ops, proj) in [
        ("qubits", System(4, Qubit()), [X(1) * Z(3) + 2 * Y(2), Sp(1) * Sm(4)],
         [Z(1) * X(3), Z(1) + Z(2) + Z(4), X(1) * X(3) + Z(1)]),
        ("fermions", System(4, Fermion()), [dag(C)(1) * C(3) + dag(C)(3) * C(1), C(2), N(4)],
         [F(1) * F(3), N(1) + N(2) + N(4), dag(C)(1) * C(3) + dag(C)(3) * C(1)]),
        ("a qubit among fermions", System([Fermion(), Qubit(), Electron(), Fermion()]),
         [dag(C)(1) * Cup(3) + dag(Cup)(3) * C(1), C(4), X(2) * Ntot(3)],
         [X(2) * F(4), N(1) + Ntot(3) + N(4), dag(C)(1) * C(4) + dag(C)(4) * C(1)])]
    rng = Xoshiro(20261103)
    ψ, ρ = random_state(Pure, sys, 4; rng), random_state(Mixed, sys, 4; rng)
    eψ, eρ = TMS.ExactState(ψ), TMS.ExactState(ρ)
    @test eψ.state ≈ raw_vector(ψ)
    @test eρ.state ≈ ket_bra_matrix(ρ)
    for (e, st) in [(eψ, ψ), (eρ, ρ)]
        for o in ops
            @test expect(e, o; complex = true) ≈ expect(st, o; complex = true) atol = 1e-12
            @test measure(e, o + dag(o)) ≈ measure(st, o + dag(o)) atol = 1e-12
            # through expect, a non-Hermitian operator too
            @test variance(e, o) ≈ variance(st, o) atol = 1e-12
            @test measure(e, Variance(o + dag(o))) ≈ variance(st, o + dag(o)) atol = 1e-12
        end
        @test expect1(e, F) ≈ expect1(st, F) atol = 1e-12
        @test trace(e) ≈ trace(st)
        @test norm(e) ≈ norm(st)
        @test purity(e) ≈ purity(st)
        # the Schmidt spectrum, the operator space one on a mixed state, of all the singular
        # values where the MPS has those of its link
        for cut in 1:4
            @test entanglement_entropy(e, cut) ≈ entanglement_entropy(st, cut) atol = 1e-10
            sp, ref = entanglement_spectrum(e, cut), entanglement_spectrum(st, cut)
            k = min(length(sp), length(ref))
            @test sp[1:k] ≈ ref[1:k] atol = 1e-12
            @test sum(sp[k+1:end]; init = 0.0) + sum(ref[k+1:end]; init = 0.0) < 1e-12
        end
        @test measure(e, EntanglementEntropy(2)) ≈ entanglement_entropy(st, 2) atol = 1e-10
        # the results, their probabilities and the states each leaves
        for o in [proj; Level(3)]
            ref = probabilities(st, o)
            @test first.(probabilities(e, o)) ≈ first.(ref)
            @test last.(probabilities(e, o)) ≈ last.(ref) atol = 1e-10
            for (λ, p) in ref
                if p > 1e-6
                    c, r = collapse(e, o => λ), collapse(st, o => λ)
                    @test c.probability ≈ r.probability atol = 1e-10
                    @test c.state.state ≈ TMS.ExactState(r.state).state atol = 1e-10
                end
            end
        end
        # the signs of the fermionic sites traced out, on the left of those kept or not
        for kept in [[1], [4], [2, 3], [1, 3], [2, 4], [1, 2, 4], 1:4]
            @test reduced_density_matrix(e, kept) ≈ reduced_density_matrix(st, kept) atol = 1e-12
        end
        @test reduced_density_matrix(e, Any[1, 3]) ≈ reduced_density_matrix(st, [1, 3]) atol =
              1e-12
        @test von_neumann_entropy(e, [1, 3]) ≈ von_neumann_entropy(st, [1, 3]) atol = 1e-10
        @test renyi2(e, [2, 4]) ≈ renyi2(st, [2, 4]) atol = 1e-10
        @test mutual_info_renyi2(e, [1, 3]) ≈ mutual_info_renyi2(st, [1, 3]) atol = 1e-10
        # a part without fermionic sites, as log_negativity asks
        bosonic = [ i for i in 1:4 if TMS.next_fermion(sys, i) ≠ i ]
        if !isempty(bosonic)
            a, b = [first(bosonic) == 1 ? 3 : 1], [first(bosonic)]
            @test log_negativity(e, a, b) ≈ log_negativity(st, a, b) atol = 1e-10
        end
    end
    @test mix(eψ).state ≈ ket_bra_matrix(mix(ψ))
    @test mix(eρ) === eρ
    # draws following the probabilities of the product states, of a site, of an operator
    d = size(eρ.state, 1)
    for (e, w) in [(eψ, abs2.(eψ.state)), (eρ, real(LA.diag(eρ.state)))]
        p = w / sum(w)
        draws = [ sample(e; rng) for _ in 1:4000 ]
        dims = [ TMS.dim(s) for s in sys.sites ]
        index(x) = 1 + sum(x[i] * prod(dims[i+1:end]; init = 1) for i in 1:4)
        counts = [ count(x -> index(x) == k, draws) for k in 1:d ]
        @test all(abs.(counts / 4000 - p) .≤ 5 * sqrt.(p .* (1 .- p) / 4000) .+ 1e-12)
        q = last.(probabilities(e, 2))
        counts = [ count(==(x), [ sample(e, 2; rng) for _ in 1:4000 ]) for x in 0:length(q)-1 ]
        @test all(abs.(counts / 4000 - q) .≤ 5 * sqrt.(q .* (1 .- q) / 4000) .+ 1e-12)
    end
    # the product state of the second level of site 2 and the first of the others
    levels = [ [ k == (i == 2 ? 2 : 1) ? 1.0 : 0.0 for k in 1:TMS.dim(sys[i]) ] for i in 1:4 ]
    @test sample(TMS.ExactState{Pure}(sys, levels)) == [0, 1, 0, 0]
    @test only(entanglement_by_sector(eψ, 2))[2].entropy ≈ entanglement_entropy(ψ, 2) atol = 1e-10
    @test_throws "takes a pure state" entanglement_by_sector(eρ, 2)
    # between two states, and a density matrix that is not Hermitian
    ψ2, ρ2 = random_state(Pure, sys, 4; rng), random_state(Mixed, sys, 4; rng)
    eψ2, eρ2 = TMS.ExactState(ψ2), TMS.ExactState(ρ2)
    @test inner(eψ, eψ2) ≈ inner(ψ, ψ2)
    @test LA.dot(eρ, eρ2) ≈ inner(ρ, ρ2)
    @test fidelity(eψ, eψ2) ≈ fidelity(ψ, ψ2)
    @test fidelity(eψ, eρ) ≈ fidelity(ψ, ρ)
    @test fidelity(eρ, eψ) ≈ fidelity(ψ, ρ)
    @test hs_fidelity(eρ, eρ2) ≈ hs_fidelity(ρ, ρ2)
    # that of Uhlmann on two mixed ones, which is the pure one on pure states, symmetric, 1 for
    # one state and the overlap of the probabilities of two commuting ones
    @test fidelity(mix(eψ), mix(eψ2)) ≈ fidelity(ψ, ψ2) rtol = 1e-10
    @test fidelity(mix(eψ), TMS.ExactState(3 * ρ)) ≈ fidelity(ψ, ρ) rtol = 1e-10
    @test fidelity(eρ, eρ2) ≈ fidelity(eρ2, eρ)
    @test fidelity(eρ, eρ) ≈ 1
    p, q = rand(rng, size(eρ.state, 1)), rand(rng, size(eρ.state, 1))
    @test fidelity(TMS.ExactState(sys, LA.diagm(p)), TMS.ExactState(sys, LA.diagm(q))) ≈
          sum(sqrt.(p .* q))^2 / (sum(p) * sum(q))
    nh = ρ + im * ρ2
    enh = TMS.ExactState(nh)
    @test hermiticity(enh) ≈ hermiticity(nh)
    @test hermiticity(eψ) == 1
    @test hermitianize(enh).state ≈ ket_bra_matrix(hermitianize(nh))
    @test hermitianize(eψ) === eψ
    # the limits of the MPS accepted, without effect
    @test hermitianize(enh; limits = Limits(maxdim = 1)).state ≈ hermitianize(enh).state
    @test mix(eψ; limits = Limits(maxdim = 1)).state ≈ mix(eψ).state
    @test dag(enh).state ≈ ket_bra_matrix(dag(nh))
    # unnormalised, measured through their trace
    for st in [2 * ψ, 3 * ρ]
        e = TMS.ExactState(st)
        @test trace(e) ≈ trace(st)
        @test expect(e, first(ops); complex = true) ≈ expect(st, first(ops); complex = true) atol =
              1e-12
        @test purity(e) ≈ purity(st)
        @test reduced_density_matrix(e, [1, 3]) ≈ reduced_density_matrix(st, [1, 3]) atol = 1e-12
    end
    @test normalize(eρ).state ≈ ket_bra_matrix(normalize(ρ))
    @test normalize(eψ).state ≈ raw_vector(normalize(ψ))
    # and back, normalised
    @test raw_vector(State(eψ)) ≈ eψ.state / norm(eψ.state)
    @test ket_bra_matrix(State(eρ)) ≈ eρ.state / LA.tr(eρ.state)
end

@testset "A state of definite charges" begin
    sys = System(4, Qubit(conserve = N))
    ψ = random_state(State{Pure}(sys, ["0", "1", "0", "1"]), 4)
    e = TMS.ExactState(ψ)
    # an operator of no definite flux measured by its part of no flux
    @test expect(e, X(1) * X(2); complex = true) ≈ expect(ψ, X(1) * X(2); complex = true)
    @test abs(expect(e, X(1); complex = true)) < 1e-12
    # a whole vector accepted with definite charges, refused otherwise
    v = zeros(16)
    v[[2, 3]] .= 1
    @test TMS.ExactState(sys, v).state ≈ v / norm(v)
    w = zeros(16)
    w[[1, 2]] .= 1
    @test_throws "no definite charge" TMS.ExactState(sys, w)
    # kept whole, its smallest component too
    u = zeros(16)
    u[2], u[3] = 1, 1e-14
    @test TMS.ExactState(sys, u).state[3] ≈ 1e-14 / norm(u)
    @test_throws "divides by the trace" TMS.ExactState(sys, zeros(16))
    # a projector of no definite charge refused, measured once weakened; a count kept charged
    @test_throws "collapse weaken(state, ()) instead" collapse(e, X(1) => 1.0)
    @test collapse(weaken(e, ()), X(1) => 1.0).state.state ≈
          TMS.ExactState(collapse(weaken(ψ, ()), X(1) => 1.0).state).state atol = 1e-10
    @test weaken(e).system == weaken(ψ).system
    total = N(1) + N(2) + N(3)
    @test last.(probabilities(e, total)) ≈ last.(probabilities(ψ, total)) atol = 1e-10
    @test collapse(e, total => 1).state.state ≈
          TMS.ExactState(collapse(ψ, total => 1).state).state atol = 1e-10
    # the sectors of the charge of the sites on the left of a cut
    for cut in 1:4
        s, ref = entanglement_by_sector(e, cut), entanglement_by_sector(ψ, cut)
        @test collect(keys(s)) == collect(keys(ref))
        for (c, x) in ref
            @test s[c].weight ≈ x.weight atol = 1e-12
            @test s[c].entropy ≈ x.entropy atol = 1e-10
            k = min(length(s[c].spectrum), length(x.spectrum))
            @test s[c].spectrum[1:k] ≈ x.spectrum[1:k] atol = 1e-12
        end
    end
    # a strong symmetry, whose mixed index full_matrix reads without weakening it
    ss = System(3, Fermion(conserve = strong(N)))
    ρs = mix(random_state(State{Pure}(ss, ["Occ", "Emp", "Occ"]), 3))
    es = TMS.ExactState(ρs)
    @test es.state ≈ ket_bra_matrix(ρs)
    hop = dag(C)(1) * C(2) + dag(C)(2) * C(1)
    @test expect(es, hop) ≈ expect(ρs, hop) atol = 1e-12
    @test reduced_density_matrix(es, [1, 3]) ≈ reduced_density_matrix(ρs, [1, 3]) atol = 1e-12
end

@testset "Entropies of a pure exact state from its Schmidt values" begin
    # of 13 sites out of 14, whose reduced density matrix, of dimension 8192, is not formed
    ψ = random_state(Pure, System(14, Qubit()), 4; rng = Xoshiro(20261108))
    e = TMS.ExactState(ψ)
    @test von_neumann_entropy(e, 1:13) ≈ entanglement_entropy(ψ, 13) atol = 1e-10
    @test von_neumann_entropy(e, 2:14) ≈ von_neumann_entropy(e, [1]) atol = 1e-10
    @test renyi2(e, 1:13) ≈ renyi2(ψ, 1:13) atol = 1e-10
    @test_throws "beyond max_dense = 16777216" reduced_density_matrix(e, 1:13)
    @test_throws "divides by the trace" von_neumann_entropy(TMS.ExactState(0 * ψ), [1])
end

@testset "Sectors of several quantities, of a modulus, strong" for sys in [
        System(3, Electron(conserve = (Ntot, 2Sz))),
        System(4, Fermion(conserve = parity(N))),
        System(3, Fermion(conserve = strong(N)))]
    states = sys[1] isa Electron ? ["UpDn", "Up", "Emp"] : fill("Occ", length(sys))
    states[end] = "Emp"
    ψ = random_state(Pure, sys, states, 4; rng = Xoshiro(20261106))
    e = TMS.ExactState(ψ)
    for cut in 1:length(sys)
        s, ref = entanglement_by_sector(e, cut), entanglement_by_sector(ψ, cut)
        @test collect(keys(s)) == collect(keys(ref))
        @test all(isapprox(s[c].weight, x.weight; atol = 1e-12) &&
                  isapprox(s[c].entropy, x.entropy; atol = 1e-12) for (c, x) in ref)
    end
end

@testset "The fidelity of two exact mixed states near rounding and approximations" begin
    # commuting states, against (Σ √(p q))² in BigFloat: eigenvalues of 1e-14 of the largest,
    # above the noise of the eigensolver, are kept
    D = 64
    sys = System(6, Qubit())
    δ = 64e-14
    p = [ (k == 1 ? 1 - δ : 0.0) + δ / D for k in 1:D ]
    ref = Float64(sum(sqrt.(big.(p) / D))^2)
    @test fidelity(TMS.ExactState(sys, LA.diagm(p)), TMS.ExactState(sys, LA.diagm(ones(D)))) ≈
          ref rtol = 1e-10
    # a defect of hermiticity: the Hermitian part, silently at rounding, with a warning beyond
    q3 = System(3, Qubit())
    rng = Xoshiro(20261107)
    # of full rank, so that the noise leaves its eigenvalues positive
    ρ = TMS.ExactState(q3, TMS.ExactState(random_state(Mixed, q3, 4; rng)).state + 0.1 * LA.I)
    g = randn(rng, ComplexF64, 8, 8)
    noisy = TMS.ExactState(q3, ρ.state + 1e-12 * g)
    @test @test_logs(fidelity(ρ, noisy)) ≈ fidelity(ρ, hermitianize(noisy))
    noisy = TMS.ExactState(q3, ρ.state + 1e-6 * g)
    @test @test_logs((:warn, r"Hermitian part"), fidelity(ρ, noisy)) ≈
          fidelity(ρ, hermitianize(noisy))
    # a negative eigenvalue within approximations is taken as zero, beyond it refused
    @test fidelity(TMS.ExactState(System(2, Qubit()), LA.diagm([1., 0.5, -1e-10, 1])),
                   TMS.ExactState(System(2, Qubit()), LA.diagm(ones(4)))) ≈
          (1 + sqrt(0.5) + 1)^2 / (2.5 * 4) rtol = 1e-9
end

# the paths a representation gets from its reduced_density_matrix and its purity, taken on a State
# and compared with its own
@testset "Generic paths of the measurements on $name" for (name, sys, h) in [
        ("qubits", System(4, Qubit()), X(1) * X(2) + Z(3)),
        ("fermions of no definite parity", System(4, Fermion()),
         dag(C)(1) * C(2) + dag(C)(2) * C(1) + N(4) + C(3)),
        ("fermions conserving N", System(4, Fermion(conserve = N)),
         dag(C)(1) * C(2) + dag(C)(2) * C(1) + N(4))]
    rng = Xoshiro(20261104)
    occ = ["Occ", "Emp", "Occ", "Emp"]
    charged = TMS.is_charged(sys)
    ψ = charged ? random_state(Pure, sys, occ, 4; rng) : random_state(Pure, sys, 4; rng)
    ρ = charged ? random_state(Mixed, sys, occ, 4; rng) : random_state(Mixed, sys, 4; rng)
    AS = TMS.AbstractState
    for kept in [[1], [4], [1, 2], [2, 3], [1, 3], [2, 4], [1, 2, 4]]
        for st in [ψ, mix(ψ), ρ]
            @test invoke(TMS.part_renyi2, Tuple{AS, Vector{Int}}, st, kept) ≈
                  TMS.part_renyi2(st, kept) atol = 1e-10
        end
        # the parts of a mixed state by the generic path, of the pure one by its own
        @test TMS.parts_renyi2(mix(ψ), kept) ≈ TMS.parts_renyi2(ψ, kept) atol = 1e-10
    end
    # the pure state through its MPO, against expect
    @test invoke(variance, Tuple{AS, Any}, ψ, h) ≈ variance(ψ, h) atol = 1e-10
end

@testset "From a whole vector or density matrix" begin
    q2 = System(2, Qubit())
    e = TMS.ExactState(q2, [1, 0, 0, 1])
    @test e isa TMS.ExactState{Pure}
    @test e.state ≈ [1, 0, 0, 1] / √2
    m = TMS.ExactState(q2, [1 0 0 0 ; 0 0 0 0 ; 0 0 0 0 ; 0 0 0 1])
    @test m isa TMS.ExactState{Mixed}
    @test LA.tr(m.state) ≈ 1
    @test_throws "a vector of 4 elements, not an array of size (3,)" TMS.ExactState(q2, [1, 0, 0])
    @test_throws "a density matrix of 4 × 4" TMS.ExactState(q2, ones(4, 3))
    # no bound by default, one given to protect oneself
    @test_throws "has 64 elements, beyond max_dense = 32" TMS.ExactState{Mixed}(System(3, Qubit()),
                                                                               "Up"; max_dense = 32)
    @test_throws "has 64 elements, beyond max_dense = 32" mix(
        TMS.ExactState{Pure}(System(3, Qubit()), "Up"); max_dense = 32)
    # the same call on a State, where it has no use
    @test mix(State{Pure}(System(3, Qubit()), "Up"); max_dense = 32) isa State{Mixed}
    @test size(reduced_density_matrix(State{Pure}(System(3, Qubit()), "Up"), 1:3;
                                      max_dense = nothing)) == (8, 8)
    # a whole array is taken whatever its size, the array of the wrong size refused
    @test_throws "a vector of 18446744073709551616 elements, not an array of size (0,)" TMS.ExactState(
        System(64, Qubit()), ComplexF64[])
    @test_throws "divides by the trace" TMS.ExactState(q2, zeros(4))
    traceless = Matrix(LA.Diagonal([1., 1, -1, -1]))
    @test_throws "zero up to rounding beside its norm" TMS.ExactState(q2, traceless)
    @test_throws "zero up to rounding beside its norm" TMS.ExactState(q2, traceless + 1e-17 * LA.I)
    @test LA.tr(TMS.ExactState(q2, 1e-20 * LA.I(4)).state) ≈ 1
end

@testset "Refusals of the measurements and conversion back" begin
    ψ = State{Pure}(System(2, Qubit()), "Up")
    @test_throws "expect divides by the trace" expect(TMS.ExactState(0 * ψ), Z(1))
    e = TMS.ExactState(ψ)
    @test_throws "is a generator" expect(e, Dissipator(Sm)(1))
    @test_throws "placed on sites" expect(e, X)
    m = mix(e)
    @test_throws "between a pure and a mixed" inner(e, m)
    @test_throws "between a pure and a mixed" inner(m, e)
    @test_logs (:warn, r"Hermitian part") fidelity(m, TMS.ExactState((1 + im) * mix(ψ)))
    negative = TMS.ExactState(System(2, Qubit()), LA.diagm([1., 0.5, -0.5, 1]))
    @test_throws "has the eigenvalue -0.25" fidelity(m, negative)
    @test_throws "negative probability" sample(negative)
    @test_throws "of probability 0" collapse(e, Z(1) => -1)
    @test_throws "meaningless on pure" dag(e)
    @test_throws "different systems" inner(e, TMS.ExactState{Pure}(System(3, Qubit()), "Up"))
    @test_throws "in increasing order" reduced_density_matrix(e, [2, 1])
    @test_throws "which the state does not have" reduced_density_matrix(e, [3])
    @test reduced_density_matrix(e, Int[]) == ones(1, 1)
    @test_throws "must be between 1 and 2" entanglement_entropy(e, 0)
    @test_throws "divides by the norm" entanglement_spectrum(TMS.ExactState(0 * ψ), 1)
    # on what the systems conserve, which an exact state does not depend on
    q = TMS.ExactState{Pure}(System(2, Qubit(conserve = N)), "Up")
    @test inner(q, e) ≈ 1
    big = TMS.ExactState(random_state(Pure, System(4, Qubit()), 4))
    @test maxlinkdim(State(big; limits = Limits(maxdim = 1))) == 1
end

@testset "Saved and loaded" begin
    in_tmp() do
        sys = System([Fermion(), Qubit()])
        for e in [TMS.ExactState(random_state(Mixed, sys, 2)), TMS.ExactState{Pure}(sys, "0")]
            save_state("s.h5", "e", e)
            l = load_state("s.h5", "e")
            @test l isa typeof(e)
            @test l.state == e.state
            @test l.system == e.system
            rm("s.h5")
        end
    end
end
