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
    @test only(entanglement_by_sector(eψ, 2))[2].entropy ≈ entanglement_entropy(ψ, 2) atol =
          1e-10
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
    # what rounding leaves in a sector the charges forbid, dropped as dense_state does
    u[1] = 1e-15
    eu = TMS.ExactState(sys, u)
    @test eu.state[1] == 0
    @test eu.state[3] ≈ 1e-14 / norm(u)
    @test length(entanglement_by_sector(eu, 2)) == length(entanglement_by_sector(State(eu), 2))
    ρu = TMS.ExactState(sys, u * u')
    # a mixture of sectors allowed, block diagonal: what joins two sectors dropped
    @test ρu.state[1, 2:end] == zeros(15) && ρu.state[2:end, 1] == zeros(15)
    @test real(ρu.state[1, 1]) > 0
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

@testset "Forbidden sectors dropped, of a strong quantity and of a modulus" begin
    occ = real(LA.diag(matrix(N, Fermion())))
    particles(b, n) = sum(occ[k + 1] for k in digits(b - 1; base = 2, pad = n))
    # strong: the blocks joining two charges and the block of another charge
    ss = System(3, Fermion(conserve = strong(N)))
    ρ0 = ket_bra_matrix(mix(random_state(State{Pure}(ss, ["Occ", "Emp", "Occ"]), 3)))
    two, one = findfirst(b -> particles(b, 3) == 2, 1:8), findfirst(b -> particles(b, 3) == 1, 1:8)
    ρ = copy(ρ0)
    ρ[two, one] = ρ[one, two] = ρ[one, one] = 1e-15
    es = TMS.ExactState(ss, ρ).state
    @test es[two, one] == es[one, two] == es[one, one] == 0
    @test es ≈ ρ0 atol = 1e-14
    # parity: 1 and 3 particles kept, 2 dropped
    ps = System(3, Fermion(conserve = parity(N)))
    i1, i2, i3 = (findfirst(b -> particles(b, 3) == k, 1:8) for k in 1:3)
    v = zeros(8)
    v[i1], v[i3], v[i2] = 1, 0.5, 1e-15
    w = TMS.ExactState(ps, v).state
    @test w[i3] ≈ 0.5 / norm([1, 0.5]) && w[i2] == 0
    v[i2] = 1
    @test_throws "no definite charge" TMS.ExactState(ps, v)
end

@testset "Entropies of a pure exact state from its Schmidt values" begin
    # of 13 sites out of 14, whose reduced density matrix, of dimension 8192, is not formed
    ψ = random_state(Pure, System(14, Qubit()), 4; rng = Xoshiro(20261108))
    e = TMS.ExactState(ψ)
    @test von_neumann_entropy(e, 1:13) ≈ entanglement_entropy(ψ, 13) atol = 1e-10
    # positions of another type converted, the reduced density matrix still not formed
    @test von_neumann_entropy(e, collect(1.0:13.0)) ≈ entanglement_entropy(ψ, 13) atol = 1e-10
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
    # the others renormalised, which keeps the fidelity of a state with itself at one
    near = TMS.ExactState(System(2, Qubit()), LA.diagm([1., 0.5, -1e-9, 0.5]))
    @test fidelity(near, near) ≈ 1 atol = 1e-14
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
    # a bound given is kept, which by default is the size of the state
    em = TMS.ExactState{Mixed}(System(3, Qubit()), "Up")
    @test_throws "beyond max_dense = 32" reduced_density_matrix(em, 1:3; max_dense = 32)
    @test size(reduced_density_matrix(em, 1:3; max_dense = 64)) == (8, 8)
    @test_throws "beyond max_dense = 32" measure(em, VonNeumannEntropy(1:3; max_dense = 32))
    @test TMS.dense_bound(em) == 2^24
    @test_throws "has 64 elements, beyond max_dense = 32" mix(
        TMS.ExactState{Pure}(System(3, Qubit()), "Up"); max_dense = 32)
    # the same call on a State, where it has no use
    @test mix(State{Pure}(System(3, Qubit()), "Up"); max_dense = 32) isa State{Mixed}
    @test size(reduced_density_matrix(State{Pure}(System(3, Qubit()), "Up"), 1:3;
                                      max_dense = nothing)) == (8, 8)
    # a whole array is taken whatever its size, the array of the wrong size refused
    msg = "a vector of 18446744073709551616 elements, not an array of size (0,)"
    @test_throws msg TMS.ExactState(System(64, Qubit()), ComplexF64[])
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
    # of a matrix that is not Hermitian, its Hermitian part
    bad = TMS.ExactState((1 + im) * mix(ψ))
    @test (@test_logs (:warn, r"Hermitian part") fidelity(m, bad)) ≈
          fidelity(m, hermitianize(bad))
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

@testset "Operators, channels and generators applied on $name" for (name, sys, ops) in [
        ("qubits", System(3, Qubit()), [X(1) * Z(3) + 2 * Y(2), Sp(1) * Sm(3)]),
        ("fermions", System(3, Fermion()), [dag(C)(1) * C(3) + dag(C)(3) * C(1), C(2)]),
        ("a qubit among fermions", System([Fermion(), Qubit(), Electron()]),
         [dag(C)(1) * Cup(3) + dag(Cup)(3) * C(1), X(2) * Ntot(3)])]
    rng = Xoshiro(20261109)
    ψ, ρ = random_state(Pure, sys, 4; rng), random_state(Mixed, sys, 4; rng)
    eψ, eρ = TMS.ExactState(ψ), TMS.ExactState(ρ)
    # as apply applies them to a State, not normalised: an operator, K ρ K† on a density matrix
    for o in ops
        @test apply(o, eψ).state ≈ raw_vector(apply(o, ψ)) atol = 1e-12
        @test apply(o, eρ).state ≈ ket_bra_matrix(apply(o, ρ)) atol = 1e-12
    end
    # a channel and a generator on a density matrix, a channel of one Kraus operator on a vector
    for c in [Gate(ops[1]), dephasing_channel(0.3)(1), Dissipator(ops[1])]
        @test apply(c, eρ).state ≈ ket_bra_matrix(apply(c, ρ)) atol = 1e-12
    end
    @test apply(Gate(ops[1]), eψ).state ≈ raw_vector(apply(ops[1], ψ)) atol = 1e-12
    @test_throws "more than one Kraus operator" apply(dephasing_channel(0.3)(1), eψ)
    @test_throws "is a generator" apply(Dissipator(ops[1]), eψ)
    @test_throws "rather than its MPO" apply(make_mpo(ψ, ops[1]), eψ)
    # the limits of the MPS without use
    @test apply(ops[1], eψ; limits = Limits(maxdim = 1)).state ≈ apply(ops[1], eψ).state
end

@testset "Evolution of an exact state" begin
    # a Hamiltonian, against the exponential of its matrix; a hop past a site on fermions
    # conserving N
    for (sys, h, st) in [(System(3, Qubit()), X(1) * X(2) + Z(2) * Z(3) + 0.5 * Z(1), "+"),
                         (System(3, Fermion(conserve = N)),
                          dag(C)(1) * C(3) + dag(C)(3) * C(1) + N(2), ["Occ", "Emp", "Emp"])]
        ψ = TMS.ExactState{Pure}(sys, st)
        u = exp(-0.7im * matrix(h, sys))
        @test evolve(-im * h, 0.7, ψ; nsteps = 2).state ≈ u * ψ.state atol = 1e-11
        # the method and the limits of the MPS without use
        @test evolve(-im * h, 0.7, ψ; nsteps = 2, algo = ApproxW(order = 2),
                     limits = Limits(maxdim = 1)).state ≈ u * ψ.state atol = 1e-11
    end
    # a Lindbladian, against the exponential of its superoperator built here with kron:
    # vec(A ρ B) = (transpose(B) ⊗ A) vec(ρ)
    sys = System(2, Qubit())
    x, z, sm, o = (matrix(a, Qubit()) for a in (X, Z, Sm, Id))
    hm, lm = kron(x, x) + 0.3 * kron(z, o), sqrt(0.4) * kron(sm, o)
    id = kron(o, o)
    side(a, b) = kron(transpose(b), a)
    lind = -im * (side(hm, id) - side(id, hm)) + side(lm, lm') -
           0.5 * (side(lm' * lm, id) + side(id, lm' * lm))
    ρ = TMS.ExactState{Mixed}(sys, ["Up", "+"])
    exact = reshape(exp(1.1 * lind) * vec(ρ.state), 4, 4)
    evolver = -im * (X(1) * X(2) + 0.3 * Z(1)) + Dissipator(sqrt(0.4) * Sm)(1)
    @test evolve(evolver, 1.1, ρ; nsteps = 3).state ≈ exact atol = 1e-11
    # depending on time, against fourth-order Runge-Kutta in steps of 1e-3, of error 1e-12
    hd(t) = cos(2t) * kron(x, o) + kron(z, z)
    ψ = TMS.ExactState{Pure}(sys, "Up")
    v = ψ.state
    δ = 1e-3
    for k in 0:999
        t = k * δ
        f(t, w) = -im * hd(t) * w
        k1 = f(t, v); k2 = f(t + δ / 2, v + δ / 2 * k1)
        k3 = f(t + δ / 2, v + δ / 2 * k2); k4 = f(t + δ, v + δ * k3)
        v += δ / 6 * (k1 + 2k2 + 2k3 + k4)
    end
    terms = [-im * X(1) => t -> cos(2t), -im * Z(1) * Z(2)]
    # the error bounded over the whole evolution, whatever its number of steps
    for nsteps in (1, 2, 50)
        @test evolve(terms, 1.0, ψ; nsteps).state ≈ v atol = 1e-10
    end
    # the evolver itself, not its PreMPO
    @test_throws "rather than its PreMPO" evolve(PreMPO(State{Pure}(sys, "Up"), -im * X(1)), 1.0,
                                                 ψ; nsteps = 1)
end

# a phase of one's own creating an exact state, for the phases of the library to act on
Base.@kwdef struct CreateExact <: TMS.AbstractPhase
    name::String = "exact"
    time_start = nothing
    final_measurements = []
    state::TMS.ExactState
end

TMS.run_phase(p::CreateExact, sim::TMS.Simulation, ::Int) = TMS.Simulation(sim, p.state)

@testset "Searches and thermal states of an exact state" begin
    sys = System(4, Qubit())
    h = X(1) * X(2) + X(2) * X(3) + X(3) * X(4) + 0.7 * sum(Z(i) for i in 1:4)
    λ = LA.eigvals(LA.Hermitian(matrix(h, sys)))
    ψ = TMS.ExactState{Pure}(sys, "+")
    # the lowest eigenvalue, then the next one orthogonal to it, the method and the limits of
    # the MPS without use and not needed
    e, g = ground_state(h, ψ)
    @test e ≈ λ[1] atol = 1e-11
    @test real(expect(g, h)) ≈ λ[1] atol = 1e-11
    e1, g1 = ground_state(h, ψ; orthogonal_to = [g], algo = Dmrg(nsweeps = 1),
                          limits = Limits(maxdim = 1))
    @test e1 ≈ λ[2] atol = 1e-10
    @test abs(inner(g, g1)) < 1e-10
    # a degenerate lowest level, reached state after state from the same start, whose
    # perturbation is drawn anew for each number of states to avoid
    p = (X(1) + Z(2)) * (X(2) + Z(3)) * (X(3) + Z(4))
    hd = p + dag(p)
    λd = minimum(LA.eigvals(LA.Hermitian(matrix(hd, sys))))
    found = TMS.ExactState{Pure}[]
    for _ in 1:3
        ed, gd = ground_state(hd, ψ; orthogonal_to = found)
        @test ed ≈ λd atol = 1e-10
        @test all(r -> abs(inner(r, gd)) < 1e-10, found)
        push!(found, gd)
    end
    # a weight below the gap, as for a State
    @test_logs (:warn, r"does not exceed the gap") ground_state(h, ψ; orthogonal_to = [g],
                                                                 weight = 1e-3)
    # in the sector of the state given, on fermions conserving N, against dmrg
    f = System(4, Fermion(conserve = N))
    hop = sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3) + 0.5 * N(1)
    start = State{Pure}(f, ["Occ", "Emp", "Occ", "Emp"])
    ef, _ = ground_state(hop, TMS.ExactState(start))
    em, _ = ground_state(hop, random_state(start, 4; rng = Xoshiro(1));
                         algo = Dmrg(nsweeps = 10))
    @test ef ≈ em atol = 1e-8
    # the steady state of two driven qubits decaying, whose Z is known
    Ω, γ = 0.8, 0.5
    L = sum(-im * Ω / 2 * X(i) + Dissipator(sqrt(γ) * Sm)(i) for i in 1:2)
    v, ρ = steady_state(L, TMS.ExactState{Mixed}(System(2, Qubit()), "FullyMixed"))
    @test abs(v) < 1e-20
    @test LA.tr(ρ.state) ≈ 1
    @test real(expect1(ρ, Z)) ≈ fill(-γ^2 / (γ^2 + 2Ω^2), 2) atol = 1e-11
    # e^{-βH}/Z from the fully mixed state, and log Z - Σ log d
    m = matrix(h, sys)
    w = exp(-1.3 * m)
    l, th = thermalize(h, 1.3, TMS.ExactState{Mixed}(sys, "FullyMixed"); nsteps = 3)
    @test l ≈ log(LA.tr(w) / 16)
    @test th.state ≈ w / LA.tr(w) atol = 1e-11
    # refused on the other representation, and given a PreMPO
    @test_throws "searches a pure state" ground_state(h, mix(ψ))
    @test_throws "searches a density matrix" steady_state(L, ψ)
    @test_throws "needs a mixed representation" thermalize(h, 1.0, ψ; nsteps = 1)
    @test_throws "rather than its PreMPO" ground_state(PreMPO(State{Pure}(sys, "+"), h), ψ)
    # the phases call the searches of the exact state once, measured at their end
    in_tmp() do
        phases = [CreateExact(state = ψ), GroundState(hamiltonian = h, algo = Dmrg(nsweeps = 2),
                                                      measurements = "e" => [:energy])]
        sim = runTMS(SimData(; name = "g", phases))
        @test sim.state isa TMS.ExactState{Pure}
        @test sim.energy ≈ λ[1] atol = 1e-11
        phases = [CreateExact(state = TMS.ExactState{Mixed}(System(2, Qubit()), "FullyMixed")),
                  SteadyState(lindbladian = L, algo = Dmrg(nsweeps = 2))]
        sim = runTMS(SimData(; name = "s", phases))
        @test real(expect(sim.state, Z(1))) ≈ -γ^2 / (γ^2 + 2Ω^2) atol = 1e-11
    end
end

@testset "Searches and evolutions of an exact state hard for Krylov" begin
    x, y, z, sm = (matrix(a, Qubit()) for a in (X, Y, Z, Sm))
    # steady states of a driven chain decaying slowly, against the kernel of its superoperator
    # built here with kron, solved with the trace as a constraint
    for (n, γ) in [(4, 0.1), (4, 0.01), (5, 0.1)]
        site(m, i) = on_qubit(m, i, n)
        id = site([1. 0.; 0. 1.], 1)
        side(a, b) = kron(transpose(b), a)
        hm = sum(0.5 * site(x, i) for i in 1:n) + sum(site(z, i) * site(z, i + 1) for i in 1:n-1)
        lind = -im * (side(hm, id) - side(id, hm))
        for i in 1:n
            l = sqrt(γ) * site(sm, i)
            lind += side(l, l') - 0.5 * (side(l' * l, id) + side(id, l' * l))
        end
        lind[1, :] = vec(id)
        ref = reshape(lind \ [1; zeros(4^n - 1)], 2^n, 2^n)
        L = sum(-im * 0.5 * X(i) + Dissipator(sqrt(γ) * Sm)(i) for i in 1:n) -
            im * sum(Z(i) * Z(i + 1) for i in 1:n-1)
        _, ρ = steady_state(L, TMS.ExactState{Mixed}(System(n, Qubit()), "FullyMixed"))
        @test ρ.state ≈ ref atol = 1e-9
    end
    # e^{-βH}/Z at large β in one step, and log Z - Σ log d
    sys = System(4, Qubit())
    h = X(1) * X(2) + X(2) * X(3) + X(3) * X(4) + 0.7 * sum(Z(i) for i in 1:4)
    λ, V = LA.eigen(LA.Hermitian(matrix(h, sys)))
    β = 200
    w = exp.(-β * (λ .- λ[1]))
    l, th = thermalize(h, β, TMS.ExactState{Mixed}(sys, "FullyMixed"); nsteps = 1)
    @test l ≈ -β * λ[1] + log(sum(w) / 16) atol = 1e-8
    @test th.state ≈ V * LA.Diagonal(w / sum(w)) * V' atol = 1e-10
    # a long evolution in one step
    ψ = TMS.ExactState{Pure}(sys, "+")
    @test evolve(-im * h, 1000, ψ; nsteps = 1).state ≈
          V * (exp.(-1000im * λ) .* (V' * ψ.state)) atol = 1e-10
    # the norm decreasing along the exponential, by a constant shift of the energies
    shifted = h + 30 * Id(1)
    l20, th20 = thermalize(h, 20, TMS.ExactState{Mixed}(sys, "FullyMixed"); nsteps = 1)
    ls, ths = thermalize(shifted, 20, TMS.ExactState{Mixed}(sys, "FullyMixed"); nsteps = 1)
    @test ls ≈ l20 - 600 atol = 1e-10
    @test ths.state ≈ th20.state atol = 1e-10
    v = V * (exp.(-20 * (λ .+ 30)) .* (V' * ψ.state))
    @test norm(evolve(-shifted, 20, ψ; nsteps = 1).state - v) / norm(v) < 1e-10
    # a norm beyond the floats refused, an evolution of no duration leaving the state as it is
    @test_throws "beyond the floats" evolve(-h, 300, ψ; nsteps = 1)
    @test evolve([-im * h, -im * X(1) => t -> cos(t)], 0.0, ψ; nsteps = 1).state == ψ.state
    # a field rotating fast: H(t) = ω Z / 2 + Ω (cos ωt X + sin ωt Y) / 2 is solved by
    # e^{-iωtZ/2} e^{-iΩtX/2}
    ω, Ω, t = 200, 1, 1.0
    terms = [-im * ω / 2 * Z(1), -im * Ω / 2 * X(1) => s -> cos(ω * s),
             -im * Ω / 2 * Y(1) => s -> sin(ω * s)]
    φ = TMS.ExactState{Pure}(System(1, Qubit()), "Up")
    @test evolve(terms, t, φ; nsteps = 1).state ≈
          exp(-im * ω * t / 2 * z) * exp(-im * Ω * t / 2 * x) * φ.state atol = 1e-10
    # orthogonal to pure states only, given once each
    e, g = ground_state(h, ψ)
    @test_throws "orthogonal to pure states" ground_state(h, ψ; orthogonal_to = [mix(g)])
    @test_throws "orthogonal to pure states" ground_state(h, ψ; orthogonal_to = [mix(State(g))])
    # a state the others span, exactly or nearly, adding nothing to avoid
    e1, g1 = ground_state(h, ψ; orthogonal_to = [g])
    @test ground_state(h, ψ; orthogonal_to = [g, g]).energy ≈ e1 atol = 1e-10
    d = normalize(randn(Xoshiro(3), ComplexF64, 16))
    for ε in (1e-12, 1e-9)
        near = TMS.ExactState(sys, g.state + ε * d)
        @test ground_state(h, ψ; orthogonal_to = [g, near]).energy ≈ e1 atol = 1e-10
    end
    both = TMS.ExactState(sys, (g.state + g1.state) / √2)
    @test ground_state(h, ψ; orthogonal_to = [g, g1, both]).energy ≈
          ground_state(h, ψ; orthogonal_to = [g, g1]).energy atol = 1e-10
    # the Thermalize and Gates phases
    in_tmp() do
        phases = [CreateExact(state = TMS.ExactState{Mixed}(sys, "FullyMixed")),
                  Thermalize(hamiltonian = h, beta = 1.3, nsteps = 3)]
        sim = runTMS(SimData(; name = "t", phases))
        m = exp(-1.3 * matrix(h, sys))
        @test sim.state.state ≈ m / LA.tr(m) atol = 1e-11
        @test sim.beta == 1.3
        @test sim.log_trace ≈ log(LA.tr(m) / 16) atol = 1e-12
        # the error bounded over the whole phase, whatever its number of steps
        phases = [CreateExact(state = TMS.ExactState{Mixed}(sys, "FullyMixed")),
                  Thermalize(hamiltonian = h, beta = 20, nsteps = 1000)]
        sim = runTMS(SimData(; name = "t", phases))
        @test sim.log_trace ≈ -20 * λ[1] + log(sum(exp.(-20 * (λ .- λ[1]))) / 16) atol = 1e-10
        # a pure state refused, even when its steps round to none
        phases = [CreateExact(state = ψ), Thermalize(hamiltonian = h, beta = 0.01,
                                                     beta_step = 0.1)]
        @test_throws "needs a mixed representation" runTMS(SimData(; name = "p", phases))
        two = TMS.ExactState{Pure}(System(2, Qubit()), "Up")
        sim = runTMS(SimData(; name = "g", phases = [CreateExact(state = two),
                                                      Gates(gates = X(1) * X(2))]))
        @test sim.state.state ≈ kron(x, x) * two.state
    end
end

@testset "Gates on an exact state beside a State, stopped between the steps of the State" begin
    # the steps of the State kept, the exact state given the gates at once at the first one; a
    # run stopped after each step ends as the uninterrupted one
    in_tmp() do
        circuit = CX(1, 5) * Rzz(0.4)(2, 3) * CX(4, 6) * prod(H(i) for i in 1:6) * CZ(2, 6)
        s6 = System(6, Qubit())
        phases = [CreateExact(state = TMS.ExactState{Pure}(s6, "Up")),
                  CreateState(type = Pure, system = s6, state = "Up", on = :b),
                  Gates(gates = circuit, on = [:state, :b])]
        sd(name) = SimData(; name, phases)
        ref = runTMS(sd("ref"))
        sim = run_through(sd("chk"); max_time = -1)
        @test ref.state.state ≈ Vector(ref.b) atol = 1e-12
        @test sim.state.state == ref.state.state
        @test Vector(sim.b) ≈ Vector(ref.b) atol = 1e-12
        nsteps = length(last(TMS.apply_steps(circuit, ref.b)))
        @test nsteps > 3
        @test count("Resuming from the checkpoint", read("chk/log", String)) > nsteps
    end
end

@testset "Searches of an exact state leave the sectors of the symmetries not declared" begin
    # an Ising chain conserves the parity of its X, undeclared: from a state of one parity, the
    # levels of the other reached, against the matrix built here with kron
    n = 4
    sys = System(n, Qubit())
    x, z = matrix(X, Qubit()), matrix(Z, Qubit())
    hd = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - 0.3 * sum(X(i) for i in 1:n)
    hm = -sum(on_qubit(z, i, n) * on_qubit(z, i + 1, n) for i in 1:n-1) -
         0.3 * sum(on_qubit(x, i, n) for i in 1:n)
    λ = LA.eigvals(LA.Hermitian(hm))
    ψ = TMS.ExactState{Pure}(sys, "+")
    e, g = ground_state(hd, ψ)
    @test e ≈ λ[1] atol = 1e-12
    @test ground_state(hd, ψ; orthogonal_to = [g]).energy ≈ λ[2] atol = 1e-12
    odd = TMS.ExactState{Pure}(sys, ["+", "+", "+", "-"])
    @test ground_state(hd, odd).energy ≈ λ[1] atol = 1e-12
    # among degenerate ground states, the one the state given picks
    _, up = ground_state(-sum(Z(i) * Z(i + 1) for i in 1:n-1), TMS.ExactState{Pure}(sys, "Up"))
    @test abs2(up.state[1]) > 0.9999
end

@testset "Backends of the states of a run" begin
    in_tmp() do
        sys = System(3, Qubit())
        ψ = State{Pure}(sys, "+")
        e = TMS.ExactState(ψ)
        # created, from a description or a state, the type given or that of the run
        sim = runTMS(SimData(name = "c", phases = [
            CreateState{Pure}(sys, "+"; backend = TMS.ExactState)]))
        @test sim.state isa TMS.ExactState{Pure} && sim.state.state ≈ e.state
        sim = runTMS(SimData(name = "s", phases = [CreateState(ψ; backend = TMS.ExactState)]))
        @test sim.state.state ≈ e.state
        sim = runTMS(SimData(name = "e", phases = [CreateState(e)]))
        @test sim.state isa State{Pure} && Vector(sim.state) ≈ e.state
        mixture = Mixture([0.25 => "Up", 0.75 => "Dn"])
        sim = runTMS(SimData(name = "m", phases = [
            CreateState{Mixed}(sys, mixture; backend = TMS.ExactState)]))
        @test sim.state.state ≈ TMS.ExactState{Mixed}(sys, mixture).state
        # a whole program on the backend of the run, against the same on MPS converted, the
        # mixing keeping the backend of the state
        circuit = CX(1, 2) * H(3)
        phases = [CreateState{Pure}(sys, "+"), Gates(gates = circuit), ToMixed(),
                  Gates(gates = X(2))]
        mps = runTMS(SimData(; name = "mps", phases))
        exact = runTMS(SimData(; name = "exact", backend = TMS.ExactState, phases))
        @test exact.state isa TMS.ExactState{Mixed}
        @test exact.state.state ≈ TMS.ExactState(mps.state).state atol = 1e-12
        @test occursin("States created as ExactState by default", read("exact/log", String))
        # converted, backend and representation at once, or back to an MPS truncated once
        sim = runTMS(SimData(name = "cm", phases = [
            CreateState(ψ), ConvertState{Mixed}(backend = TMS.ExactState)]))
        @test sim.state isa TMS.ExactState{Mixed} && sim.state.state ≈ e.state * e.state'
        # cos(0.3) |000⟩ + sin(0.3) |110⟩, cut to its largest Schmidt value, not normalized again
        v = zeros(8)
        v[1], v[7] = cos(0.3), sin(0.3)
        sim = runTMS(SimData(name = "cp", phases = [
            CreateState(TMS.ExactState(sys, v); backend = TMS.ExactState),
            ConvertState{Pure}(backend = State, limits = Limits(maxdim = 1))]))
        @test sim.state isa State{Pure} && maxlinkdim(sim.state) == 1
        @test abs(Vector(sim.state)[1]) ≈ cos(0.3)
        # charged fermions, against their vector contracted here
        f = System(3, Fermion(conserve = N))
        φ = random_state(State{Pure}(f, ["Occ", "Emp", "Occ"]), 3; rng = Xoshiro(7))
        sim = runTMS(SimData(name = "f", phases = [
            CreateState(φ), ConvertState{Mixed}(backend = TMS.ExactState)]))
        @test sim.state.state ≈ dense_vector(φ) * dense_vector(φ)' atol = 1e-12
        # a state given randomized once a State, an exact state saved and loaded with limits
        sim = runTMS(SimData(name = "rs", seed = 1, phases = [CreateState(e; randomize = 2)]))
        @test sim.state isa State{Pure} && norm(sim.state) ≈ 1
        runTMS(SimData(name = "we", phases = [CreateState(e; backend = TMS.ExactState),
                                              SaveState(file = "e.h5")]))
        sim = runTMS(SimData(name = "le", phases = [
            LoadState(file = "../we/e.h5", limits = Limits(maxdim = 1))]))
        @test sim.state isa TMS.ExactState{Pure} && sim.state.state ≈ e.state
        # an exact run stopped after each step resumes as the uninterrupted one
        sd(name) = SimData(; name, backend = TMS.ExactState, phases)
        ref = runTMS(sd("ref"))
        sim = run_through(sd("chk"); max_time = -1)
        @test sim.state.state == ref.state.state
        @test count("Resuming from the checkpoint", read("chk/log", String)) > 0
        # the positional constructors of 2.0
        @test SimData("p", "", 0., [], "%.13g", nothing, [CreateState(ψ)]).backend === State
        @test CreateState("c", nothing, [], Pure, sys, "Up", 0, :state).backend === nothing
        @test CreateState{Pure}("c", nothing, [], Pure, sys, "Up", 0, :state).on === :state
        @test LoadState("l", nothing, [], "s.h5", "state", Limits(), :state).backend === nothing
        # loaded as written, or converted
        runTMS(SimData(name = "w", phases = [CreateState(ψ), SaveState(file = "s.h5")]))
        sim = runTMS(SimData(name = "l", phases = [LoadState(file = "../w/s.h5")]))
        @test sim.state isa State{Pure}
        sim = runTMS(SimData(name = "lx", phases = [
            LoadState(file = "../w/s.h5", backend = TMS.ExactState)]))
        @test sim.state isa TMS.ExactState{Pure} && sim.state.state ≈ e.state
        # refused: a backend with its representation or not a type, a mixed state made pure
        @test_throws "without its representation" CreateState{Pure}(
            sys, "Up"; backend = TMS.ExactState{Pure})
        @test_throws "a type of states" SimData(backend = 3, phases = [CreateState(ψ)])
        for b in (TMS.AbstractState, Union{}, Union{State, TMS.ExactState})
            @test_throws "a type of states" CreateState(ψ; backend = b)
        end
        @test_throws "cannot be made pure" runTMS(SimData(name = "p", phases = [
            CreateState{Mixed}(sys, "Up"; backend = TMS.ExactState),
            ConvertState{Pure}(backend = State)]))
    end
end

@testset "Options ignored in a run and partial traces of an exact state" begin
    in_tmp() do
        sys = System(3, Qubit())
        h = X(1) * X(2) + Z(3)
        evolve_with(name; kwargs...) = runTMS(SimData(; name, backend = TMS.ExactState, phases = [
            CreateState{Pure}(sys, "Up"),
            Evolve(; evolver = -im * h, duration = 0.5, nsteps = 2, kwargs...)]))
        evolve_with("given"; algo = ApproxW(order = 2), limits = Limits(maxdim = 4))
        @test occursin("ExactState ignores algo and limits", read("given/log", String))
        evolve_with("default")
        @test !occursin("ignores", read("default/log", String))
        # the link dimension of randomize, and the limits of Gates once for its two states
        sys2 = System(2, Qubit())
        runTMS(SimData(name = "gl", seed = 1, backend = TMS.ExactState, phases = [
            CreateState{Pure}(sys2, "Up"; randomize = 4),
            CreateState{Pure}(sys2, "Up"; on = :b),
            Gates(gates = X(1), limits = Limits(maxdim = 2), on = [:state, :b]),
            SaveState(file = "s.h5"),
            LoadState(file = "s.h5", backend = TMS.ExactState, limits = Limits(maxdim = 2))]))
        log = read("gl/log", String)
        @test occursin("ExactState ignores the link dimension of randomize", log)
        @test count("ExactState ignores limits", log) == 2
    end
    # the density matrix of the sites kept, against partial_trace of the State, on qubits and on
    # charged fermions, by keep and by trace_out
    q = mix(random_state(Pure, System(4, Qubit()), 4; rng = Xoshiro(3)))
    f = mix(random_state(State{Pure}(System(4, Fermion(conserve = N)),
                                     ["Occ", "Emp", "Occ", "Emp"]), 4; rng = Xoshiro(4)))
    for ρ in (q, f), args in ((; keep = [1, 3]), (; trace_out = [2]), (; keep = [4]))
        e = partial_trace(TMS.ExactState(ρ); args...)
        @test e.system == partial_trace(ρ; args...).system
        @test e.state ≈ TMS.ExactState(partial_trace(ρ; args...)).state atol = 1e-12
    end
    # of the trace of the state, as for a State, a traceless one included
    ρa, ρb = mix(State{Pure}(System(3, Qubit()), "Up")), mix(State{Pure}(System(3, Qubit()), "+"))
    for (ex, st) in ((2 * TMS.ExactState(ρa), 2 * ρa),
                     (TMS.ExactState(ρa) - TMS.ExactState(ρb), ρa - ρb))
        @test partial_trace(ex; keep = [1]).state ≈
              TMS.ExactState(partial_trace(st; keep = [1])).state atol = 1e-12
    end
    # converted to an MPS of its norm or its trace
    @test norm(State(2 * TMS.ExactState(State{Pure}(System(3, Qubit()), "+")))) ≈ 2
    @test trace(State(2 * TMS.ExactState(ρa))) ≈ 2
    δ = TMS.ExactState(ρa) - TMS.ExactState(ρb)
    @test TMS.ExactState(State(δ)).state ≈ δ.state atol = 1e-12
    @test_throws "needs a mixed representation" partial_trace(TMS.ExactState(State{Pure}(
        System(2, Qubit()), "Up")); keep = [1])
    @test_throws "every site" partial_trace(TMS.ExactState(q); trace_out = 1:4)
    # the phase in a run
    in_tmp() do
        sim = runTMS(SimData(name = "p", backend = TMS.ExactState, phases = [
            CreateState{Mixed}(System(3, Qubit()), ["Up", "+", "Dn"]),
            PartialTrace(keep = [2, 3])]))
        @test sim.state.state ≈ TMS.ExactState{Mixed}(System(2, Qubit()), ["+", "Dn"]).state
    end
end

@testset "Random exact states and their link dimension" begin
    # uniform on the sphere: the mean weight of a basis vector is 1/d
    q = System(3, Qubit())
    rng = Xoshiro(11)
    w = [ abs2(random_state(TMS.ExactState{Pure}, q, 4; rng).state[1]) for _ in 1:2000 ]
    @test sum(w) / length(w) ≈ 1 / 8 atol = 0.01
    # its second moment, 2 / (d (d + 1)) on the sphere, about 4.7 standard deviations wide
    @test sum(w .^ 2) / length(w) ≈ 2 / 72 atol = 0.005
    # a mixed state of the Hilbert-Schmidt measure, of mean purity 2d / (d² + 1)
    p = [ real(LA.tr(random_state(TMS.ExactState{Mixed}, System(2, Qubit()), 4; rng).state^2))
          for _ in 1:2000 ]
    @test sum(p) / length(p) ≈ 8 / 17 atol = 0.01
    # in the sector of the local states, pure and mixed, the mixed one block diagonal over the
    # sectors its purification reaches, against the number of particles counted here
    f = System(3, Fermion(conserve = N))
    occ = real(LA.diag(matrix(N, Fermion())))
    count_n = [ sum(occ[k + 1] for k in digits(b - 1; base = 2, pad = 3)) for b in 1:8 ]
    ψ = random_state(TMS.ExactState{Pure}, f, ["Occ", "Emp", "Occ"], 4; rng)
    @test norm(ψ.state) ≈ 1
    @test all(iszero(ψ.state[b]) for b in 1:8 if count_n[b] ≠ 2)
    @test random_state(ψ, 4; rng).state[count_n .≠ 2] == zeros(count(≠(2), count_n))
    ρ = random_state(TMS.ExactState{Mixed}, f, ["Occ", "Emp", "Occ"], 4; rng).state
    nm = LA.Diagonal(count_n)
    @test LA.tr(ρ) ≈ 1 && ρ ≈ ρ' && norm(ρ * nm - nm * ρ) < 1e-12
    @test minimum(LA.eigvals(LA.Hermitian(ρ))) > -1e-12
    σ = random_state(TMS.ExactState{Mixed}, q, 4; rng).state
    @test LA.tr(σ) ≈ 1 && minimum(LA.eigvals(LA.Hermitian(σ))) > -1e-12
    # refused as for a State
    @test_throws "conserves something" random_state(TMS.ExactState{Pure}, f, 4)
    @test_throws "strongly" random_state(TMS.ExactState{Mixed}, System(2, Fermion(conserve =
        strong(N))), ["Occ", "Emp"], 4)
    @test_throws "local states" random_state(TMS.ExactState{Mixed}, q, FullyMixed(), 4)
    # drawn in a run
    in_tmp() do
        sim = runTMS(SimData(name = "r", seed = 2, backend = TMS.ExactState, phases = [
            CreateState{Pure}(f, ["Occ", "Emp", "Occ"]; randomize = 4)]))
        @test sim.state isa TMS.ExactState{Pure} && all(iszero(sim.state.state[b]) for b in 1:8
                                                        if count_n[b] ≠ 2)
    end
    # the link dimension of the whole space
    @test maxlinkdim(TMS.ExactState{Pure}(System(4, Qubit()), "Up")) == 4
    @test maxlinkdim(TMS.ExactState{Mixed}(q, "Up")) == 4
    @test maxlinkdim(TMS.ExactState{Pure}(System([Qubit(), Qudit(3), Qubit()]), "0")) == 2
    @test maxlinkdim(TMS.ExactState{Pure}(System(1, Qubit()), "Up")) == 1
end

@testset "Arithmetic, truncation and vector of an exact state" begin
    sys = System(3, Qubit())
    a, b = State{Pure}(sys, "Up"), State{Pure}(sys, "+")
    ea, eb = TMS.ExactState(a), TMS.ExactState(b)
    # as for a State, without truncation
    @test Vector(2ea - eb / 3 + ea * 1im) ≈ Vector(2a - b / 3 + a * 1im)
    @test Vector(-ea) == -Vector(ea)
    ρ = mix(ea) + mix(eb)
    @test ρ.state ≈ Vector(ea) * Vector(ea)' + Vector(eb) * Vector(eb)'
    @test truncate(ea; limits = Limits(maxdim = 1)) === ea
    # a copy of the vector
    v = Vector(ea)
    v[1] = 0
    @test Vector(ea)[1] == 1
    @test_throws "beyond max_dense = 4" Vector(ea; max_dense = 4)
    @test_throws "different systems" ea + TMS.ExactState{Pure}(System(3, Qubit(conserve = N)),
                                                               "Up")
end

@testset "Products, matrix elements and comparisons with a State of exact states" begin
    rng = Xoshiro(5)
    for (sa, sb) in ((System(2, Qubit()), System(1, Qubit())),
                     (System(2, Fermion()), System(2, Fermion())))
        a, b, c = (random_state(Pure, s, 2; rng) for s in (sa, sb, sa))
        # a mixed state of both parities on fermions, whose product with another has signs
        ρa = mix(a) + 0.3 * mix(random_state(Pure, sa, 2; rng))
        for (x, y) in ((a, b), (ρa, mix(b)), (a, mix(b)))
            product = TMS.ExactState(x) ⊗ TMS.ExactState(y)
            @test product.state ≈ TMS.ExactState(x ⊗ y).state atol = 1e-14
        end
        fermions = sa[1] isa Fermion
        op = fermions ? dag(C)(1) * C(2) + 0.5 * N(1) : X(1) * Y(2) + Z(1)
        @test inner(TMS.ExactState(a), op, TMS.ExactState(c)) ≈ inner(a, op, c) atol = 1e-14
        g = fermions ? Dissipator(C)(1) : Gate(X(1))
        @test dot(TMS.ExactState(ρa), g, TMS.ExactState(mix(c))) ≈ dot(ρa, g, mix(c)) atol = 1e-14
        # against a State, converted exactly
        @test inner(TMS.ExactState(a), c) ≈ inner(a, c) atol = 1e-14
        @test TMS.fidelity(a, TMS.ExactState(c)) ≈ TMS.fidelity(a, c) atol = 1e-14
        @test TMS.hs_fidelity(TMS.ExactState(ρa), mix(c)) ≈ TMS.hs_fidelity(ρa, mix(c)) atol = 1e-14
    end
    ρ = TMS.ExactState{Mixed}(System(2, Qubit()), "Up")
    @test_throws "is ambiguous" inner(ρ, X(1), ρ)
end

@testset "Sites of an exact state permuted" begin
    # against the State, fermionic signs on the kets and the bras included
    rng = Xoshiro(9)
    for (sys, st) in ((System(4, Fermion(conserve = N)), ["Occ", "Emp", "Occ", "Occ"]),
                      (System([Fermion(), Qubit(), Electron()]), ["Occ", "Up", "Up"]))
        n = length(sys)
        ψ = random_state(State{Pure}(sys, st), 2; rng)
        ρ = mix(ψ) + 0.2 * mix(random_state(State{Pure}(sys, st), 2; rng))
        for p in (collect(n:-1:1), circshift(1:n, 1)), x in (ψ, ρ)
            moved = permute_sites(TMS.ExactState(x), p)
            @test moved.system == permute_sites(x, p).system
            @test moved.state ≈ TMS.ExactState(permute_sites(x, p)).state atol = 1e-13
        end
    end
    @test_throws "a permutation of the positions" permute_sites(
        TMS.ExactState{Pure}(System(2, Qubit()), "Up"), [1, 1])
end

@testset "Measurements against a reference state in a run on exact states" begin
    # the same program on MPS and on exact states, the reference named or given, against
    # each other: product states and gates, exact on both
    sys = System(3, Qubit())
    ref = State{Pure}(sys, "+")
    ms = ["F" => Fidelity(:ref), "O" => Overlap(:ref), "OZ" => Overlap(:ref, Z),
          "OZX" => Overlap(:ref, Z(1) * X(2)), "Fg" => Fidelity(ref), "OXg" => Overlap(ref, X(3))]
    phases = [CreateState(ref), CopyState(:state => :ref),
              Gates(gates = X(1) * H(2), final_measurements = "m" => ms)]
    in_tmp() do
        runTMS(SimData(; name = "mps", phases))
        runTMS(SimData(; name = "exact", backend = TMS.ExactState, phases))
        values(file) = [ parse.(Float64, split(l)[2:end]) for l in readlines(file)[2:end] ]
        @test values("exact/m") ≈ values("mps/m") atol = 1e-12
    end
end

@testset "A product of sums applied factor by factor" begin
    # a product whose factors have several terms each is applied factor by factor, not developed:
    # a chain of 11 CZ on 12 qubits took 1.8 GB applied whole
    q12 = System(12, Qubit())
    x = ExactState(random_state(Pure, q12, 4; rng = Xoshiro(31)))
    cz = prod(CZ(i, i + 1) for i in 1:11)
    one_by_one = foldl((s, i) -> apply(CZ(i, i + 1), s), 11:-1:1; init = x)
    apply(cz, x)
    @test (@allocated apply(cz, x)) < 50 * 2^20
    @test apply(cz, x).state ≈ one_by_one.state
    @test expect(x, cz; complex = true) ≈
          LA.dot(x.state, one_by_one.state) / LA.dot(x.state, x.state)

    # on fermions, against matrices built by hand, pure and mixed
    f5 = System(5, Fermion())
    c(k) = jw_placed(f5, k, matrix(C, Fermion()))
    n(k) = c(k)' * c(k)
    op = (C(1) + dag(C)(2)) * (N(2) + C(3)) * (dag(C)(3) + N(4)) * (C(4) + dag(C)(5))
    m = (c(1) + c(2)') * (n(2) + c(3)) * (c(3)' + n(4)) * (c(4) + c(5)')
    @test isnothing(TMS.developed(f5, TMS.gathered(f5, op)))
    ψ = ExactState(random_state(Pure, f5, 4; rng = Xoshiro(32)))
    v = ψ.state
    @test apply(op, ψ).state ≈ m * v
    @test expect(ψ, op; complex = true) ≈ (v' * m * v) / (v' * v)
    ρ = mix(ψ)
    r = ρ.state
    @test apply(op, ρ).state ≈ m * r * m'
    @test expect(ρ, op; complex = true) ≈ LA.tr(m * r) / LA.tr(r)

    # a product of definite flux whose factors have none, X(1) Swap(1, 2) X(2) conserving N,
    # applied as on a State
    qn = System(3, Qubit(conserve = N))
    xs = ExactState(State{Pure}(qn, ["Up", "Dn", "Up"]))
    @test apply(X(1) * Swap(1, 2) * X(2), xs).state ≈ apply(Swap(1, 2), xs).state
    # on a charged system, a sum of parts changing the charges by different amounts is refused
    fc = System(4, Fermion(conserve = N))
    xc = ExactState(State{Pure}(fc, ["Occ", "Emp", "Occ", "Emp"]))
    @test_throws "no definite flux" apply((C(1) + dag(C)(2)) * (N(3) + C(4)) + C(1), xc)
end

@testset "Channels applied to an exact density matrix by their structure" begin
    # a gate, an evolver and a dissipator act through their operators, never developed, any
    # other channel by its terms: the same as the action of all the terms, on qubits, on
    # fermions with their strings, and on a charged system, whose flux is checked
    byterms(x, a) = TMS.dense_acted(Mixed, TMS.weaken(x.system, ()),
                                    TMS.lower(TMS.weaken(x.system, ()), a), x.state)
    q4 = System(4, Qubit())
    xq = mix(ExactState(random_state(Pure, q4, 4; rng = Xoshiro(3))))
    f4 = System(4, Fermion(conserve = N))
    xf = mix(ExactState(random_state(State{Pure}(f4, ["Occ", "Emp", "Occ", "Emp"]), 4;
                                     rng = Xoshiro(4))))
    hop = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    for (x, a) in [(xq, Gate(CZ(1, 2) * CZ(2, 3) * CZ(3, 4))), (xq, 0.3 * Gate(X(1) + Z(3))),
                   (xq, -im * X(2) * X(3) + Dissipator(Sm(1) * Sm(4))),
                   (xq, depolarizing_channel(0.2)(2) * Gate(H(1))),
                   (xq, Dephase(Z ⊗ Z)(1, 3) + Reset("Up")(4)),
                   (xf, Gate(exp(-0.3im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))(1, 3)) * Gate(C(2) * dag(C)(4))),
                   (xf, -im * hop + Dissipator(C(2)) + Dissipator(N(1) * N(4)))]
        @test apply(a, x).state ≈ byterms(x, a)
        @test inner(x, a, x) ≈ LA.dot(x.state, byterms(x, a))
    end
    @test_throws "no definite flux" apply(Gate(C(1) + dag(C)(2)), xf)
    # but a product of gates whose product of operators has one, C(1) (C(1) + C(2)†) = C(1) C(2)†
    g = Gate(C(1)) * Gate(C(1) + dag(C)(2))
    @test apply(g, xf).state ≈ byterms(xf, g)
    # a gate of a chain of eight CZ, which developed into 4^7 terms on its density matrix
    q8 = System(8, Qubit())
    x8 = mix(ExactState(random_state(Pure, q8, 2; rng = Xoshiro(5))))
    cz = prod(CZ(i, i + 1) for i in 1:7)
    u = matrix(cz, q8)
    @test apply(Gate(cz), x8).state ≈ u * x8.state * u'
    # and the dissipator of a chain of five evolved: K ρ K - ρ, K² = 1, gives
    # e^{-t}(cosh(t) ρ + sinh(t) K ρ K)
    q6 = System(6, Qubit())
    x6 = mix(ExactState(random_state(Pure, q6, 2; rng = Xoshiro(6))))
    cz6 = prod(CZ(i, i + 1) for i in 1:5)
    u6 = matrix(cz6, q6)
    t = 0.3
    r = x6.state
    @test evolve(Dissipator(cz6), t, x6; nsteps = 1).state ≈
        exp(-t) * (cosh(t) * r + sinh(t) * u6 * r * u6')
    # by the structure on fermions, charged or not: P develops into no terms, and a dissipator
    # of three terms is made of the sides of its operator; a run of gates on the charged system
    # whose product has a definite flux where a factor has none
    hop(i, j) = dag(C)(i) * C(j) + dag(C)(j) * C(i)
    P = hop(1, 2) * hop(3, 4) * (N(1) + N(4))
    fu = System(4, Fermion())
    xu = mix(ExactState(random_state(Pure, fu, 4; rng = Xoshiro(7))))
    g3 = Gate(C(1)) * Gate(C(1) + dag(C)(2)) * Gate(N(3) + N(4))
    for (x, a) in [(xu, Gate(P)), (xu, -im * P + Dissipator(P)),
                   (xu, Dissipator(C(1) + C(2) + C(3))), (xf, Gate(P)),
                   (xf, -im * P + Dissipator(P)), (xf, Dissipator(C(1) + C(2) + C(3))), (xf, g3)]
        @test apply(a, x).state ≈ byterms(x, a)
    end
    # and a thermal state of a Hamiltonian that does not develop
    hp = P + dag(P)
    w = exp(-0.7 * matrix(hp, fu))
    _, th = thermalize(hp, 0.7, TMS.ExactState{Mixed}(fu, "FullyMixed"); nsteps = 1)
    @test th.state ≈ w / LA.tr(w) atol = 1e-11
end
