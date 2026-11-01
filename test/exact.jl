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

@testset "Measurements on $name" for (name, sys, ops) in [
        ("qubits", System(4, Qubit()), [X(1) * Z(3) + 2 * Y(2), Sp(1) * Sm(4)]),
        ("fermions", System(4, Fermion()), [dag(C)(1) * C(3) + dag(C)(3) * C(1), C(2), N(4)]),
        ("a qubit among fermions", System([Fermion(), Qubit(), Electron(), Fermion()]),
         [dag(C)(1) * Cup(3) + dag(Cup)(3) * C(1), C(4), X(2) * Ntot(3)])]
    rng = Xoshiro(20261103)
    ψ, ρ = random_state(Pure, sys, 4; rng), random_state(Mixed, sys, 4; rng)
    eψ, eρ = TMS.ExactState(ψ), TMS.ExactState(ρ)
    @test eψ.state ≈ raw_vector(ψ)
    @test eρ.state ≈ ket_bra_matrix(ρ)
    for (e, st) in [(eψ, ψ), (eρ, ρ)]
        for o in ops
            @test expect(e, o; complex = true) ≈ expect(st, o; complex = true) atol = 1e-12
            @test measure(e, o + dag(o)) ≈ measure(st, o + dag(o)) atol = 1e-12
        end
        @test expect1(e, F) ≈ expect1(st, F) atol = 1e-12
        @test trace(e) ≈ trace(st)
        @test norm(e) ≈ norm(st)
    end
    # unnormalised, measured through their trace
    for st in [2 * ψ, 3 * ρ]
        e = TMS.ExactState(st)
        @test trace(e) ≈ trace(st)
        @test expect(e, first(ops); complex = true) ≈ expect(st, first(ops); complex = true) atol =
              1e-12
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
    # a strong symmetry, whose mixed index full_matrix reads without weakening it
    ss = System(3, Fermion(conserve = strong(N)))
    ρs = mix(random_state(State{Pure}(ss, ["Occ", "Emp", "Occ"]), 3))
    es = TMS.ExactState(ρs)
    @test es.state ≈ ket_bra_matrix(ρs)
    hop = dag(C)(1) * C(2) + dag(C)(2) * C(1)
    @test expect(es, hop) ≈ expect(ρs, hop) atol = 1e-12
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
    @test_throws "beyond the 4096 TMS builds" TMS.ExactState(State{Mixed}(System(13, Qubit()), "Up"))
    @test_throws "beyond the 16777216 TMS builds" TMS.ExactState(System(64, Qubit()), ComplexF64[])
    @test_throws "divides by the trace" TMS.ExactState(q2, zeros(4))
end

@testset "Refusals of the measurements and conversion back" begin
    ψ = State{Pure}(System(2, Qubit()), "Up")
    @test_throws "expect divides by the trace" expect(TMS.ExactState(0 * ψ), Z(1))
    e = TMS.ExactState(ψ)
    @test_throws "is a generator" expect(e, Dissipator(Sm)(1))
    @test_throws "placed on sites" expect(e, X)
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
