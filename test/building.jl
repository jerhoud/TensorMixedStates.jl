# Building the basic objects: System, State and random states, and the limits of a
# truncation, plus Aqua. What a state measures belongs to sites.jl.
#
# The phases that make the state, run by runTMS, are in runner.jl; the noise and the Krylov
# parameters of the solvers in searches.jl.

@testset "Aqua" begin
    Aqua.test_all(TensorMixedStates)
end

@testset "System building" begin
    @test_ok System(1, Qubit())
    @test_ok System(3, Qubit())
    @test_ok System([Qubit()])
    @test_ok System([Qubit(), Qubit(), Qubit()])
    # a product is the system of its sites, a charged system and a plain one included; two
    # whose quantities could not live on one system are refused
    weak = System(2, Fermion(conserve = N))
    strong_one = System(2, Fermion(conserve = strong(N)))
    @test weak ⊗ System(2, Qubit()) == System([fill(Fermion(conserve = N), 2); fill(Qubit(), 2)])
    @test_throws "strongly on one site and weakly on another" strong_one ⊗ weak
    @test_ok weak ⊗ System(2, Fermion(conserve = N))
end

@testset "Systems of the same sites are one system" begin
    # the indices derive from the content, so two systems of the same sites have the same
    # indices, and the halves of S ⊗ S differ by their position
    @test System(3, Qubit()) == System(3, Qubit())
    @test hash(System(3, Qubit())) == hash(System(3, Qubit()))
    @test siteinds(Pure, System(3, Qubit())) == siteinds(Pure, System(3, Qubit()))
    @test siteinds(Mixed, System(3, Qubit())) == siteinds(Mixed, System(3, Qubit()))
    @test System(3, Qubit()) ≠ System(3, Qubit(conserve = N))
    @test System([Boson(4)]) ≠ System([Boson(5)])
    # a system keeps its own copy of the sites it is given
    v = TMS.AbstractSite[Qubit(), Qubit()]
    s = System(v)
    v[1] = Boson(3)
    @test s[1] == Qubit()
    # and its key, from which the indices take their identities, reads every site: Julia's hash
    # of a long vector skips some, so that two systems shared their indices
    a = [isodd(i) ? Qubit() : Boson(3) for i in 1:100000]
    b = copy(a)
    b[50000] = Boson(4)
    @test System(a).key ≠ System(b).key
    S, U = System(2, Qubit()), System(1, Qubit())
    for sys in (S ⊗ S ⊗ S, S ⊗ U ⊗ S, S ⊗ (S ⊗ S))
        @test allunique(siteinds(Pure, sys))
        @test allunique(siteinds(Mixed, sys))
        @test expect(apply(X(3), State{Pure}(sys, "Up")), Z(3)) ≈ -1
    end
end

@testset "State building" begin
    @test_pm State{type}(System(1, Qubit()), "Up")
    @test_pm State{type}(3, Qubit(), "Up")
    @test_pm State{type}([Qubit(), Fermion(), Boson(4)], ["Up", "Occ", "2"])
    @test_pm State{type}(System(3, Qubit()), [1., 0.])
    @test_ok State{Pure}(System(3, Qubit()), [[1., 0.], [0, 1], [0, im]])
    @test_ok State{Mixed}(System(3, Qubit()), [[1 0 ; 0 0 ], [1 1; 1 1], [0 0 ; 0 1]])
    @test_ok State{Mixed}(System(3, Qubit()), "FullyMixed")
end

@testset "A sum of states keeps its bond dimension" begin
    # the eigenvalues rounding leaves were kept as states by the sum: nine product states of
    # twenty qubits had a bond dimension of 2304
    rng = Xoshiro(3)
    sys = System(20, Qubit())
    cfgs = [ rand(rng, ["Up", "Dn"], 20) for _ in 1:12 ]
    cs = randn(rng, 12)
    s = sum(cs[k] * State{Pure}(sys, cfgs[k]) for k in 1:12)
    @test maxlinkdim(s) ≤ 12
    @test [ inner(State{Pure}(sys, c), s) for c in cfgs ] ≈ cs
    @test maxlinkdim(s - cs[1] * State{Pure}(sys, cfgs[1])) ≤ 13
end

@testset "A difference sees what the two states differ by" begin
    # a direct sum, then a truncation judged on the result. The sum of ITensorMPS
    # cut on the density matrices of the terms, and norm(a - b) saw no relative difference
    # below about 1e-7
    rng = Xoshiro(1)
    sys = System(20, Qubit())
    a = random_state(Pure, sys, 8; rng)
    c = random_state(Pure, sys, 8; rng)
    for δ in (1e-6, 1e-8, 1e-10, 1e-12)
        b = +(a, δ * c; limits = Limits(cutoff = 0))
        @test norm(a - b) ≈ δ rtol = 1e-3
    end
    # and the difference of two equal states on charged sites is the null state, which ITensors
    # cannot truncate, its tensors having no block left
    fermions = State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Emp", "Occ", "Emp"])
    for s in (State{Pure}(System(3, Qubit(conserve = N)), ["Up", "Dn", "Up"]),
              random_state(fermions, 4))
        @test norm(s - s) < 1e-12
        @test norm(s + (-1) * s) < 1e-12
        @test norm(mix(s) - mix(s)) < 1e-12
        @test norm(s - 1.0000001 * s) ≈ 1e-7 * norm(s) rtol = 1e-6
    end
end

@testset "A sum of states of one site" begin
    # the direct sum needs links, which a state of one site does not have
    sys = System(1, Qubit())
    u, d = State{Pure}(sys, "Up"), State{Pure}(sys, "Dn")
    @test norm(u + u) ≈ 2
    @test norm(u - u) < 1e-12
    @test norm(mix(u) - mix(u)) < 1e-12
    @test real(expect(normalize(u + d), X(1))) ≈ 1
    c = State{Pure}(System(1, Fermion(conserve = N)), "Occ")
    @test norm(c - c) < 1e-12
end

@testset "The default limits cut the noise of rounding" begin
    rng = Xoshiro(3)
    sys = System(20, Qubit())
    s = sum(randn(rng) * State{Pure}(sys, rand(rng, ["Up", "Dn"], 20)) for _ in 1:9)
    # a field left nothing takes the default of the function, 1e-26 for a sum
    @test Limits().cutoff === nothing
    @test maxlinkdim(apply(prod(controlled(Z)(i, i + 1) for i in 1:19), s)) == maxlinkdim(s)
    @test maxlinkdim(evolve(-im * sum(Z(i) for i in 1:20), 0.1, s; nsteps = 1)) == maxlinkdim(s)
    # a weight of 1e-28 is below the default cutoff, one of 1e-18 is not, and a cutoff of zero
    # keeps everything
    up, dn = State{Pure}(System(6, Qubit()), "Up"), State{Pure}(System(6, Qubit()), "Dn")
    @test maxlinkdim(up + 1e-14 * dn) == 1
    @test maxlinkdim(up + 1e-9 * dn) == 2
    @test maxlinkdim(+(up, 1e-14 * dn; limits = Limits(cutoff = 0))) ≥ 2
    # an MPO is applied by "naive" by default, cut at 1e-26, and by "densitymatrix" at eps(),
    # below which its density matrices are rounding
    q2 = State{Pure}(System(2, Qubit()), "Up")
    mpo = make_mpo(q2, Id(1) + 1e-9 * X(1) * X(2))
    @test maxlinkdim(apply(mpo, q2)) == 2
    @test maxlinkdim(apply(mpo, q2; apply_algo = "densitymatrix")) == 1
end

@testset "Random states" begin
    sys = System(6, Qubit())
    @test maxlinkdim(random_state(Pure, sys, 8)) == 8
    # the mixed one goes through a purification, whose link dimension it squares, and nothing
    # is truncated
    for d in [1, 5, 8, 17, 50]
        @test maxlinkdim(random_state(Mixed, sys, d)) == isqrt(d)^2
    end
    # randomizing an existing state reaches the requested link dimension and leaves the state
    # it was given untouched
    st = State{Pure}(sys, "Up")
    rd = random_state(st, 8)
    @test maxlinkdim(rd) == 8
    @test norm(rd) ≈ 1
    @test maxlinkdim(st) == 1
    # the element type can be chosen, and defaults to ComplexF64
    @test eltype(random_state(Pure, sys, 4).state[1]) == ComplexF64
    @test eltype(random_state(Pure, sys, 4; eltype = Float64).state[1]) == Float64
    @test eltype(random_state(st, 4; eltype = Float64).state[1]) == Float64
    # limits holding vectors, for the steps of an evolution or the sweeps of a search, are
    # refused elsewhere, where ITensors raised a MethodError
    @test_throws "are for the evolutions and searches only" truncate(
        random_state(Pure, System(4, Qubit()), 4); limits = Limits(cutoff = [1e-10, 1e-12]))
    # a complex state is refused clearly, where ITensors raised an InexactError
    @test_throws "of complex amplitudes is not randomised" random_state(
        State{Pure}(System(3, Qubit()), [0.6, 0.8im]), 2; eltype = Float64)
    # a generator of one's own draws the same state twice
    rng(k) = Xoshiro(k)
    for R in (Pure, Mixed)
        draw() = random_state(R, sys, 4; rng = rng(5))
        @test norm(draw() - draw()) < 1e-14
    end
end

@testset "A random mixed state is a state of its system" begin
    # traced from its purification with nothing truncated, it has a trace of one and no
    # negative eigenvalue, and lives on the system it was drawn for
    sys = System(4, Qubit())
    for d in [4, 10, 12]
        ρ = random_state(Mixed, sys, d)
        @test ρ.system == sys
        m = ket_bra_matrix(ρ)
        @test LA.tr(m) ≈ 1
        @test minimum(LA.eigvals(LA.Hermitian((m + m') / 2))) > -1e-12
        @test_ok fidelity(ρ, State{Pure}(sys, "Up"))
    end
    sq = System(4, Fermion(conserve = N))
    conf = ["Occ", "Emp", "Occ", "Emp"]
    ρq = random_state(Mixed, sq, conf, 4)
    @test ρq.system == sq
    @test_ok fidelity(ρq, State{Pure}(sq, conf))
    # the local state of the purification may be given once, as for State
    @test trace(random_state(Mixed, System(3, Qubit()), [1., 0.], 4)) ≈ 1
    @test trace(random_state(Mixed, System(3, Qubit()), Level => 1, 4)) ≈ 1
    # a vector of local states of the wrong length is refused with the sizes of the system
    @test_throws "system (3) and states (2)" random_state(
        Mixed, System(3, Qubit()), ["Up", "Dn"], 4)
    # randomizing a state of one site leaves the one it was given alone: its tensor was
    # written in place when the element type already matched, and its cache went stale
    for (elt, s) in [(Float64, State{Pure}(System(1, Qubit()), "Up")),
                     (ComplexF64, random_state(Pure, System(1, Qubit()), 1))]
        t = copy(s.state[1])
        z = expect(s, Z(1))
        random_state(s, 2; eltype = elt)
        @test norm(s.state[1] - t) < 1e-14
        @test expect(s, Z(1)) ≈ z
    end
end

@testset "A number is not a local state" begin
    # a basis state is named, "1" or Level => 1, a number standing for nothing
    sys = System(3, Qubit())
    @test norm(State{Pure}(sys, Level => 1) - State{Pure}(sys, [0., 1.])) < 1e-12
    @test norm(State{Mixed}(sys, Level => 1) - State{Mixed}(sys, [0., 1.])) < 1e-12
    @test norm(State{Pure}(sys, "1") - State{Pure}(sys, [0., 1.])) < 1e-12
    @test_throws "a number does not denote a state" State{Pure}(sys, 1)
    @test_throws "a number does not denote a state" State{Pure}(sys, Any[0, 1, 0])
    @test norm(State{Pure}(sys, Level .=> [0, 1, 0]) - State{Pure}(sys, ["Up", "Dn", "Up"])) < 1e-12
end

@testset "Limits" begin
    # a cutoff is a real number, an integer included
    @test Limits(cutoff = 0).cutoff === 0.0
    @test Limits(cutoff = [0, 1e-10]).cutoff == [0.0, 1e-10]
    # a state taken to zero can still be truncated
    @test Limits(mindim = 0).mindim == 1
    @test Limits(mindim = [0, 2]).mindim == [1, 2]
    # a dimension is an integer of any type, and a float is refused by name
    @test Limits(maxdim = Int32(10)).maxdim === 10
    @test Limits(maxdim = Int32[10, 20], mindim = UInt8(2)).maxdim == [10, 20]
    @test_throws "maxdim is an integer" Limits(maxdim = 1e3)
    @test_throws "mindim is an integer" Limits(mindim = [1.5])
    # a field left nothing takes the default of the function given the limits, field by field
    @test Limits() == Limits(cutoff = nothing, maxdim = nothing, mindim = nothing)
    @test TMS.truncation(Limits(maxdim = 5)) ==
          (cutoff = TMS.rounding_tol^2, maxdim = 5, mindim = 1)
    @test TMS.truncation(Limits(maxdim = 5); cutoff = 0.).cutoff === 0.0
    @test TMS.truncation(Limits(cutoff = 1e-8); cutoff = 0.).cutoff === 1e-8
    # superposition is exact by default, where a cutoff of eps() given, which the limits by
    # default could not be told from, cuts a weight of 1e-18
    s4, terms = System(4, Qubit()), [1 => "Up", 1e-9 => "Dn"]
    @test maxlinkdim(superposition(s4, terms)) == 2
    @test maxlinkdim(superposition(s4, terms; limits = Limits(cutoff = eps()))) == 1
    @test maxlinkdim(superposition(s4, terms; limits = Limits(maxdim = 1))) == 1
    q = State{Pure}(System(3, Qubit()), "Up")
    z = 0 * q
    @test norm(z - z) == 0
    @test norm(apply((Sp ⊗ Sp)(1, 2), q)) == 0
    sys = System(6, Qubit())
    st = random_state(Pure, sys, 8)
    # `mindim` is the floor the cutoff is not allowed to cross, and `maxdim` keeps the last
    # word when the two ask for opposite things
    @test maxlinkdim(truncate(st; limits = Limits(cutoff = 1e-1))) < 4
    @test maxlinkdim(truncate(st; limits = Limits(cutoff = 1e-1, mindim = 4))) == 4
    @test maxlinkdim(truncate(st; limits = Limits(maxdim = 2, mindim = 4))) == 2
    # the sum of two states goes through the same truncation
    @test maxlinkdim(+(st, st; limits = Limits(cutoff = 1e-1, mindim = 4))) == 4
end

@testset "States of systems of the same sites compare" begin
    # a state built on a system, saved and read back, or built on another system of the same
    # sites, contracts with the others; one of other sites is refused
    for R in (Pure, Mixed)
        a = random_state(R, System(3, Qubit()), 4)
        b = State{R}(System(3, Qubit()), "Up")
        @test_ok inner(a, b)
        @test inner(b, State{R}(System(3, Qubit()), "Up")) ≈ 1
        @test_throws "on different systems" inner(a, State{R}(System(3, Qubit(conserve = N)), "Up"))
    end
    @test_throws "on different systems" inner(State{Pure}(System([Qubit(), Boson(4)]), "0"),
                                              State{Pure}(System([Qubit(), Boson(5)]), "0"))
end

@testset "An operator of several sites on a repeated site" begin
    for op in (() -> Swap(1, 1), () -> (X ⊗ Y)(2, 2), () -> controlled(X)(3, 3),
               () -> (X ⊗ Id ⊗ Z)(1, 2, 1))
        @test_throws "repeats a site" op()
    end
    @test_ok (X ⊗ Y)(1, 2)
end

@testset "States on a charged system" begin
    sys = System(3, Qubit(conserve = 2Sz))
    # a configuration belongs to one sector, so it builds in either representation
    @test_ok State{Pure}(sys, ["Up", "Dn", "Up"])
    @test_ok State{Mixed}(sys, ["Up", "Dn", "Up"])
    # a density matrix carries no charge whatever sector it lives in
    @test iszero(flux(State{Mixed}(sys, ["Up", "Dn", "Up"]).state))
    @test iszero(flux(State{Mixed}(sys, ["Dn", "Dn", "Dn"]).state))
    @test !iszero(flux(State{Pure}(sys, ["Dn", "Dn", "Dn"]).state))
    # a state spread over two sectors has no charge of its own, and is told so
    @test_throws "spreads over several charges of 2Sz" State{Pure}(sys, "+")
    @test_throws "spreads over several charges of 2Sz" State{Mixed}(sys, "+")
    @test_throws "spreads over several charges" State{Mixed}(sys, [0.5 0.5; 0.5 0.5])
    # a mixture over sectors is representable although a coherence is not
    @test_ok State{Mixed}(sys, [0.3 0.; 0. 0.7])
    # a site declaring nothing, among sites that do, takes a trivial charge, and every one of
    # its states stays available
    m = System([Qubit(conserve = 2Sz), Qubit(), Qubit(conserve = 2Sz)])
    @test_ok State{Pure}(m, ["Up", "+", "Dn"])
    @test_ok State{Mixed}(m, ["Up", "+", "Dn"])
    @test_throws "spreads over several charges" State{Pure}(m, ["+", "Up", "Dn"])
end

@testset "States under a strong symmetry" begin
    sys = System(3, Fermion(conserve = strong(N)))
    conf = ["Occ", "Emp", "Occ"]
    # the state lives in one sector on the ket side and in the same one on the bra side
    @test flux(State{Mixed}(sys, conf).state) == IT.QN(("N", 2), ("N*", -2))
    @test flux(mix(State{Pure}(sys, conf)).state) == IT.QN(("N", 2), ("N*", -2))
    # a mixture over two sectors has no charge of its own and is refused, where the same state
    # is representable when the quantity is conserved weakly
    @test_throws "spreads over several charges" State{Mixed}(sys, [0.3 0.; 0. 0.7])
    @test_ok State{Mixed}(System(3, Fermion(conserve = N)), [0.3 0.; 0. 0.7])
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

@testset "The trace is kept with the state" begin
    # computed once, under the lock of the state, with the weak form a state conserving
    # something strongly is measured on, whose system its system makes once for all its states
    ρ = mix(State{Pure}(System(3, Fermion(conserve = strong(N))), ["Occ", "Emp", "Occ"]))
    @test isempty(ρ.preobs.trace) && isempty(ρ.preobs.weak)
    @test trace(ρ) ≈ 1
    @test length(ρ.preobs.trace) == 1
    @test only(ρ.preobs.weak).system === weaken(ρ).system
    @test weaken(mix(State{Pure}(ρ.system, ["Emp", "Occ", "Occ"]))).system === weaken(ρ).system
    ψ = random_state(Pure, System(4, Qubit()), 4)
    @test trace(ψ) ≈ norm(ψ)^2
    @test only(ψ.preobs.trace) === trace(ψ)
    # a state built from it starts with caches of its own
    @test isempty(normalize(ψ).preobs.trace)
    # tr of LinearAlgebra is the trace, where it asked for the size of the state
    @test tr(ψ) === trace(ψ) && tr(ρ) === trace(ρ)
end

@testset "A random state needs a sector" begin
    sys = System(4, Fermion(conserve = N))
    conf = ["Occ", "Emp", "Occ", "Emp"]
    # there is no sector to draw a pure state in, so the form taking a system alone is
    # refused by a message naming the one that works
    @test_throws "no sector to draw it in" random_state(Pure, sys, 4)
    @test_throws "random_state(state, linkdims)" random_state(Pure, sys, 4)
    # randomising a state one already has keeps it in the sector it was in
    p = State{Pure}(sys, conf)
    @test flux(random_state(p, 4).state) == flux(p.state)
    # a mixed one is drawn from the states its purification starts from
    @test_throws "name the local states its purification starts from" random_state(Mixed, sys, 4)
    m = random_state(Mixed, sys, conf, 16)
    @test trace(m) ≈ 1
    @test hermiticity(m) ≈ 1
    @test real(purity(m)) < 1
    @test flux(m.state) == flux(mix(p).state)
    # conserving strongly leaves no room for it
    strongly(n) = System(n, Fermion(conserve = strong(N)))
    @test_throws "cannot hold" random_state(Mixed, strongly(4), conf, 16)
    @test_throws "conserving something strongly" random_state(Mixed, strongly(2), 4)
end

@testset "The ladder of symmetries" begin
    # weakening walks strong, then weak, then nothing, and stops there. Every rung holds the
    # same physics, so every rung must give the numbers a system conserving nothing gives,
    # including a four point correlator, which is what a chain sized on a guess got wrong
    conf = ["Occ", "Emp", "Occ", "Emp"]
    build(site) = begin
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        s = p(conf) + 0.5 * p(reverse(conf))
        return mix(s / norm(s))
    end
    numbers(s) = [ real(trace(s)), real(expect(s, N(2))),
                   real(expect(s, dag(C)(1) * dag(C)(2) * C(3) * C(4))) ]
    purities(s) = [ real(purity(s)), real(hermiticity(s)) ]
    reference = numbers(build(Fermion()))

    x = build(Fermion(conserve = strong(N)))
    @test symmetries(x.system) == TensorMixedStates.Conserved([("N", true)])
    @test numbers(x) ≈ reference
    @test purities(x) ≈ purities(build(Fermion()))

    x = weaken(x)
    @test symmetries(x.system) == TensorMixedStates.Conserved([("N", false)])
    @test numbers(x) ≈ reference
    @test purities(x) ≈ purities(build(Fermion()))

    x = weaken(x)
    @test !TensorMixedStates.is_charged(x.system)
    @test numbers(x) ≈ reference
    @test purities(x) ≈ purities(build(Fermion()))
    @test norm(x - build(Fermion())) < 1e-12

    # nothing left to weaken gives the state back, so it is always safe to call
    @test weaken(x) === x

    # a target says what must still be conserved, in the vocabulary `conserve` takes
    y = build(Fermion(conserve = strong(N)))
    @test weaken(y, symmetries(y.system)) === y
    @test norm(weaken(y, ()) - build(Fermion())) < 1e-12
    @test numbers(weaken(y, N)) ≈ reference
    @test norm(weaken(weaken(y, N), ()) - build(Fermion())) < 1e-12

    # a quantity may be dropped or asked for less strongly, never invented nor strengthened
    w = build(Fermion(conserve = N))
    @test_throws "cannot make it strong" weaken(w, strong(N))
    @test_throws "cannot start conserving" weaken(w, Ntot)
    @test_throws "cannot make it strong" weaken(w.system, strong(N))
    @test_throws "cannot start conserving" weaken(w.system, Ntot)

    # and a system reports what it conserves in the form the target takes
    @test repr(symmetries(System(2, Electron(conserve = (strong(Ntot), 2Sz))))) ==
        "(2Sz, strong(Ntot))"
    @test repr(symmetries(System(2, Fermion()))) == "()"

    # a site conserving nothing, put first, used to hide what the others conserve: weakening
    # then gave the state back untouched, still charged, and refused `N` as a target
    mixed_sites(site) = begin
        sys = System([Qubit(); fill(site, 4)])
        p(v) = State{Pure}(sys, ["Up"; v])
        s = p(conf) + 0.5 * p(reverse(conf))
        return mix(s / norm(s))
    end
    shifted(s) = [ real(trace(s)), real(expect(s, N(3))),
                   real(expect(s, dag(C)(2) * dag(C)(3) * C(4) * C(5))) ]
    plain = mixed_sites(Fermion())
    reference = shifted(plain)
    z = mixed_sites(Fermion(conserve = strong(N)))
    @test repr(symmetries(z.system)) == "strong(N)"
    @test shifted(weaken(z, N)) ≈ reference
    @test norm(weaken(weaken(z, N), ()) - plain) < 1e-12
    z = weaken(z)
    @test repr(symmetries(z.system)) == "N"
    @test shifted(z) ≈ reference
    z = weaken(z)
    @test !TensorMixedStates.is_charged(z.system)
    @test shifted(z) ≈ reference
    @test real(purity(z)) ≈ real(purity(plain))
    @test norm(z - plain) < 1e-12
end

@testset "A pure state dropping part of what it conserves" begin
    # a site left conserving nothing is cut into a single block on a system still charged,
    # where the relabelled index keeps one block per basis state
    sys = System([Qubit(conserve = N), Boson(3, conserve = named(N, "Nb")),
                  Boson(3, conserve = named(N, "Nb")), Qubit(conserve = N)])
    p(v) = State{Pure}(sys, v)
    ψ = normalize(p(["1", "2", "0", "0"]) + 0.5 * p(["0", "1", "1", "1"]) -
                  0.7im * p(["1", "0", "2", "0"]))
    plain = weaken(ψ, ())
    for target in (N, named(N, "Nb"))
        w = weaken(ψ, target)
        @test symmetries(w.system) == symmetries(weaken(sys, target))
        @test norm(weaken(w, ()) - plain) < 1e-12
        ρw = weaken(mix(ψ), target)
        @test norm(weaken(ρw, ()) - mix(plain)) < 1e-12
        @test fidelity(ρw, w) ≈ 1
    end
end
