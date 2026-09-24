# Building the basic objects, and package wide quality checks.
#
# Goes here: anything that only checks that an object can be constructed and has the
# announced shape, for System, State, Simulation and RandomState, plus Aqua, and the
# arguments that are refused when an operator meets a system. What a state *measures*
# belongs to sites.jl, not here.

@testset "Aqua" begin
    Aqua.test_all(TensorMixedStates)
end

@testset "System building" begin
    @test_ok System(1, Qubit())
    @test_ok System(3, Qubit())
    @test_ok System([Qubit()])
    @test_ok System([Qubit(), Qubit(), Qubit()])
    # a product keeps the indices of its factors, so a charged system and a plain one cannot
    # be put together, nor two whose quantities could not live on one system
    weak = System(2, Fermion(conserve = N))
    strong_one = System(2, Fermion(conserve = strong(N)))
    @test_throws "carrying charges and one carrying none" weak ⊗ System(2, Qubit())
    @test_throws "strongly on one site and weakly on another" strong_one ⊗ weak
    @test_ok weak ⊗ System(2, Fermion(conserve = N))
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

@testset "Random states" begin
    sys = System(6, Qubit())
    @test maxlinkdim(RandomState{Pure}(sys, 8)) == 8
    # the mixed one goes through a purification, whose link dimension it squares, and has
    # to land on the dimension that was asked for and not on the square below it
    for d in [1, 5, 8, 17, 50]
        @test maxlinkdim(RandomState{Mixed}(sys, d)) == d
    end
    # randomizing an existing state must reach the requested link dimension
    # and must leave the state it was given untouched
    st = State{Pure}(sys, "Up")
    rd = RandomState(st, 8)
    @test maxlinkdim(rd) == 8
    @test norm(rd) ≈ 1
    @test maxlinkdim(st) == 1
    # the element type can be chosen, and defaults to ComplexF64
    @test eltype(RandomState{Pure}(sys, 4).state[1]) == ComplexF64
    @test eltype(RandomState{Pure}(Float64, sys, 4).state[1]) == Float64
    @test eltype(RandomState(Float64, st, 4).state[1]) == Float64
end

@testset "Limits" begin
    # a cutoff is a real number, an integer included, which a field of union type refused
    @test Limits(cutoff = 0).cutoff === 0.0
    @test Limits(cutoff = [0, 1e-10]).cutoff == [0.0, 1e-10]
    # a state taken to zero can still be truncated: a gate of several sites and the sum of two
    # zero states used to fail inside ITensors
    @test Limits(mindim = 0).mindim == 1
    @test Limits(mindim = [0, 2]).mindim == [1, 2]
    q = State{Pure}(System(3, Qubit()), "Up")
    z = apply(Sp(2), q)
    @test norm(z - z) == 0
    @test norm(apply((Sp ⊗ Sp)(1, 2), q)) == 0
    sys = System(6, Qubit())
    st = RandomState{Pure}(sys, 8)
    # `mindim` is the floor the cutoff is not allowed to cross, and `maxdim` keeps the
    # last word when the two ask for opposite things
    @test maxlinkdim(truncate(st; limits = Limits(cutoff = 1e-1))) < 4
    @test maxlinkdim(truncate(st; limits = Limits(cutoff = 1e-1, mindim = 4))) == 4
    @test maxlinkdim(truncate(st; limits = Limits(maxdim = 2, mindim = 4))) == 2
    # the sum of two states goes through the same truncation
    @test maxlinkdim(+(st, st; limits = Limits(cutoff = 1e-1, mindim = 4))) == 4
end

@testset "Simulation building" begin
   @test_ok Simulation(nothing)
   @test_pm Simulation(State{type}(System(3, Qubit()), "Up"))
end

@testset "Putting a state on another system" begin
    # a System draws indices of its own, so the same sites twice give two systems whose
    # states cannot be contracted together. This is what makes them comparable
    for R in (Pure, Mixed)
        sys1 = System(3, Qubit())
        sys2 = System(3, Qubit())
        a = RandomState{R}(sys1, 4)
        b = State(sys2, a)
        @test b.system === sys2
        # the state itself is untouched
        @test expect1(a, Z) ≈ expect1(b, Z)
        @test trace(a) ≈ trace(b)
        @test first(entanglement_entropy(a, 2)) ≈ first(entanglement_entropy(b, 2))
        # and it is comparable with a state of sys2, which a was not
        @test_throws "do not share their System" inner(a, State{R}(sys2, "Up"))
        @test_ok inner(b, State{R}(sys2, "Up"))
    end
    # the sites must match, and a parametric site must match on its parameter too
    sys = System([Qubit(), Boson(4), Qubit()])
    @test_ok State(System([Qubit(), Boson(4), Qubit()]), State{Pure}(sys, "0"))
    @test_throws "system of other sites" State(System(3, Qubit()), State{Pure}(sys, "0"))
    @test_throws "system of other sites" State(System([Qubit(), Boson(5), Qubit()]),
                                               State{Pure}(sys, "0"))
    @test_throws "system of other sites" State(System(2, Qubit()),
                                               State{Pure}(System(3, Qubit()), "Up"))
    # the very same system is a no-op rather than an error
    @test_ok State(sys, State{Pure}(sys, "0"))
end

@testset "Site indices out of the system" begin
    # nothing between writing X(10) and contracting its tensor compares that number with
    # the size of the system, and the three paths an indexed operator can take each reach
    # a different array first, so each used to report a BoundsError on an internal vector
    sys = System(4, Qubit())
    st = State{Pure}(sys, "Up")
    stm = mix(st)

    # the three entries, and both bounds
    for op in (X(5), X(0), X(-1))
        @test_throws "does not have" expect(st, op)
        @test_throws "does not have" measure(st, op)
        @test_throws "does not have" make_mpo(st, op)
        @test_throws "does not have" apply(op, st)
    end
    # the message names the factor at fault and the size of the system
    @test_throws "X(5) acts on site 5" expect(st, X(5))
    @test_throws "it has 4 sites" expect(st, X(5))
    # inside a product and inside a sum, the offending factor is the one named
    @test_throws "Z(9) acts on site 9" expect(st, Z(1) * Z(9))
    @test_throws "Z(7) acts on site 7" make_mpo(st, Z(1) * Z(2) + Z(3) * Z(7))
    # a multi site operator is checked on each of its sites
    @test_throws "Swap(1,9) acts on site 9" apply(Swap(1, 9), st)
    # a time dependent evolver is a vector of terms, each one checked
    @test_throws "acts on site 6" PreMPO(st, [-im * X(1), -im * X(6)])
    # and a mixed representation goes through the same entries
    @test_throws "acts on site 8" expect(stm, X(8))
    @test_throws "acts on site 8" make_mpo(stm, Dissipator(Sm)(8))
    # an operator acting on several sites at once, defined by a matrix or a function of one,
    # has no one site factors to place, and says so rather than failing inside
    m2 = Operator{2}("M2", [1. 0. 0. 0. ; 0. 0. 1. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 1.],
                     involution_op)
    for op in (m2(1, 2), exp(-0.3im * (X ⊗ X))(2, 3))
        @test_throws "acts on several sites at once" expect(st, op)
        @test_throws "acts on several sites at once" make_mpo(st, op)
        @test_ok apply(op, st)
    end
    # partial_trace used to skip a position it did not find when tracing, leaving the state
    # whole, and to raise a BoundsError when keeping it
    @test_throws "given site 5, which the state does not have" partial_trace(stm, [5])
    @test_throws "given site 0" partial_trace(stm, [0, 1]; keepers = true)
    @test_throws "does not have" renyi2(stm, [2, 7])

    # what is inside the system is untouched
    @test_ok expect(st, X(4) * Z(1))
    @test_ok make_mpo(st, sum(Z(i) * Z(i+1) for i in 1:3))
    @test_ok apply(Swap(1, 4), st)

    # a failed measurement used to leave the preobs cache longer than the state, with
    # undefined entries, so a second call on the same state met an UndefRefError rather
    # than the error it deserved. Refusing before the cache is touched settles it
    s2 = State{Pure}(sys, "+")
    @test_throws "does not have" expect(s2, X(9))
    @test length(s2.preobs.left) ≤ length(s2)
    @test_throws "does not have" expect(s2, X(9))
    @test expect(s2, X(3)) ≈ 1
end

@testset "States on a charged system" begin
    IT = TensorMixedStates.ITensors
    sys = System(3, Qubit(conserve = 2Sz))

    # a configuration belongs to one sector, so it builds in either representation
    @test_ok State{Pure}(sys, ["Up", "Dn", "Up"])
    @test_ok State{Mixed}(sys, ["Up", "Dn", "Up"])

    # a density matrix carries no charge whatever sector it lives in, the two ends of
    # |m><m| cancelling in the difference the mixed index holds. This is what lets a state
    # of any particle number be represented, and what forbids a coherence between two
    @test iszero(flux(State{Mixed}(sys, ["Up", "Dn", "Up"]).state))
    @test iszero(flux(State{Mixed}(sys, ["Dn", "Dn", "Dn"]).state))
    @test !iszero(flux(State{Pure}(sys, ["Dn", "Dn", "Dn"]).state))

    # a state spread over two sectors has no charge of its own, and is told so by a message
    # naming it rather than by the `Fluxes not all equal` of ITensors
    @test_throws "spreads over several charges of 2Sz" State{Pure}(sys, "+")
    @test_throws "spreads over several charges of 2Sz" State{Mixed}(sys, "+")
    @test_throws "spreads over several charges" State{Mixed}(sys, [0.5 0.5; 0.5 0.5])

    # a mixture over sectors is representable although a coherence is not: the first is
    # block diagonal, the second is exactly what the charges forbid
    @test_ok State{Mixed}(sys, [0.3 0.; 0. 0.7])

    # a site declaring nothing, among sites that do, takes a trivial charge, and every one
    # of its states stays available, `"+"` included
    m = System([Qubit(conserve = 2Sz), Qubit(), Qubit(conserve = 2Sz)])
    @test_ok State{Pure}(m, ["Up", "+", "Dn"])
    @test_ok State{Mixed}(m, ["Up", "+", "Dn"])
    @test_throws "spreads over several charges" State{Pure}(m, ["+", "Up", "Dn"])
end

@testset "States under a strong symmetry" begin
    IT = TensorMixedStates.ITensors
    strong = TensorMixedStates.strong
    sys = System(3, Fermion(conserve = strong(N)))
    conf = ["Occ", "Emp", "Occ"]

    # the state lives in one sector on the ket side and in the same one on the bra side.
    # This is the rule a pure state already obeys, transposed to a density matrix
    @test flux(State{Mixed}(sys, conf).state) == IT.QN(("N", 2), ("N*", -2))
    @test flux(mix(State{Pure}(sys, conf)).state) == IT.QN(("N", 2), ("N*", -2))

    # a mixture over two sectors has no charge of its own and is refused, where the same
    # state is representable when the quantity is conserved weakly
    @test_throws "spreads over several charges" State{Mixed}(sys, [0.3 0.; 0. 0.7])
    @test_ok State{Mixed}(System(3, Fermion(conserve = N)), [0.3 0.; 0. 0.7])
end

@testset "A random state needs a sector" begin
    strong = TensorMixedStates.strong
    sys = System(4, Fermion(conserve = N))
    conf = ["Occ", "Emp", "Occ", "Emp"]

    # CreateState draws a mixed one from its purification, the only form a charged system
    # allows, where it used to refuse any random mixed state given a `state`
    sim = runTMS(SimData(phases = [CreateState{Mixed}(4, Fermion(conserve = N), conf;
                                                      randomize = 4)]); output = devnull)
    @test sim.state isa State{Mixed}
    @test maxlinkdim(sim.state) ≤ 4

    # there is no sector to draw a pure state in, so the form taking a system alone is
    # refused by a message naming the one that works
    @test_throws "no sector to draw it in" RandomState{Pure}(sys, 4)
    @test_throws "RandomState(state, linkdims)" RandomState{Pure}(sys, 4)

    # randomising a state one already has keeps it in the sector it was in
    p = State{Pure}(sys, conf)
    @test flux(RandomState(p, 4).state) == flux(p.state)

    # a mixed one is drawn from the states its purification starts from. Tracing half of
    # that purification leaves a genuine mixture over the sectors around the one named,
    # which is what a weak symmetry allows and a pure state cannot be
    @test_throws "name the states its purification starts from" RandomState{Mixed}(sys, 4)
    m = RandomState{Mixed}(sys, conf, 16)
    @test trace(m) ≈ 1
    @test hermiticity(m) ≈ 1
    @test real(trace2(m)) < 1
    @test flux(m.state) == flux(mix(p).state)

    # conserving strongly leaves no room for it, what the partial trace leaves spreading
    # over several sectors
    @test_throws "cannot hold" RandomState{Mixed}(System(4, Fermion(conserve = strong(N))),
                                                  conf, 16)
end

@testset "The ladder of symmetries" begin
    # weakening walks strong, then weak, then nothing, and stops there. Every rung holds the
    # same physics, so every rung must give the numbers a system conserving nothing gives —
    # including a four point correlator, which is what a chain sized on a guess got wrong
    conf = ["Occ", "Emp", "Occ", "Emp"]
    build(site) = begin
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        s = p(conf) + 0.5 * p(reverse(conf))
        return mix(s / norm(s))
    end
    numbers(s) = [ real(trace(s)), real(expect(s, N(2))), real(trace2(s)),
                   real(hermiticity(s)),
                   real(expect(s, dag(C)(1) * dag(C)(2) * C(3) * C(4))) ]
    reference = numbers(build(Fermion()))

    x = build(Fermion(conserve = strong(N)))
    @test symmetries(x.system) == TensorMixedStates.Conserved([("N", true)])
    @test numbers(x) ≈ reference

    x = weaken(x)
    @test symmetries(x.system) == TensorMixedStates.Conserved([("N", false)])
    @test numbers(x) ≈ reference

    x = weaken(x)
    @test !TensorMixedStates.is_charged(x.system)
    @test numbers(x) ≈ reference

    # nothing left to weaken gives the state back, so it is always safe to call
    @test weaken(x) === x

    # a target says what must still be conserved, in the vocabulary `conserve` takes
    y = build(Fermion(conserve = strong(N)))
    @test weaken(y, symmetries(y.system)) === y
    @test numbers(weaken(y, ())) ≈ reference
    @test numbers(weaken(y, N)) ≈ reference

    # a quantity may be dropped or asked for less strongly, never invented nor strengthened
    w = build(Fermion(conserve = N))
    @test_throws "cannot make it strong" weaken(w, strong(N))
    @test_throws "cannot start conserving" weaken(w, Ntot)
    @test_throws "cannot make it strong" weaken(w.system, strong(N))
    @test_throws "cannot start conserving" weaken(w.system, Ntot)

    # and a system reports what it conserves in the form the target takes
    @test repr(symmetries(System(2, Electron(conserve = (strong(Ntot), 2Sz))))) ==
        "(strong(Ntot), 2Sz)"
    @test repr(symmetries(System(2, Fermion()))) == "()"

    # a site conserving nothing, put first, used to hide what the others conserve: weakening
    # then gave the state back untouched, still charged, and refused `N` as a target
    mixed_sites(site) = begin
        sys = System([Qubit(); fill(site, 4)])
        p(v) = State{Pure}(sys, ["Up"; v])
        s = p(conf) + 0.5 * p(reverse(conf))
        return mix(s / norm(s))
    end
    shifted(s) = [ real(trace(s)), real(expect(s, N(3))), real(trace2(s)),
                   real(expect(s, dag(C)(2) * dag(C)(3) * C(4) * C(5))) ]
    reference = shifted(mixed_sites(Fermion()))

    z = mixed_sites(Fermion(conserve = strong(N)))
    @test repr(symmetries(z.system)) == "strong(N)"
    @test shifted(weaken(z, N)) ≈ reference
    z = weaken(z)
    @test repr(symmetries(z.system)) == "N"
    @test shifted(z) ≈ reference
    z = weaken(z)
    @test !TensorMixedStates.is_charged(z.system)
    @test shifted(z) ≈ reference
end

