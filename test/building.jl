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

@testset "A product of systems has an index per site" begin
    # S ⊗ U ⊗ S shares indices with its left part without being the same object, and the right
    # operand is then renewed as in S ⊗ S: a gate could not find its site otherwise
    S, U = System(2, Qubit()), System(1, Qubit())
    for sys in (S ⊗ S ⊗ S, S ⊗ U ⊗ S, S ⊗ (S ⊗ S))
        @test allunique(sys.pure_indices)
        @test allunique(sys.mixed_indices)
        @test expect(apply(X(3), State{Pure}(sys, "Up")), Z(3)) ≈ -1
    end
    # disjoint systems keep their indices
    @test (S ⊗ U).pure_indices == [S.pure_indices; U.pure_indices]
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
    rng = TensorMixedStates.Random.Xoshiro(3)
    sys = System(20, Qubit())
    cfgs = [ rand(rng, ["Up", "Dn"], 20) for _ in 1:12 ]
    cs = randn(rng, 12)
    s = sum(cs[k] * State{Pure}(sys, cfgs[k]) for k in 1:12)
    @test maxlinkdim(s) ≤ 12
    @test [ inner(State{Pure}(sys, c), s) for c in cfgs ] ≈ cs
    # what the truncation of the sum left of the first term, of a weight of the order of the
    # cutoff, may or may not be cut: the bound is that of a sum of states of dimensions 12 and 1
    @test maxlinkdim(s - cs[1] * State{Pure}(sys, cfgs[1])) ≤ 13
end

@testset "The default limits cut the noise of rounding" begin
    # gates and an evolution that create no entanglement kept the singular values of rounding
    rng = TensorMixedStates.Random.Xoshiro(3)
    sys = System(20, Qubit())
    s = sum(randn(rng) * State{Pure}(sys, rand(rng, ["Up", "Dn"], 20)) for _ in 1:9)
    @test Limits().cutoff == eps()
    @test maxlinkdim(apply(prod(controlled(Z)(i, i + 1) for i in 1:19), s)) == maxlinkdim(s)
    @test maxlinkdim(tdvp(-im * sum(Z(i) for i in 1:20), 0.1, s)) == maxlinkdim(s)
    # a weight of 1e-18 is below the default cutoff, and a cutoff of zero keeps everything
    up, dn = State{Pure}(System(6, Qubit()), "Up"), State{Pure}(System(6, Qubit()), "Dn")
    dn = State(up.system, dn)
    @test maxlinkdim(up + 1e-9 * dn) == 1
    @test maxlinkdim(+(up, 1e-9 * dn; limits = Limits(cutoff = 0))) ≥ 2
end

@testset "Random states" begin
    sys = System(6, Qubit())
    @test maxlinkdim(RandomState{Pure}(sys, 8)) == 8
    # the mixed one goes through a purification, whose link dimension it squares, and nothing
    # is truncated: it lands on the square not above the dimension asked for, a state the
    # first truncation of an evolution to that dimension keeps whole
    for d in [1, 5, 8, 17, 50]
        @test maxlinkdim(RandomState{Mixed}(sys, d)) == isqrt(d)^2
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

@testset "A random mixed state is a state of its system" begin
    # traced from its purification with nothing truncated, it has a trace of one and no
    # negative eigenvalue, where truncating to a dimension that is not a square cut into a
    # spectrum with no small tail. It lives on the system it was drawn for, which inner and
    # the fidelities require
    LA = TensorMixedStates.LinearAlgebra
    function dense(ρ)
        n = length(ρ.system)
        v = Array(prod(ρ.state), reverse(ρ.system.mixed_indices)...)
        perm = vcat([2k - 1 for k in 1:n], [2k for k in 1:n])
        return reshape(permutedims(reshape(v, ntuple(_ -> 2, 2n)...), perm), 2^n, 2^n)
    end
    sys = System(4, Qubit())
    for d in [4, 10, 12]
        ρ = RandomState{Mixed}(sys, d)
        @test ρ.system === sys
        m = dense(ρ)
        @test LA.tr(m) ≈ 1
        @test minimum(LA.eigvals(LA.Hermitian((m + m') / 2))) > -1e-12
        @test_ok fidelity(ρ, State{Pure}(sys, "Up"))
    end
    sq = System(4, Fermion(conserve = N))
    conf = ["Occ", "Emp", "Occ", "Emp"]
    ρq = RandomState{Mixed}(sq, conf, 4)
    @test ρq.system === sq
    @test_ok fidelity(ρq, State{Pure}(sq, conf))
    # the local state of the purification may be given once, as an amplitude vector or by its
    # index, as for State
    @test trace(RandomState{Mixed}(System(3, Qubit()), [1., 0.], 4)) ≈ 1
    @test trace(RandomState{Mixed}(System(3, Qubit()), 1, 4)) ≈ 1

    # randomizing a state of one site leaves the one it was given alone: its tensor was
    # written in place when the element type already matched, and its cache went stale
    for (elt, s) in [(Float64, State{Pure}(System(1, Qubit()), "Up")),
                     (ComplexF64, RandomState{Pure}(System(1, Qubit()), 1))]
        t = copy(s.state[1])
        z = expect(s, Z(1))
        RandomState(elt, s, 2)
        @test norm(s.state[1] - t) < 1e-14
        @test expect(s, Z(1)) ≈ z
    end
end

@testset "A local state given by its index" begin
    # repeated on every site, it is that basis state: fill(1, n) is a Vector{Int}, which was
    # taken for the amplitudes of a single site
    sys = System(3, Qubit())
    @test norm(State{Pure}(sys, 1) - State{Pure}(sys, [0., 1.])) < 1e-12
    @test norm(State{Mixed}(sys, 1) - State{Mixed}(sys, [0., 1.])) < 1e-12
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

@testset "A noise given as an integer" begin
    # a noise is a real number, an integer one included, stored as the Float64 dmrg takes
    @test GroundState(hamiltonian = Z(1), limits = Limits(), nsweeps = 2, noise = 0).noise === 0.
    @test GroundState(hamiltonian = Z(1), limits = Limits(), nsweeps = 2, noise = [1, 0]).noise == [1., 0.]
    @test SteadyState(lindbladian = Dissipator(Sm)(1), limits = Limits(), nsweeps = 2,
                      noise = 0).noise === 0.
    @test SteadyState(lindbladian = Dissipator(Sm)(1), limits = Limits(), nsweeps = 2,
                      noise = [1, 0]).noise == [1., 0.]
end

@testset "The parameters of a Krylov method" begin
    # a tolerance is a real number, an integer one included. A field left to nothing is not
    # passed on, so that the method keeps its default, and `dim` is the `krylovdim` of
    # KrylovKit, which dmrg takes with the prefix of ITensorMPS
    @test Krylov(tol = 1).tol === 1.
    kw = TensorMixedStates.krylov_kwargs
    @test kw(Krylov()) == (;)
    @test kw(Krylov(dim = 8, maxiter = 3, tol = 1e-10)) == (krylovdim = 8, maxiter = 3, tol = 1e-10)
    @test kw(Krylov(dim = 8), "eigsolve_") == (eigsolve_krylovdim = 8,)
end

@testset "The phases that make the state" begin
    # the time_start of SimData is the time of the simulation from its start, which CreateState
    # keeps unless given one of its own
    sim = runTMS(SimData(time_start = 2.5, phases = [CreateState{Pure}(2, Qubit(), "Up")]);
                 output = devnull)
    @test sim.time == 2.5
    sim = runTMS(SimData(phases = [CreateState{Pure}(2, Qubit(), "Up"; time_start = 1.)]);
                 output = devnull)
    @test sim.time == 1.
    # a State cannot be randomised into a mixed state, having no purification to draw from
    st = State{Pure}(System(2, Qubit()), "Up")
    @test_throws "cannot randomize a State into a mixed state" runTMS(SimData(phases = [
        CreateState(type = Mixed(), state = st, randomize = 4)]); output = devnull)
    # ToMixed holds a state already mixed to its limits as well
    ρ = mix(RandomState{Pure}(System(4, Qubit()), 4))
    sim = runTMS(SimData(phases = [CreateState(type = Mixed(), state = ρ),
                                   ToMixed(limits = Limits(maxdim = 2))]); output = devnull)
    @test maxlinkdim(sim.state) ≤ 2
    # a system conserving something strongly is refused a random mixed state by the message
    # that says why, and not sent to another form that refuses it too
    @test_throws "conserving something strongly" RandomState{Mixed}(System(2, Fermion(conserve = strong(N))), 4)
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
    # and it points to the constructor that splits it, given its sites
    @test_throws "create it with the sites it acts on" make_mpo(st, m2(1, 2))
    # a matrix of the wrong size for its site is refused by a message naming the operator and
    # the site, on the three paths. It failed on a DimensionMismatch from reshape
    b3 = Operator{1}("B3", [1. 0. 0. ; 0. 1. 0. ; 0. 0. 1.], plain_op)
    @test_throws "B3 is given by a 3×3 matrix and cannot act on Qubit()" expect(st, b3(1))
    @test_throws "whose dimension is 2" make_mpo(st, b3(1))
    @test_throws "whose dimension is 2" apply(b3(1), st)
    @test_throws "whose dimension is 4" apply(Operator{2}("B9", ones(9, 9), plain_op)(1, 2), st)
    # partial_trace used to skip a position it did not find when tracing, leaving the state
    # whole, and to raise a BoundsError when keeping it
    @test_throws "given site 5, which the state does not have" partial_trace(stm, [5])
    @test_throws "given site 0" partial_trace(stm, [0, 1]; keep = true)
    @test_throws "does not have" renyi2(stm, [2, 7])
    # an operator placed on no site, or a superoperator, is refused by expect by name rather
    # than by a MethodError about iterating it
    @test_throws "such as X(1) rather than X" expect(st, X)
    @test_throws "such as Swap(1,2) rather than Swap" expect(st, Swap)
    @test_throws "is a superoperator" expect(stm, Gate(X))

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

@testset "An operator of several sites on a repeated site" begin
    # its matrix acts on distinct sites and says nothing of a site taken twice: a gate put one
    # index into its tensor twice. The product on one site is written as a product
    for op in (() -> Swap(1, 1), () -> (X ⊗ Y)(2, 2), () -> controlled(X)(3, 3),
               () -> (X ⊗ Id ⊗ Z)(1, 2, 1))
        @test_throws "repeats a site" op()
    end
    @test_ok (X ⊗ Y)(1, 2)
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

@testset "A phase refuses at once what it could not run" begin
    # refused when the phase ran, a wrong field could not be corrected and resumed, the
    # checkpoint then belonging to a simulation of other phases
    @test_throws "of order 5 is not implemented" ApproxW(order = 5)
    @test_throws "w=1 or 2 (not 3)" ApproxW(order = 2, w = 3)
    @test_throws "apply_algo is" ApproxW(order = 2, apply_algo = "fit")
    @test_throws "mpo_algo is" SteadyState(lindbladian = Dissipator(Sm)(1), limits = Limits(),
                                           nsweeps = 2, mpo_algo = "zipp")
    @test_throws "cannot be zero" Evolve(duration = 1., time_step = 0, algo = Tdvp(),
                                         evolver = -im * X(1))
    @test_ok ApproxW(order = 4, w = 1, apply_algo = "naive")
end

@testset "A pure state dropping part of what it conserves" begin
    # a site left conserving nothing is cut into a single block on a system still charged,
    # where the relabelled index keeps one block per basis state: weaken refused to put the
    # one on the other, and so did the Fidelity of a pure reference on a weakened state
    sys = System([Qubit(conserve = N), Boson(3, conserve = named(N, "Nb")),
                  Boson(3, conserve = named(N, "Nb")), Qubit(conserve = N)])
    p(v) = State{Pure}(sys, v)
    ψ = normalize(p(["1", "2", "0", "0"]) + 0.5 * p(["0", "1", "1", "1"]) -
                  0.7im * p(["1", "0", "2", "0"]))
    numbers(s) = [ real(expect(s, N(k))) for k in 1:4 ]
    for target in (N, named(N, "Nb"))
        w = weaken(ψ, target)
        @test symmetries(w.system) == symmetries(weaken(sys, target))
        @test numbers(w) ≈ numbers(ψ)
        ρw = weaken(mix(ψ), target)
        @test fidelity(ρw, State(ρw.system, w)) ≈ 1
        @test last(only(measure(ρw, Fidelity(ψ)))) ≈ 1
    end
end

