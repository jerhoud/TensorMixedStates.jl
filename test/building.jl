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
