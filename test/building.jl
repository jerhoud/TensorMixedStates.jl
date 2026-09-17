# Building the basic objects, and package wide quality checks.
#
# Goes here: anything that only checks that an object can be constructed and has the
# announced shape, for System, State, Simulation and RandomState, plus Aqua. What a
# state *measures* belongs to sites.jl, not here.

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
    @test maxlinkdim(RandomState{Mixed}(sys, 8)) ≤ 8
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

@testset "Simulation building" begin
   @test_ok Simulation(nothing)
   @test_pm Simulation(State{type}(System(3, Qubit()), "Up"))
end
