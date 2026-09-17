# Getting states in and out.
#
# Goes here: setting a state or part of a state, saving it to a file and reading it
# back, and any other serialisation, with a round trip check whenever possible.

@testset "SetState" begin
    sys = System(3, Qubit())
    st = State{Mixed}(sys, ["Dn", "FullyMixed", "+"])  # arbitrary starting local states

    # SetState overwrites the site regardless of its previous (mixed) content
    st2 = apply(SetState("Up")(1), st)
    @test expect1(st2, Z)[1] ≈ 1

    st3 = apply(SetState("Dn")(2), st)
    @test expect1(st3, Z)[2] ≈ -1

    # other sites are left untouched
    @test expect1(st3, Z)[[1, 3]] ≈ expect1(st, Z)[[1, 3]]

    # SetState also accepts a named mixed state ("FullyMixed" resolves to a density matrix)
    st4 = apply(SetState("FullyMixed")(1), st)
    @test expect1(st4, Z)[1] ≈ 0 atol=1e-12

    # SetState is trace-preserving
    @test trace(st4) ≈ 1

    # SetState also accepts an explicit (mixed) density matrix, not just a named state
    p0 = 0.3
    st5 = apply(SetState([p0 0. ; 0. 1 - p0])(1), st)
    @test expect1(st5, Z)[1] ≈ 2p0 - 1
    @test trace(st5) ≈ 1
end

@testset "Saving and loading states" begin
    dir = mktempdir()
    file = joinpath(dir, "states.h5")

    # pure state round trip
    sys = System(4, Qubit())
    stp = State{Pure}(sys, ["Up", "Dn", "+", "i"])
    save_state(file, "pure", stp)
    lp = load_state(file, "pure")
    @test lp isa State{Pure}
    @test length(lp) == 4
    @test expect1(lp, Z) ≈ expect1(stp, Z)
    @test expect1(lp, X) ≈ expect1(stp, X)

    # mixed state round trip, in the same file
    stm = mix(stp)
    save_state(file, "mixed", stm)
    lm = load_state(file, "mixed")
    @test lm isa State{Mixed}
    @test expect1(lm, Z) ≈ expect1(stm, Z)
    @test trace(lm) ≈ 1

    # both states are still readable
    @test expect1(load_state(file, "pure"), Z) ≈ expect1(stp, Z)

    # saving under a name already used replaces the state
    save_state(file, "pure", State{Pure}(System(2, Qubit()), "Dn"))
    lp2 = load_state(file, "pure")
    @test length(lp2) == 2
    @test expect1(lp2, Z) ≈ [-1, -1]

    # sites with parameters are rebuilt identically
    sysp = System([Qubit(), Boson(4), Spin(3/2), Qboson(0.1, 3)])
    stq = State{Pure}(sysp, ["Up", "2", "1/2", "1"])
    save_state(file, "params", stq)
    lq = load_state(file, "params")
    @test lq.system.sites == sysp.sites

    # a loaded state can still be used for further computations
    @test trace(mix(lq)) ≈ 1

    # SaveState and LoadState phases, the check fails if the state is not restored
    file2 = joinpath(dir, "phases.h5")
    @test_ok test_phases([
        CreateState{Pure}(3, Qubit(), ["Up", "Dn", "Up"]),
        SaveState(file = file2, statename = "st"),
        CreateState{Pure}(1, Qubit(), "Dn"),
        LoadState(file = file2, statename = "st",
            final_measures = check(Z, [1, -1, 1])),
    ])
end
