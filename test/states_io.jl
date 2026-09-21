# Getting data in and out.
#
# Goes here: setting a state or part of a state, saving it to a file and reading it
# back, any other serialisation, and the output side of a simulation, with a round trip
# check whenever possible.

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

@testset "DataFrames extension" begin
    # a Data target gathers its measurements in the data field of the simulation,
    # one entry per name, and it works even when the output is redirected
    sim = runTMS(SimData(name = "datatoframe", phases = [
        CreateState{Pure}(2, Qubit(), ["Z+", "Z-"]),
        Evolve(algo = Tdvp(), duration = 0.2, time_step = 0.1, evolver = -im * X(1),
               measures = Data("obs") => [Norm, Trace]),
    ]); output = devnull)
    @test collect(keys(sim.data)) == ["obs"]
    @test sort(collect(keys(sim.data["obs"]))) == ["Norm", "Trace"]

    # DataToFrame lives in the DataFrames extension, so this also checks that the
    # extension loads at all. Several measures are joined on the time column.
    df = DataToFrame(sim.data["obs"])
    @test df isa DataFrame
    @test sort(names(df)) == ["Norm", "Trace", "time"]      # the order is not stable
    @test size(df) == (2, 3)
    @test df.time ≈ [0.1, 0.2]
    @test df.Norm ≈ [1, 1]
    @test df.Trace ≈ [1, 1]

    # a single measure needs no join and comes back as it is
    sim1 = runTMS(SimData(name = "datatoframe", phases = [
        CreateState{Pure}(2, Qubit(), ["Z+", "Z-"]),
        Evolve(algo = Tdvp(), duration = 0.2, time_step = 0.1, evolver = -im * X(1),
               measures = Data("one") => Norm),
    ]); output = devnull)
    df1 = DataToFrame(sim1.data["one"])
    @test df1 isa DataFrame
    @test sort(names(df1)) == ["Norm", "time"]
    @test size(df1) == (2, 2)
end

@testset "Output of a complex simulation time" begin
    # a complex time reaches the output through the function level interface, and the
    # checkpoint stores its two parts, so it has to be writable. `Printf` refuses it, and
    # the columns a complex measurement takes are the form to follow
    st = State{Pure}(System(2, Qubit()), "Up")
    line(t) = (io = IOBuffer(); output(Simulation(st; time = t), io, "h", [1.0]);
               String(take!(io)))
    @test line(0.3 + 0.2im) == "h\t     0.3\t     0.2\t             1\n"
    # and a real time keeps exactly the format it had. The default is a float, but an
    # integer given explicitly must still be written with the time format, unlike a
    # measured value, which keeps its own form when it is not a float
    @test line(0.) == "h\t       0\t             1\n"
    @test line(0) == "h\t       0\t             1\n"
    @test line(0.25) == "h\t    0.25\t             1\n"
end

@testset "CreateState with a State object" begin
    # `type` is what the phase was asked for, so a State handed to it in the other
    # representation must be converted and not silently kept
    sys = System(3, Qubit())
    pure = State{Pure}(sys, "Up")
    run1(phase) = runTMS(SimData(; phases = [phase]); output = devnull).state
    @test run1(CreateState(type = Pure(), state = pure)) isa State{Pure}
    @test run1(CreateState(type = Mixed(), state = pure)) isa State{Mixed}
    @test_throws "cannot be turned back into a pure state" run1(
        CreateState(type = Pure(), state = mix(pure)))
    # and randomizing needs a system to randomize over
    @test_throws "needs a system to create a random state" run1(
        CreateState(type = Pure(), randomize = 10))
end

@testset "Standard streams are not closed" begin
    # "stdout", "stderr" and "" are output destinations like any other, but the streams
    # they name belong to the process and a simulation must leave them open on its way
    # out. The run is redirected so that a regression cannot take the output of the test
    # run down with it.
    mktempdir() do dir
        cd(dir) do
            phases = [
                CreateState{Pure}(2, Qubit(), "Up"),
                Gates(gates = X(1),
                      final_measures = ["stdout" => Z, "stderr" => Norm, "" => Linkdim]),
            ]
            still_open = open("captured", "w") do io
                redirect_stdout(io) do
                    redirect_stderr(io) do
                        runTMS(SimData(; name = "streams", phases))
                        (isopen(stdout), isopen(stderr))
                    end
                end
            end
            @test still_open == (true, true)
        end
    end
end

@testset "Output survives a failing phase" begin
    # a phase that fails must not take away what was collected before it: a json
    # destination is only written when the simulation closes its files, so that has to
    # happen on the way out of an exception too
    mktempdir() do dir
        cd(dir) do
            boom = StateFunc("boom", _ -> error("phase failure on purpose"))
            phases = [
                CreateState{Pure}(2, Qubit(), "Up"),
                Gates(gates = X(1), final_measures = ["data.json" => Norm]),
                Gates(gates = X(1), final_measures = ["data.json" => boom]),
            ]
            @test_throws ErrorException runTMS(SimData(; name = "failing", phases))
            @test isfile(joinpath("failing", "data.json"))
            @test occursin("Norm", read(joinpath("failing", "data.json"), String))
        end
    end
end

@testset "Loading a state onto a system" begin
    dir = mktempdir()
    file = joinpath(dir, "onto.h5")
    sys = System(3, Qubit())
    a = RandomState{Pure}(sys, 4)
    save_state(file, "a", a)

    # read back as it comes, the state lands on a system built from the file and cannot be
    # compared with the one it was saved from
    b = load_state(file, "a")
    @test b.system !== sys
    @test_throws "do not share their System" inner(a, b)

    # read onto an existing system, it lands where it can be compared and is the same state
    c = load_state(file, "a"; system = sys)
    @test c.system === sys
    @test fidelity(a, c) ≈ 1
    @test inner(a, c) ≈ norm(a)^2
    @test expect1(c, Z) ≈ expect1(a, Z)

    # the sites of the system given must be the ones of the state
    @test_throws "system of other sites" load_state(file, "a"; system = System(2, Qubit()))

    # and the same for a mixed state
    m = mix(a)
    save_state(file, "m", m)
    lm = load_state(file, "m"; system = sys)
    @test lm.system === sys
    @test hs_fidelity(m, lm) ≈ 1
end
