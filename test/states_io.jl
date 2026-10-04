# Getting data in and out.
#
# Goes here: setting a state or part of a state, saving it to a file and reading it
# back, any other serialisation, and the output side of a simulation, with a round trip
# check whenever possible.

# Sites whose fields are not all numbers. Version 1 of the state file format wrote every
# field as a `Float64` and could carry none of these; a site declaring which quantum numbers
# it conserves needs exactly that, so the round trip is checked below. The fields cover every
# kind a state file accepts.
struct Kindly <: AbstractSite
    conserve::Symbol
    on::Bool
    label::String
    n::Int
    x::Float64
    void::Nothing
end

TensorMixedStates.dim(::Kindly) = 2

# a field of a kind no state file can carry, to check that the refusal names it
struct Unkindly <: AbstractSite
    range::UnitRange{Int}
end

# a site type of a module inside another one, which a state file named by the last name of its
# module alone and could not find again
module Enclosing
    using TensorMixedStates
    struct Enclosed <: AbstractSite end
    TensorMixedStates.dim(::Enclosed) = 2
end

TensorMixedStates.dim(::Unkindly) = 2

# fields a state file wrote as their printed form and could not read back: a float of another
# width prints as "0.1f0", and an integer beyond `Int` overflowed
struct Widely <: AbstractSite
    y::Float32
    big::UInt64
end

TensorMixedStates.dim(::Widely) = 2

# a site whose inner constructor replaces the one taking its fields, by which a state file and
# weaken rebuild a site
struct Guarded <: AbstractSite
    conserve::String
    Guarded(; conserve = ()) = new(conserve_string(new(""), conserve))
end

TensorMixedStates.dim(::Guarded) = 2

@def_operators(Guarded(), [
    selfadjoint_op => [
        N = [0. 0. ; 0. 1.],
    ],
])

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

@testset "Weakening beside a site whose conserve field is not a string" begin
    # such a site conserves nothing, and weakening, which every measurement of a state
    # conserving something strongly goes through, leaves it as it is
    sys = System([Kindly(:sz, true, "a label", 3, 1.5, nothing), Fermion(conserve = strong(N))])
    ρ = mix(State{Pure}(sys, Any[[1., 0.], "Occ"]))
    @test trace(ρ) ≈ 1
    @test expect(ρ, N(2)) ≈ 1
    @test_ok weaken(ρ, ())
end

@testset "A site rebuilt from its fields must keep the constructor taking them" begin
    sys = System([Guarded(conserve = strong(N)), Guarded(conserve = strong(N))])
    st = State{Pure}(sys, Any[[0., 1.], [1., 0.]])
    @test expect(st, N(1)) ≈ 1
    @test_throws "Guarded must keep its constructor Guarded(conserve) for its conserved" weaken(st, ())
    mktempdir() do dir
        file = joinpath(dir, "guarded.h5")
        save_state(file, "g", st)
        @test_throws "Guarded must keep its constructor Guarded(conserve) for a state file" load_state(file, "g")
    end
end

@testset "Measurement sets and their rows" begin
    mktempdir() do dir
        cd(dir) do
            # a Measure among the measurements stands for its measurements, and a string is a
            # measurement in final_measures too, a line with its label and the time
            sim = runTMS(SimData(name = "sets", phases = [
                CreateState{Pure}(2, Qubit(), "Up"),
                Gates(gates = X(1), final_measures = [Data("d") => [Measure(X(1), Z(1)), Z(2)],
                                                      "out.dat" => "label", Data("s") => "label"])]))
            d = sim.data["d"]
            @test sort(collect(keys(d))) == ["X(1)", "Z(1)", "Z(2)"]
            @test only(d["Z(1)"]["data"]) ≈ -1
            @test only(d["Z(2)"]["data"]) ≈ 1
            @test startswith(readline("sets/out.dat"), "label\t")
            @test haskey(sim.data["s"], "label")

            # a frame has a row per measurement set: the time repeats over a circuit or once it
            # is set back, and joining on it paired Z(1), measured in the first phase, with
            # X(1), measured in the second
            sim = runTMS(SimData(name = "rows", phases = [
                CreateState{Pure}(2, Qubit(), "Up"),
                Gates(gates = X(1), final_measures = Data("d") => [Z(1), Z(2)]),
                Gates(gates = X(2), final_measures = Data("d") => [X(1), Z(2)])]))
            d = sim.data["d"]
            @test d["Z(2)"]["events"] == [1, 2]
            df = data_to_frame(d)
            @test size(df, 1) == 2
            @test all(df[!, "Z(2)"] .≈ [1, -1])
            @test ismissing(df[2, "Z(1)"])
            @test ismissing(df[1, "X(1)"])
        end
    end
end

@testset "SetState under a strong symmetry" begin
    # resetting a site moves the charge of one side of the density matrix only, which a strong
    # symmetry forbids: it kept the block of charge zero alone, a state of trace zero
    ss = System(2, Fermion(conserve = strong(N)))
    ρs = mix(State{Pure}(ss, ["Occ", "Emp"]))
    @test_throws "strongly forbids" apply(SetState("Emp")(1), ρs)
    @test_throws "strongly forbids" make_mpo(ρs, SetState("Emp")(1))
    # conserved weakly, it resets the site as it does without charges
    ρw = mix(State{Pure}(System(2, Fermion(conserve = N)), ["Occ", "Emp"]))
    r = apply(SetState("Emp")(1), ρw)
    @test trace(r) ≈ 1
    @test expect(r, N(1)) ≈ 0 atol = 1e-12
end

@testset "Loading a partial trace of a charged system" begin
    # it keeps charged indices on the sites that conserve nothing, and loading has to rebuild
    # the pure ones in that mode
    file = joinpath(mktempdir(), "traced.h5")
    sq = System([Qubit(), Fermion(conserve = N)])
    ρ = partial_trace(mix(State{Pure}(sq, ["+", "Occ"])), [1]; keepers = true)
    save_state(file, "traced", ρ)
    lt = load_state(file, "traced")
    @test trace(lt) ≈ 1
    @test expect(lt, X(1)) ≈ expect(ρ, X(1))
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

    # a site type of a module inside another one is rebuilt as well, its module being found
    # from the root one
    nested = State{Pure}(System(2, Enclosing.Enclosed()), ["0", "1"])
    save_state(file, "nested", nested)
    @test load_state(file, "nested").system.sites == nested.system.sites

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

@testset "A state file with its conserved quantities unsorted" begin
    # written before they were sorted, in the order of their declaration: the state read back
    # sat on sites no longer equal to those built now, and could not be put on their system
    old = Electron("Ntot:0,1,1,2;2Sz:0,1,-1,0")
    file = joinpath(mktempdir(), "old.h5")
    save_state(file, "st", State{Pure}(System(2, old), "Up"))
    sys = System(2, Electron(conserve = (Ntot, 2Sz)))
    st = load_state(file, "st")
    @test st.system.sites == sys.sites
    @test_ok load_state(file, "st"; system = sys)
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

    # data_to_frame lives in the DataFrames extension, so this also checks that the
    # extension loads at all. Several measures are joined on the time column.
    df = data_to_frame(sim.data["obs"])
    @test df isa DataFrame
    @test names(df) == ["time", "Norm", "Trace"]            # in the order of their names
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
    df1 = data_to_frame(sim1.data["one"])
    @test df1 isa DataFrame
    @test sort(names(df1)) == ["Norm", "time"]
    @test size(df1) == (2, 2)

    # a measurement named time or event keeps its column, renamed, beside the time of its rows
    rows(v) = Dict("events" => [1, 2], "times" => [0.1, 0.2], "data" => v)
    dfc = data_to_frame(Dict("time" => rows([5., 6.]), "event" => rows([7., 8.]), "Norm" => rows([1., 1.])))
    @test sort(names(dfc)) == ["Norm", "event_1", "time", "time_1"]
    @test dfc.time ≈ [0.1, 0.2]
    @test dfc.time_1 ≈ [5, 6]
    @test dfc.event_1 ≈ [7, 8]
    @test sort(names(data_to_frame(Dict("time" => rows([5., 6.]))))) == ["time", "time_1"]

    # one call of output is one row, two pairs going to one Data made two, and nothing
    # gathered is an empty table, where outerjoin of no table failed
    simd = Simulation(State{Pure}(System(2, Qubit()), "Up"); output = devnull)
    output(simd, [Data("d") => X(1), Data("d") => Z(1)])
    @test size(data_to_frame(simd.data["d"])) == (1, 3)
    @test size(data_to_frame(Dict())) == (0, 1)
end

@testset "Output of a complex simulation time" begin
    # a complex time reaches the output through the function level interface, and the
    # checkpoint stores its two parts, so it has to be writable. `Printf` refuses it, and
    # the columns a complex measurement takes are the form to follow
    st = State{Pure}(System(2, Qubit()), "Up")
    line(t) = (io = IOBuffer(); output(Simulation(st; time = t, output = io), "data" => Trace);
               String(take!(io)))
    @test line(0.3 + 0.2im) == "Trace\t     0.3\t     0.2\t             1\n"
    # and a real time keeps exactly the format it had. The default is a float, but an
    # integer given explicitly must still be written with the time format, unlike a
    # measured value, which keeps its own form when it is not a float
    @test line(0.) == "Trace\t       0\t             1\n"
    @test line(0) == "Trace\t       0\t             1\n"
    @test line(0.25) == "Trace\t    0.25\t             1\n"
end

@testset "Output of complex values" begin
    # a complex value takes two columns, its real part then its imaginary part, an imaginary
    # one a single column under the name Im(...), and a value holding several is written
    # number by number
    lines(st, m) = (io = IOBuffer(); output(Simulation(st; output = io), "data" => m);
                    split(String(take!(io)), '\n'; keepempty = false))
    fields(line) = split(line, '\t')
    ph = apply(Phase(0.7)(1), State{Pure}(System(2, Qubit()), "+"))
    l = only(lines(ph, Sp(1)))
    @test fields(l)[1] == "Sp(1)"
    @test parse.(Float64, fields(l)[3:end]) ≈ [cos(0.7), sin(0.7)] / 2 atol = 1e-7
    ls = lines(ph, (Sp, Sm))
    @test ls[1] == "SpSm"
    @test length.(fields.(ls[2:3])) == [6, 6]
    l = only(lines(State{Pure}(System(2, Qubit()), ["+", "i"]), im * X(1) * Y(2)))
    @test fields(l)[1] == "Im(im*X(1)*Y(2))"
    @test length(fields(l)) == 3
    # the header, the time, the two complex values of Sp, the reference and the distance
    l = only(lines(State{Pure}(System(2, Qubit()), "+"), Check("c", Sp, [0.5, 0.5])))
    @test length(fields(l)) == 2 + 4 + 2 + 1
    @test !occursin('[', l)
    # a real reference and the distance keep one column each next to a complex value
    l = only(lines(State{Pure}(System(2, Qubit()), "+"), Check("c", Sp(1), 0.5)))
    @test length(fields(l)) == 2 + 2 + 1 + 1
    @test first.(fields.(lines(ph, [[X(1), Z(2)]]))) == ["X(1)", "Z(2)"]

    # a warning of `measure` goes to the log of the simulation, the same stream here, and
    # none reaches the logger of the caller, which still gets the warnings of anything else
    ls = @test_logs lines(ph, RealValue(Sp(1)))
    @test startswith(ls[1], "WARNING: large imaginary part: time 0.0, Sp(1)")
    @test_logs (:warn, "from the user") lines(ph, StateFunc("user", _ -> (@warn "from the user"; 1.0)))

    # a Data destination holds a complex value as it is, and a json file writes it as
    # {"re": …, "im": …} whatever JSON.jl would make of it, and a matrix as the list of its
    # rows, as the text file writes it. A complex simulation time is written the same way
    mktempdir() do dir
        cd(dir) do
            ms = [X(1), Sp(1), (Sp, Sm)]
            sim = runTMS(SimData(name = "cplx", phases = [
                CreateState{Pure}(2, Qubit(), "+"),
                Gates(gates = Phase(0.7)(1), final_measures = [Data("d") => ms, "out.json" => ms])]))
            d = sim.data["d"]
            @test only(d["X(1)"]["data"]) isa Float64
            @test only(d["Sp(1)"]["data"]) ≈ exp(0.7im) / 2
            m = only(d["SpSm"]["data"])
            @test m isa Matrix{ComplexF64}
            @test nonmissingtype(eltype(data_to_frame(d)[!, "Sp(1)"])) == ComplexF64
            js = TensorMixedStates.JSON.parsefile("cplx/out.json")
            sp = only(js["Sp(1)"]["data"])
            @test complex(sp["re"], sp["im"]) ≈ exp(0.7im) / 2
            e = only(js["SpSm"]["data"])[1][2]
            @test complex(e["re"], e["im"]) ≈ m[1, 2]

            runTMS(SimData(name = "ctime", phases = [
                CreateState{Pure}(2, Qubit(), "Up"),
                Gates(gates = X(1), time_start = 0.5im, final_measures = "out.json" => Z(1))]))
            t = only(TensorMixedStates.JSON.parsefile("ctime/out.json")["Z(1)"]["times"])
            @test complex(t["re"], t["im"]) == 0.5im
        end
    end
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
                      final_measures = ["stdout" => Z, "stderr" => Norm, "" => MaxLinkdim]),
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

@testset "Site fields of every kind" begin
    dir = mktempdir()
    file = joinpath(dir, "fields.h5")

    site = Kindly(:sz, true, "a label", 3, 1.5, nothing)
    st = State{Pure}(System(2, site), [1., 0.])
    save_state(file, "kinds", st)
    l = load_state(file, "kinds")

    # each field comes back as what it was, not as the Float64 version 1 would have made of it
    back = l.system.sites[1]
    @test back.conserve === :sz
    @test back.on === true
    @test back.label == "a label"
    @test back.n === 3
    @test back.x === 1.5
    @test back.void === nothing
    @test l.system.sites == st.system.sites

    # the first field is named conserve without being one, being a Symbol and not a string.
    # A site carrying such a field conserves nothing rather than making the machinery choke
    @test TensorMixedStates.conserved(site) == ""

    # a field of an unsupported kind is refused by a message naming the site and the field,
    # rather than by a MethodError raised by convert somewhere inside HDF5
    bad = State{Pure}(System(2, Unkindly(1:3)), [1., 0.])
    @test_throws "its field range is a" save_state(file, "bad", bad)

    # and the refusal comes before the file is touched, so saving a state that cannot be
    # written over a name already in use leaves what was there intact
    @test_throws "its field range is a" save_state(file, "kinds", bad)
    @test load_state(file, "kinds").system.sites == st.system.sites

    wide = Widely(0.1f0, typemax(UInt64))
    save_state(file, "wide", State{Pure}(System(2, wide), [1., 0.]))
    @test load_state(file, "wide").system.sites[1] === wide
end

@testset "Reading a version 1 state file" begin
    # written by the released 1.3.0, see reference/make_state_v1.jl. Version 1 wrote every
    # site field as a Float64, and `load_state` has to go on reading it: a checkpoint left by
    # such a version holds a state in that format, so a resume depends on it
    file = joinpath(@__DIR__, "reference", "state_v1.h5")

    p = load_state(file, "pure")
    @test p isa State{Pure}
    @test length(p) == 4
    @test p.system.sites == [Qubit(), Boson(4), Spin(3/2), Qboson(0.1, 3)]
    @test expect(p, Z(1)) ≈ 1
    @test expect(p, N(2)) ≈ 2
    @test expect(p, Sz(3)) ≈ 0.5

    m = load_state(file, "mixed")
    @test m isa State{Mixed}
    @test m.system.sites == p.system.sites
    @test trace(m) ≈ 1
    @test expect(m, Z(1)) ≈ 1

    # a file of an unknown version is still refused, by a message naming what is accepted.
    # The group needs nothing but its version attribute, the check coming first
    future = joinpath(mktempdir(), "future.h5")
    TensorMixedStates.HDF5.h5open(future, "cw") do h
        g = TensorMixedStates.HDF5.create_group(h, "s")
        TensorMixedStates.HDF5.attributes(g)["version"] = 99
    end
    @test_throws "expected one of 1, 2" load_state(future, "s")
end

@testset "Reading a version 2 state file" begin
    # written by the released 1.6.0, see reference/make_state_v2.jl: a file save_state writes
    # stays readable by every later version, the current format included
    file = joinpath(@__DIR__, "reference", "state_v2.h5")
    strong = TensorMixedStates.strong

    e = load_state(file, "electrons")
    @test e isa State{Pure}
    @test e.system.sites == fill(Electron(conserve = (Ntot, 2Sz)), 2)
    @test real(expect1(e, Nup)) ≈ [1, 0]
    @test real(expect1(e, Ndn)) ≈ [0, 1]

    ρ = load_state(file, "fermions_strong")
    @test ρ isa State{Mixed}
    @test ρ.system.sites == fill(Fermion(conserve = strong(N)), 2)
    @test real(expect1(ρ, N)) ≈ [0.8, 0.2]
    @test expect(ρ, dag(C)(1) * C(2)) ≈ 0.4
    @test trace(ρ) ≈ 1

    k = load_state(file, "kindly")
    @test k.system.sites[1] == Kindly(:none, true, "a", 3, 0.5, nothing)
    @test expect(k, Z(2)) ≈ 1
end

@testset "Saving a state that carries charges" begin
    strong = TensorMixedStates.strong
    dir = mktempdir()
    file = joinpath(dir, "charged.h5")

    # the quantities a site conserves travel in the string it keeps, the strength included,
    # so a state read back draws the same indices and measures the same numbers
    for (name, site) in (("weak", Fermion(conserve = N)), ("strong", Fermion(conserve = strong(N))))
        sys = System(3, site)
        stp = State{Pure}(sys, ["Occ", "Emp", "Occ"])
        save_state(file, name, stp)
        lp = load_state(file, name)
        @test lp isa State{Pure}
        @test lp.system.sites == sys.sites
        @test flux(lp.state) == flux(stp.state)
        @test expect1(lp, N) ≈ expect1(stp, N)

        stm = mix(stp)
        save_state(file, name * "-mixed", stm)
        lm = load_state(file, name * "-mixed")
        @test lm isa State{Mixed}
        @test flux(lm.state) == flux(stm.state)
        @test trace(lm) ≈ 1
        @test expect1(lm, N) ≈ expect1(stm, N)
    end
end

@testset "Printing Limits" begin
    # printed as the call that builds it, the fields left at their default omitted, so that
    # the text evaluates back to the same fields
    @test repr(Limits()) == "Limits()"
    @test repr(Limits(cutoff = 1e-10, maxdim = 50)) == "Limits(cutoff = 1.0e-10, maxdim = 50)"
    for l in [Limits(), Limits(cutoff = 1e-10, maxdim = 50), Limits(maxdim = [10, 20, 50]),
              Limits(cutoff = 1e-14, maxdim = 100, mindim = 10)]
        back = eval(Meta.parse(repr(l)))
        @test all(getfield(back, f) == getfield(l, f) for f in fieldnames(Limits))
    end
end

@testset "Printing Krylov" begin
    # as Limits is, the fields left to the default of the method omitted
    @test repr(Krylov()) == "Krylov()"
    @test repr(Krylov(dim = 8, tol = 1e-10)) == "Krylov(dim = 8, tol = 1.0e-10)"
    for k in [Krylov(), Krylov(dim = 8, tol = 1e-10), Krylov(maxiter = 3)]
        back = eval(Meta.parse(repr(k)))
        @test all(getfield(back, f) == getfield(k, f) for f in fieldnames(Krylov))
    end
end

@testset "The stamp records the threading" begin
    # the running time depends on these settings of the process, so a run records them
    mktempdir() do dir
        cd(dir) do
            runTMS(SimData(name = "stamped", phases = [CreateState{Pure}(2, Qubit(), "Up")]))
            stamp = readlines(joinpath("stamped", "stamp"))
            @test any(l -> startswith(l, "BLAS lib") && endswith(l, " threads"), stamp)
            @test "Julia threads $(Threads.nthreads())" in stamp
            @test any(l -> startswith(l, "Strided threads "), stamp)
            @test "Block sparse multithreading off" in stamp
            @test "CPU threads $(Sys.CPU_THREADS)" in stamp
        end
    end
end
