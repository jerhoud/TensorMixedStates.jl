# The outputs of the measurements: the text tables, to files and to streams, the JSON Lines
# files, Data and data_to_frame, the destinations of other packages, and the marks a resumed run
# reopens them at, all through Outputs, the object a simulation will keep.
#
# An assertion that needs a later layer is kept as `@test_skip`, marked with that layer. What
# goes through runTMS, the log and the redirection of the output of 1.x is in v1/states_io.jl.

using DataFrames

# a destination of another package, which records its rows, and one that refuses X(1) before
# anything is measured
struct Recorder <: Destination
    name::String
end

mutable struct Recording
    rows::Vector{Any}
    closed::Bool
end

const recordings = Dict{String, Recording}()

function Base.open(d::Recorder)
    r = Recording([], false)
    recordings[d.name] = r
    return r
end

TMS.write_row(r::Recording, time, row) = push!(r.rows, (time, row))

function Base.close(r::Recording)
    r.closed = true
    return nothing
end

struct Picky <: Destination end

struct PickyOpen end

Base.open(::Picky) = PickyOpen()

function TMS.check_list(::PickyOpen, names)
    if "X(1)" in names
        error("Picky takes no X(1)")
    end
    return nothing
end

TMS.write_row(::PickyOpen, _, _) = nothing
Base.close(::PickyOpen) = nothing

@testset "The columns of a text table" begin
    # the time then each measurement (decision 8): an operator of one site on each site, a pair
    # on each pair of sites, A(i)B(j) or c(i,j) when it is named, a complex value in two columns,
    # an integer as it is, a Check in three (layer 3, step 4)
    mktempdir() do dir
        o = TMS.Outputs(; dir)
        ψ = State{Pure}(System(2, Qubit()), "+")
        list = [X, "m" => Z, Z(1) * Z(2), (X, Z), "c" => (Z, Z), Sp(1), MaxLinkdim,
                "chk" => c -> Check(1.0, 1.5)]
        TMS.output!(o, ψ, "data.dat" => list; time = 0.)
        TMS.output!(o, ψ, "data.dat" => list; time = 0.1)
        TMS.close!(o)
        lines = readlines(joinpath(dir, "data.dat"))
        @test length(lines) == 3
        @test split(lines[1], '\t') ==
              ["time", "X(1)", "X(2)", "m(1)", "m(2)", "Z(1)*Z(2)",
               "Re(X(1)Z(1))", "Im(X(1)Z(1))", "Re(X(1)Z(2))", "Im(X(1)Z(2))",
               "Re(X(2)Z(1))", "Im(X(2)Z(1))", "Re(X(2)Z(2))", "Im(X(2)Z(2))",
               "c(1,1)", "c(1,2)", "c(2,1)", "c(2,2)", "Re(Sp(1))", "Im(Sp(1))", "MaxLinkdim",
               "chk.value", "chk.reference", "chk.difference"]
        # a table that reads back as numbers, the header naming the columns, with no #
        cells = split(lines[3], '\t')
        @test length(cells) == 24
        v = parse.(Float64, cells)
        @test v[1] ≈ 0.1
        @test v[2:3] ≈ [1, 1]
        @test v[4:5] ≈ [0, 0] atol = 1e-12
        @test v[19:20] ≈ [0.5, 0] atol = 1e-12
        @test cells[21] == "1"
        @test v[22:23] ≈ [1, 1.5]
        @test cells[24] == "0.5"
    end
    mktempdir() do dir
        o = TMS.Outputs(; dir, format = "%.3f")
        # an operator not defined on a site leaves its cell empty, a complex time takes two
        # columns, the format goes for the time and the values alike, and a string keeps to its
        # cell
        mixed = State{Pure}(System([Qubit(), Fermion()]), ["Up", "Occ"])
        TMS.output!(o, mixed, "m.dat" => ["n" => dag(C) * C, Z(1), "s" => c -> "a\tb\nc"];
                    time = 0.3 + 0.2im)
        TMS.close!(o)
        @test readlines(joinpath(dir, "m.dat")) ==
              ["Re(time)\tIm(time)\tn(1)\tn(2)\tZ(1)\ts", "0.300\t0.200\t\t1.000\t1.000\ta b c"]
    end
    # one format, for the time and the values, writing one number
    @test_throws "the format \"%d %d\" does not write one number" TMS.Outputs(; format = "%d %d")
    @test_throws "does not write one number" TMS.Outputs(; format = "no number")
end

@testset "A file holds one list" begin
    # its columns are fixed by its first line: another list is refused before anything is
    # measured, and values of another shape when they are written (decision 8)
    mktempdir() do dir
        o = TMS.Outputs(; dir)
        ψ = State{Pure}(System(2, Qubit()), "Up")
        count = Ref(0)
        f = "f" => c -> (count[] += 1; 1.0)
        TMS.output!(o, ψ, "a.dat" => [X(1)]; time = 0)
        @test_throws "a.dat holds the list X(1), and another needs another file" TMS.output!(
            o, ψ, ["b.dat" => [f], "a.dat" => [Z(1)]]; time = 1)
        @test count[] == 0
        # the same list, through the same file however it is named
        TMS.output!(o, ψ, "./a.dat" => [X(1)]; time = 2)
        TMS.output!(o, ψ, "v.dat" => ["v" => c -> c.v]; time = 0, v = [1., 2.])
        @test_throws "the values of v do not fit the columns" TMS.output!(
            o, ψ, "v.dat" => ["v" => c -> c.v]; time = 1, v = [1., 2., 3.])
        @test_throws "the values of v do not fit the columns" TMS.output!(
            o, ψ, "v.dat" => ["v" => c -> c.v]; time = 1, v = [1., 2im])
        @test_throws "the time does not fit the columns" TMS.output!(
            o, ψ, "v.dat" => ["v" => c -> c.v]; time = 1im, v = [1., 2.])
        TMS.close!(o)
        @test readlines(joinpath(dir, "a.dat")) == ["time\tX(1)", "0\t0", "2\t0"]
        @test isempty(readlines(joinpath(dir, "b.dat")))
        @test readlines(joinpath(dir, "v.dat")) == ["time\tv[1]\tv[2]", "0\t1\t2"]
    end
    # a file is emptied at its first opening, and a measurement cannot be named time
    mktempdir() do dir
        write(joinpath(dir, "old.dat"), "from before\n")
        o = TMS.Outputs(; dir)
        ψ = State{Pure}(System(1, Qubit()), "Up")
        TMS.output!(o, ψ, "old.dat" => Z(1); time = 0)
        @test_throws "cannot be named time, the column of the times" TMS.output!(
            o, ψ, "t.dat" => ["time" => Z(1)]; time = 0)
        @test_throws "cannot be named time" TMS.output!(o, ψ, Data("t") => ["time" => Z(1)])
        TMS.close!(o)
        @test readlines(joinpath(dir, "old.dat")) == ["time\tZ(1)", "0\t1"]
    end
end

@testset "Streams" begin
    # stdout and stderr take any list, a new header coming when the columns change, and belong
    # to the process, which keeps them open; "" measures all the same, so that a Check runs
    o = TMS.Outputs()
    @test TMS.destination(o, "stdout") == TMS.destination(o, "-") == TMS.Stream(stdout, "%.8g")
    @test TMS.destination(o, "stderr") == TMS.Stream(stderr, "%.8g")
    @test TMS.destination(o, "") == TMS.Stream(devnull, "%.8g")
    io = IOBuffer()
    s = TMS.Stream(io, "%.8g")
    ψ = State{Pure}(System(2, Qubit()), "Up")
    TMS.output!(o, ψ, s => X(1); time = 0)
    TMS.output!(o, ψ, s => Z(1); time = 1)
    TMS.output!(o, ψ, s => Z(1); time = 2)
    TMS.output!(o, ψ, s => ["v" => c -> c.v]; time = 3, v = [1., 2.])
    TMS.output!(o, ψ, s => ["v" => c -> c.v]; time = 4, v = [1.])
    @test String(take!(io)) ==
          "time\tX(1)\n0\t0\ntime\tZ(1)\n1\t1\n2\t1\ntime\tv[1]\tv[2]\n3\t1\t2\ntime\tv[1]\n4\t1\n"
    @test_throws "Check chk failed" TMS.output!(o, ψ, "" => ["chk" => c -> Check(1., 2.; tol = 0.1)])
    TMS.close!(o)
    @test isopen(io)
    @test isopen(stdout)
end

@testset "JSON Lines" begin
    # one object per call, the time then the measurements in their order, nested as measure gives
    # them: a complex number {"re", "im"}, missing null, a float that is not finite a string, a
    # Check an object, and another list in the same file (decision 8)
    mktempdir() do dir
        o = TMS.Outputs(; dir)
        ψ = State{Pure}(System(2, Qubit()), "+")
        list = [Z, Sp(1), (Z, Z), "nan" => c -> NaN, "inf" => c -> -Inf,
                "chk" => c -> Check([1., 2.], [1., 2.5])]
        TMS.output!(o, ψ, ["data.json" => list, "data.jsonl" => Z(1)]; time = 0.5)
        TMS.output!(o, ψ, "data.json" => [X(1)]; time = 1.5im)
        mixed = State{Pure}(System([Qubit(), Fermion()]), ["Up", "Occ"])
        TMS.output!(o, mixed, "missing.json" => ["n" => dag(C) * C])
        TMS.close!(o)
        lines = readlines(joinpath(dir, "data.json"))
        @test length(lines) == 2
        @test startswith(lines[1], "{\"time\":0.5,\"Z\":")
        d = TMS.JSON.parse(lines[1])
        @test d["Z"] ≈ [0, 0] atol = 1e-12
        @test d["Sp(1)"]["re"] ≈ 0.5
        @test d["Sp(1)"]["im"] ≈ 0 atol = 1e-12
        @test d["(Z, Z)"][1][2] ≈ 0 atol = 1e-12
        @test d["(Z, Z)"][2][2] ≈ 1
        @test d["nan"] == "NaN"
        @test d["inf"] == "-Inf"
        @test d["chk"]["reference"] == [1, 2.5]
        @test d["chk"]["difference"] ≈ 0.5
        @test occursin("\"chk\":{\"value\":[1.0,2.0],\"reference\"", lines[1])
        e = TMS.JSON.parse(lines[2])
        @test e["time"]["im"] ≈ 1.5
        @test e["X(1)"] ≈ 1
        @test TMS.JSON.parse(only(readlines(joinpath(dir, "data.jsonl"))))["Z(1)"] ≈ 0 atol = 1e-12
        n = TMS.JSON.parse(only(readlines(joinpath(dir, "missing.json"))))["n"]
        @test isnothing(n[1])
        @test n[2] ≈ 1
    end
end

@testset "Data and data_to_frame" begin
    # rows kept in memory, the time then the measurements; a table of a column per measurement
    # in the order they first come, each value as it is, missing in a row without it
    o = TMS.Outputs()
    ψ = State{Pure}(System(2, Qubit()), "Up")
    TMS.output!(o, ψ, Data("d") => [Z(1), Z(2)]; time = 0.)
    TMS.output!(o, ψ, Data("d") => [X(1), Z(2), X]; time = 0.1)
    rows = TMS.data_rows(o, "d")
    @test first.(rows[1]) == ["time", "Z(1)", "Z(2)"]
    @test last.(rows[1]) ≈ [0, 1, 1]
    df = data_to_frame(rows)
    @test names(df) == ["time", "Z(1)", "Z(2)", "X(1)", "X"]
    @test df.time == [0., 0.1]
    @test df[!, "Z(2)"] ≈ [1, 1]
    @test ismissing(df[2, "Z(1)"])
    @test ismissing(df[1, "X(1)"])
    @test df[2, "X"] ≈ [0, 0] atol = 1e-12
    # two pairs to one Data in one call are two rows, and nothing gathered is a table of times
    TMS.output!(o, ψ, [Data("e") => X(1), Data("e") => Z(1)])
    @test length(TMS.data_rows(o, "e")) == 2
    empty = data_to_frame(Vector{Pair{String, Any}}[])
    @test names(empty) == ["time"]
    @test size(empty) == (0, 1)
    # a complex value stays complex, and the rows stay there once the outputs are closed
    TMS.output!(o, ψ, Data("c") => Sp(1); time = 1)
    TMS.close!(o)
    @test last(only(TMS.data_rows(o, "c"))[2]) isa ComplexF64
end

@testset "A destination of another package" begin
    # a subtype of Destination with open, write_row and close receives the rows with the values
    # in their types; check_list refuses a list before anything is measured or written; without
    # resume_mark it has no mark, and open(d, mark) opens it anew
    mktempdir() do dir
        o = TMS.Outputs(; dir)
        ψ = State{Pure}(System(2, Qubit()), "+")
        TMS.output!(o, ψ, Recorder("a") => [X(1), Sp(1), MaxLinkdim]; time = 0.5)
        t, row = only(recordings["a"].rows)
        @test t == 0.5
        @test first.(row) == ["X(1)", "Sp(1)", "MaxLinkdim"]
        @test last(row[1]) isa Float64 && last(row[1]) ≈ 1
        @test last(row[2]) isa ComplexF64
        @test last(row[3]) === 1
        count = Ref(0)
        @test_throws "Picky takes no X(1)" TMS.output!(
            o, ψ, ["p.dat" => ["f" => c -> (count[] += 1; 1.0)], Picky() => [X(1)]]; time = 1)
        @test count[] == 0
        @test isempty(readlines(joinpath(dir, "p.dat")))
        TMS.output!(o, ψ, Picky() => [Z(1)]; time = 1)
        @test !haskey(TMS.resume_marks(o), Recorder("a"))
        @test open(Recorder("b"), 42) isa Recording
        TMS.close!(o)
        @test recordings["a"].closed
    end
    # anything else is not a destination
    @test_throws "3 is not a destination" TMS.output!(TMS.Outputs(), State{Pure}(System(1, Qubit()), "Up"),
                                                       3 => Z(1))
end

@testset "Resuming the outputs" begin
    # the marks of the destinations that have one, a resumed run reopening them there: what was
    # written after them is cut, the header is not written again, the list of a file and its
    # columns are remembered, and a Data has its rows back. A stream has no mark
    mktempdir() do dir
        o = TMS.Outputs(; dir)
        ψ = State{Pure}(System(1, Qubit()), "Up")
        list = [Z(1), "s" => c -> c.step]
        out(o, k) = TMS.output!(o, ψ, ["r.dat" => list, "r.json" => list, Data("r") => list,
                                       "" => list]; time = k, step = k)
        out(o, 1)
        out(o, 2)
        marks = TMS.resume_marks(o)
        @test Set(keys(marks)) ==
              Set([TMS.destination(o, "r.dat"), TMS.destination(o, "r.json"), Data("r")])
        out(o, 3)
        TMS.close!(o)
        o2 = TMS.reopen!(TMS.Outputs(; dir), marks)
        out(o2, 4)
        @test_throws "holds the list Z(1), s" TMS.output!(o2, ψ, "r.dat" => Z(1); time = 5)
        @test_throws "the values of s do not fit" TMS.output!(
            o2, ψ, "r.dat" => list; time = 5, step = [1, 2])
        TMS.close!(o2)
        @test readlines(joinpath(dir, "r.dat")) == ["time\tZ(1)\ts", "1\t1\t1", "2\t1\t2", "4\t1\t4"]
        @test [ TMS.JSON.parse(l)["time"] for l in readlines(joinpath(dir, "r.json")) ] == [1, 2, 4]
        @test [ last(first(r)) for r in TMS.data_rows(o2, "r") ] == [1, 2, 4]
    end
end
