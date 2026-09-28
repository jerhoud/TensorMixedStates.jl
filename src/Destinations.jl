export Data

"""
    Data(name)

represent a storage with the given name where to put measurement data
"""
struct Data
    name::String
end

############### the values of an accumulating destination ###############

"""
    json_value(x)

a value of a json destination as it is written: a complex number becomes
`{"re": …, "im": …}`, wherever it is, which JSON.jl writes that way or refuses depending on
its version, a matrix the list of its rows, as a file writes it, and a value that is not
finite the string Julia prints it as, `"Inf"`, `"-Inf"` or `"NaN"`, json having no number for
it.
"""
json_value(x::AbstractFloat) = isfinite(x) ? x : string(x)
json_value(x::Complex) = Dict("re" => json_value(real(x)), "im" => json_value(imag(x)))
json_value(x::AbstractMatrix) = [ json_value.(x[i, :]) for i in axes(x, 1) ]
json_value(x::AbstractArray) = map(json_value, x)
json_value(x::AbstractDict) = Dict(k => json_value(v) for (k, v) in x)
json_value(x) = x

"""
    checkpoint_value(x)
    restored_value(x)

a value written into a checkpoint, and read back from it. Json holds neither complex numbers
nor matrices, a matrix coming back as the vector of its columns, nor a number that is not
finite, so all three are marked and rebuilt, and an array read back is given its element type
again: a resumed run hands back the values an uninterrupted one would.
"""
checkpoint_value(x::AbstractFloat) = isfinite(x) ? x : Dict("float" => string(x))
checkpoint_value(x::Complex) = Dict("complex" => [checkpoint_value(real(x)), checkpoint_value(imag(x))])
checkpoint_value(x::AbstractMatrix) = Dict("matrix" => [ checkpoint_value(x[i, :]) for i in axes(x, 1) ])
checkpoint_value(x::AbstractArray) = map(checkpoint_value, x)
checkpoint_value(x::AbstractDict) = Dict(k => checkpoint_value(v) for (k, v) in x)
checkpoint_value(x) = x

function restored_value(x::AbstractDict)
    if haskey(x, "complex")
        r, i = restored_value.(x["complex"])
        return complex(r, i)
    elseif haskey(x, "float")
        return parse(Float64, x["float"])
    elseif haskey(x, "matrix")
        return stack(restored_value.(x["matrix"]); dims = 1)
    end
    return Dict(k => restored_value(v) for (k, v) in x)
end
restored_value(x::AbstractVector) = map(restored_value, x)
restored_value(x) = x

"""
    Series

what an accumulating destination holds, a json file or a `Data` one: for each header, the
times, the values and the calls of `output` they came from, as
`Dict("times" => [...], "data" => [...], "events" => [...])`. It is the form `sim.data`
hands a `Data` destination over in, and the one `data_to_frame` reads.

A series only ever grows, which is what lets a checkpoint record it by its lengths and write
back, later, exactly the part it held then.
"""
const Series = Dict{String, Dict}

# a series numbers the calls of `output` it takes part in, and each value records the one it
# came from, so that what was measured together can be told apart from what only shares its
# time: the time repeats over the sweeps of dmrg, a circuit, or once it is set back
next_event(s::Series) = 1 + maximum((last(d["events"]) for d in values(s)); init = 0)

function push_value!(s::Series, header, time, value, event::Int)
    d = get!(s, header, Dict("times" => [], "data" => [], "events" => Int[]))
    push!(d["times"], time)
    push!(d["data"], value)
    push!(d["events"], event)
    return nothing
end

series_lengths(s::Series) = Dict{String, Int}(h => length(d["times"]) for (h, d) in s)

# the part of a series its lengths cover, encoded for a checkpoint: a header opened since is
# left out, as the uninterrupted run had not opened it at that point either
persist_series(s::Series, mark::Dict{String, Int}) =
    Dict(h => Dict(k => checkpoint_value(s[h][k][1:n]) for k in ("times", "data", "events"))
         for (h, n) in mark)

# the series are pushed onto by the measurements to come, so they stay vectors of any element
# type, whatever the values read back
restore_series(d) =
    Series(h => Dict("times" => Any[ restored_value(x) for x in s["times"] ],
                     "data" => Any[ restored_value(x) for x in s["data"] ],
                     "events" => Int[ x for x in s["events"] ]) for (h, s) in d)

############### the rows of a text destination ###############

function output_one(file, x::AbstractFloat, format)
    Printf.format(file, format, x)
end

function output_one(file, x::Complex, format)
    output_one(file, real(x), format)
    print(file, "\t")
    output_one(file, imag(x), format)
end

# a value holding several, the part of a `Check` made on a vector observable for instance,
# is written number by number, a matrix row by row as the lines of a matrix are, rather than
# as the literal Julia prints
function output_one(file, x::AbstractArray, format)
    for (k, y) in enumerate(row_major(x))
        if k > 1
            print(file, "\t")
        end
        output_one(file, y, format)
    end
end

function output_one(file, x, _)
    print(file, x)
end

# the time column always takes the time format, unlike a measured value, which keeps its
# own printed form when it is not a float: `Linkdim` is meant to read as 8, not as 8.000
function output_time(file, t::Number, format)
    Printf.format(file, format, t)
end

# `Printf` refuses a complex number outright, so a complex simulation time is written as
# the two columns a complex measurement takes, real then imaginary
function output_time(file, t::Complex, format)
    output_time(file, real(t), format)
    print(file, "\t")
    output_time(file, imag(t), format)
end

write_row(file::IO, formats, time, header, data) =
    write_row(file, formats, time, header, [data])

function write_row(file::IO, formats, time, header, data::Vector)
    print(file, header, "\t")
    output_time(file, time, first(formats))
    for x in data
        print(file, "\t")
        output_one(file, x, last(formats))
    end
    println(file)
end

function write_row(file::IO, formats, time, header, data::Matrix)
    println(file, header)
    for l in 1:size(data, 1)
        write_row(file, formats, time, "$header:$l", data[l,:])
    end
end

############### destinations ###############

"""
    Destination

where the measurements of an `output` go. There are four kinds, and each says for itself what
it takes, how far it has got and what a checkpoint has to carry of it:

- `TextFile`: a file of the simulation, written line by line as the measurements are made;
- `Stream`: a stream of the process, `stdout`, `stderr` or `devnull`, or the one `runTMS` was
  asked to redirect everything to; the process owns it, so it is neither closed nor resumed;
- `JsonFile`: a `.json` file, accumulated in memory and written when the files are closed;
- `DataStore`: a `Data(name)` destination, accumulated in memory and handed over in
  `sim.data`.

A destination answers `emit!` (a call of `output`), `emit_line!` (a line of the log),
`reached` (how far it has got), `persist` (what a checkpoint carries of it up to there) and
`close!`.
"""
abstract type Destination end

struct TextFile <: Destination
    io::IO
end

struct Stream <: Destination
    io::IO
end

struct JsonFile <: Destination
    path::String
    series::Series
end

struct DataStore <: Destination
    series::Series
end

function emit!(d::Union{TextFile, Stream}, formats, time, values)
    for (header, value) in values
        write_row(d.io, formats, time, header, value)
    end
    # written out as they are made, so that a run can be followed, and a crash loses nothing
    # already measured
    flush(d.io)
end

function emit!(d::Union{JsonFile, DataStore}, _, time, values)
    event = next_event(d.series)
    for (header, value) in values
        push_value!(d.series, header, time, value, event)
    end
end

function emit_line!(d::Union{TextFile, Stream}, text)
    println(d.io, text)
    flush(d.io)
end

reached(d::TextFile) = (flush(d.io); position(d.io))
reached(d::Union{JsonFile, DataStore}) = series_lengths(d.series)

# a text file is on the disk up to where it had got already, and a resume cuts it back there;
# what a series held has to travel in the checkpoint, since it is only written at the end
persist(::TextFile, pos::Int) = Dict("text" => pos)
persist(d::JsonFile, lengths) = Dict("json" => persist_series(d.series, lengths))

close!(d::TextFile) = close(d.io)
close!(::Union{Stream, DataStore}) = nothing

function close!(d::JsonFile)
    # written out before the file is opened, which empties it: a value json cannot hold then
    # leaves the file of the last run rather than nothing
    text = JSON.json(json_value(d.series))
    open(d.path, "w") do io
        print(io, text)
    end
end

# the handle a destination was given out as before destinations were objects, which
# `get_sim_file` still gives
handle(d::Union{TextFile, Stream}) = d.io
handle(d::Union{JsonFile, DataStore}) = d.series

"""
    Outputs

the destinations of a simulation, found by name the first time they are asked for and the
same object afterwards, the `Data` ones apart, and the formats of what is written. This is the
only place that reads a destination name: `"stdout"` (or `"-"`), `"stderr"` and `""` are the
streams of the process, a name ending in `.json` an accumulating file, and any other a text
file. When `redirect` is given, every name is that stream, a `Data` destination excepted,
which is not written anywhere.
"""
struct Outputs
    redirect::Union{Nothing, IO}
    files::Dict{String, Destination}
    data::Dict{String, Series}
    formats::Tuple{Printf.Format, Printf.Format}
    Outputs(redirect, time_format::String, data_format::String) =
        new(redirect, Dict{String, Destination}(), Dict{String, Series}(),
            (Printf.Format(time_format), Printf.Format(data_format)))
end

function open_destination(name::AbstractString)
    if name == "stdout" || name == "-"
        return Stream(stdout)
    elseif name == ""
        return Stream(devnull)
    elseif name == "stderr"
        return Stream(stderr)
    elseif last(splitext(name)) == ".json"
        return JsonFile(name, Series())
    end
    return TextFile(open(name, "w"))
end

destination(o::Outputs, name::AbstractString) =
    isnothing(o.redirect) ? get!(() -> open_destination(name), o.files, name) : Stream(o.redirect)

destination(o::Outputs, d::Data) = DataStore(get!(Series, o.data, d.name))

"""
    output_marks(::Outputs)

how far every destination has got, the part a checkpoint records so that it can write back
later exactly what they held now: the position of each text file and the lengths of each
series. The streams of the process have no position to keep.
"""
output_marks(o::Outputs) =
    (files = Dict{String, Any}(name => reached(d) for (name, d) in o.files if !(d isa Stream)),
     data = Dict{String, Dict{String, Int}}(name => series_lengths(s) for (name, s) in o.data))

# what a checkpoint carries of the destinations, up to where they had got
persist_outputs(o::Outputs, marks) =
    Dict("files" => Dict(name => persist(o.files[name], m) for (name, m) in marks.files),
         "data" => Dict(name => persist_series(o.data[name], m) for (name, m) in marks.data))

"""
    restore_outputs!(::Outputs, persisted)

put back the destinations a checkpoint carries: a text file is cut back to where it was and
continued, a series holds again what it held. A file the checkpoint does not know is not
touched here, and it is created, emptying it, the first time it is asked for, as the
uninterrupted run creates it.
"""
function restore_outputs!(o::Outputs, persisted)
    for (name, p) in persisted["files"]
        if haskey(p, "text")
            if isfile(name) && filesize(name) > p["text"]
                open(name, "a") do io
                    Base.truncate(io, p["text"])
                end
            end
            o.files[name] = TextFile(open(name, "a"))
        else
            o.files[name] = JsonFile(name, restore_series(p["json"]))
        end
    end
    for (name, s) in persisted["data"]
        o.data[name] = restore_series(s)
    end
    return o
end

close_outputs!(o::Outputs) = foreach(close!, values(o.files))
