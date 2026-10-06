# The destinations of measurements: text files, streams, json files and Data, kept in memory,
# how values are written to each, and how they are persisted and restored with a checkpoint.

export Data

"""
    Data(name)

a destination of measurements kept in memory: what is measured into `Data(name)` is gathered
in `sim.data[name]`, which `data_to_frame` turns into a table.

# Examples

    output(sim, Data("magnetization") => [X, Z(1)])
    data_to_frame(sim.data["magnetization"])
"""
struct Data
    name::String
end

############### the values of an accumulating destination ###############

"""
    json_value(x)

`x` as a json destination writes it: a complex number as `{"re": …, "im": …}` at any depth, a
matrix as the list of its rows, a float that is not finite as the string `"Inf"`, `"-Inf"` or
`"NaN"`
"""
json_value(x::AbstractFloat) = isfinite(x) ? x : string(x)
json_value(x::Complex) = Dict("re" => json_value(real(x)), "im" => json_value(imag(x)))
json_value(x::AbstractMatrix) = [ json_value.(x[i, :]) for i in axes(x, 1) ]
json_value(x::AbstractArray) = map(json_value, x)
json_value(x::AbstractDict) = Dict(k => json_value(v) for (k, v) in x)
json_value(x) = x

"""
    checkpoint_value(x)

`x` encoded for a checkpoint, which `restored_value` decodes: complex numbers, matrices,
floats that are not finite and integers other than `Int`, which json does not give back as
they are, are wrapped in a dictionary saying what they are.
"""
checkpoint_value(x::AbstractFloat) = isfinite(x) ? x : Dict("float" => string(x))
checkpoint_value(x::Union{Base.BitInteger, BigInt}) = Dict("integer" => string(x), "type" => string(typeof(x)))
checkpoint_value(x::Int) = x
checkpoint_value(x::Complex) = Dict("complex" => [checkpoint_value(real(x)), checkpoint_value(imag(x))])
checkpoint_value(x::AbstractMatrix) = Dict("matrix" => [ checkpoint_value(x[i, :]) for i in axes(x, 1) ])
checkpoint_value(x::AbstractArray) = map(checkpoint_value, x)
checkpoint_value(x::AbstractDict) = Dict(k => checkpoint_value(v) for (k, v) in x)
checkpoint_value(x) = x

"""
    checkpoint_integers

the integer types a checkpoint writes with their name, see `checkpoint_value`
"""
const checkpoint_integers = Dict(string(T) => T for T in (Int8, Int16, Int32, Int64, Int128, UInt8,
                                                           UInt16, UInt32, UInt64, UInt128, BigInt))

"""
    restored_value(x)

a value read back from a checkpoint, decoded from what `checkpoint_value` wrote, an array
taking back the element type of its values
"""
function restored_value(x::AbstractDict)
    if haskey(x, "complex")
        r, i = restored_value.(x["complex"])
        return complex(r, i)
    elseif haskey(x, "float")
        return parse(Float64, x["float"])
    elseif haskey(x, "integer")
        return parse(checkpoint_integers[x["type"]], x["integer"])
    elseif haskey(x, "matrix")
        return stack(restored_value.(x["matrix"]); dims = 1)
    end
    return Dict(k => restored_value(v) for (k, v) in x)
end
restored_value(x::AbstractVector) = map(restored_value, x)
restored_value(x) = x

"""
    Series

what an accumulating destination, a json file or a `Data` one, holds, and `data_to_frame`
reads: for each header, the times, the values and the calls of `output` they came from, as
`Dict("times" => [...], "data" => [...], "events" => [...])`. A series only grows, so a
checkpoint records it by its lengths.
"""
const Series = Dict{String, Dict}

"""
    next_event(::Series)

the number of the next call of `output` on a series, which tells values measured together
from values that only share a time
"""
next_event(s::Series) = 1 + maximum((last(d["events"]) for d in values(s)); init = 0)

"""
    push_value!(::Series, header, time, value, event)

append to the series of `header`, opened if needed, a value measured at `time` by the call
`event` of `output`.
"""
function push_value!(s::Series, header, time, value, event::Int)
    d = get!(s, header, Dict("times" => [], "data" => [], "events" => Int[]))
    push!(d["times"], time)
    push!(d["data"], value)
    push!(d["events"], event)
    return nothing
end

"""
    series_lengths(::Series)

the number of values held under each header, which marks a series in a checkpoint
"""
series_lengths(s::Series) = Dict{String, Int}(h => length(d["times"]) for (h, d) in s)

"""
    persist_series(::Series, lengths)

the part of a series that `lengths` covers, encoded for a checkpoint, a header opened since
being left out
"""
persist_series(s::Series, mark::Dict{String, Int}) =
    Dict(h => Dict(k => checkpoint_value(s[h][k][1:n]) for k in ("times", "data", "events"))
         for (h, n) in mark)

"""
    restore_series(d)

the series a checkpoint carried, decoded, its times and values in vectors of element type
`Any`, onto which the measurements to come are pushed
"""
restore_series(d) =
    Series(h => Dict("times" => Any[ restored_value(x) for x in s["times"] ],
                     "data" => Any[ restored_value(x) for x in s["data"] ],
                     "events" => Int[ x for x in s["events"] ]) for (h, s) in d)

############### the rows of a text destination ###############

"""
    complex_columns(put, io, z, format)

write the complex number `z` as two columns, real then imaginary part, each written by
`put(io, part, format)`
"""
function complex_columns(put, file, z, format)
    put(file, real(z), format)
    print(file, "\t")
    put(file, imag(z), format)
end

"""
    output_one(io, x, format)

write a measured value in a row of a text destination: a float with `format`, a complex
number as two columns, an array number by number, a matrix row by row, anything else as
`print` does.
"""
function output_one(file, x::AbstractFloat, format)
    Printf.format(file, format, x)
end

output_one(file, x::Complex, format) = complex_columns(output_one, file, x, format)

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

"""
    output_time(io, t, format)

write the time column of a row, always with the time format, unlike a measured value that is
not a float
"""
function output_time(file, t::Number, format)
    Printf.format(file, format, t)
end

# `Printf` refuses a complex number
output_time(file, t::Complex, format) = complex_columns(output_time, file, t, format)

"""
    write_row(io, formats, time, header, data)

write a value as rows of a text destination, `formats` being the time and data formats: the
header, the time and the numbers of the value, separated by tabs. A matrix takes a line with
the header alone, then one row per line `l` of the matrix, headed `header:l`.
"""
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

where the measurements of an `output` go, of five kinds:

- `TextFile`: a file of the simulation, written line by line as the measurements are made;
- `LogFile`: the file `log`, a text file that a resumed run continues without cutting it back;
- `Stream`: `stdout`, `stderr`, `devnull` or the stream `runTMS` redirects everything to;
  it belongs to the process, so it is neither closed nor resumed;
- `JsonFile`: a `.json` file, accumulated in memory and written when the files are closed;
- `DataStore`: a `Data(name)` destination, accumulated in memory in `sim.data`.

Each kind has its methods of `emit!`, `new_event`, `reached`, `persist`, `close!` and
`handle`; `emit_line!` is for those that can hold the log.
"""
abstract type Destination end

"""
    TextFile(io)

a text file of the simulation, see `Destination`.
"""
struct TextFile <: Destination
    io::IO
end

"""
    LogFile(io)

the log of the simulation, see `Destination`.
"""
struct LogFile <: Destination
    io::IO
end

"""
    Stream(io)

a stream of the process, see `Destination`.
"""
struct Stream <: Destination
    io::IO
end

"""
    JsonFile(path, series)

a `.json` file, whose series are written to `path` when it is closed, see `Destination`.
"""
struct JsonFile <: Destination
    path::String
    series::Series
end

"""
    DataStore(series)

a `Data` destination, whose series is the one `sim.data` holds, see `Destination`.
"""
struct DataStore <: Destination
    series::Series
end

"""
    emit!(::Destination, formats, time, values; event)

take the values of one call of `output`, pairs `header => value`: a text destination writes
them as rows and flushes, an accumulating destination appends them to its series as the
event `event`, a new one by default.
"""
function emit!(d::Union{TextFile, LogFile, Stream}, formats, time, values; event = nothing)
    for (header, value) in values
        write_row(d.io, formats, time, header, value)
    end
    # so that a run can be followed, and a crash loses nothing already measured
    flush(d.io)
end

# `formats` named although unused: Julia 1.10 refuses `_` beside a keyword whose default is
# computed
function emit!(d::Union{JsonFile, DataStore}, formats, time, values;
               event = next_event(d.series))
    for (header, value) in values
        push_value!(d.series, header, time, value, event)
    end
end

"""
    emit_line!(::Destination, text)

write a line of the log to the log file or to a stream, flushed at once.
"""
function emit_line!(d::Union{LogFile, Stream}, text)
    println(d.io, text)
    flush(d.io)
end

"""
    new_event(::Destination)

the event the values of a call of `output` take in an accumulating destination, `nothing`
for a text destination, which does not number them
"""
new_event(::Union{TextFile, LogFile, Stream}) = nothing
new_event(d::Union{JsonFile, DataStore}) = next_event(d.series)

"""
    reached(::Destination)

how far a destination has got, as a checkpoint records it: the position of a file, once
flushed, the lengths of the series of an accumulating destination, `nothing` for a stream,
which is not resumed.
"""
reached(d::Union{TextFile, LogFile}) = (flush(d.io); position(d.io))
reached(::Stream) = nothing
reached(d::Union{JsonFile, DataStore}) = series_lengths(d.series)

"""
    persist(::Destination, reached)

what a checkpoint carries of a destination up to `reached`, as pairs headed by its `kind`:
the position of a text file, already on the disk, which a resume cuts it back to, nothing
more for the log, the part of the series of an accumulating destination.
"""
persist(::TextFile, pos::Int) = ("kind" => "text", "position" => pos)
persist(::LogFile, _) = ("kind" => "log",)
persist(d::JsonFile, lengths) =
    ("kind" => "json", "series" => persist_series(d.series, lengths))
persist(d::DataStore, lengths) =
    ("kind" => "data", "series" => persist_series(d.series, lengths))

"""
    close!(::Destination)

close a file, or write a json file; a stream or a `Data` destination is left as it is.
"""
close!(d::Union{TextFile, LogFile}) = close(d.io)
close!(::Union{Stream, DataStore}) = nothing

function close!(d::JsonFile)
    # serialized before the file is opened, which empties it: a value json cannot hold leaves
    # the previous file
    text = JSON.json(json_value(d.series))
    open(d.path, "w") do io
        print(io, text)
    end
end

"""
    handle(::Destination)

what `get_sim_file` returns for a destination: the stream of a text destination, the series
of an accumulating one.
"""
handle(d::Union{TextFile, LogFile, Stream}) = d.io
handle(d::Union{JsonFile, DataStore}) = d.series

"""
    in_dir(dir, name)

the file `name` taken in the directory `dir` when it is relative, an empty `dir` standing
for the working directory
"""
in_dir(dir::AbstractString, name::AbstractString) =
    isempty(dir) || isabspath(name) ? name : joinpath(dir, name)

"""
    Outputs(redirect, time_format, data_format[, dir])

the destinations of a simulation and the formats of what is written. A destination is opened
the first time it is asked for and is the same afterwards, kept under its key: the
normalized path of a name, or the `Data` itself. `"stdout"` (or `"-"`), `"stderr"` and `""`
(`devnull`) are the streams of the process, `"log"` the log, a name ending in `.json` an
accumulating file, any other a text file, taken in the directory `dir` when it is relative.
With `redirect`, every name is that stream, a `Data` destination still being kept in memory.
"""
struct Outputs
    redirect::Union{Nothing, Stream}
    destinations::Dict{Union{String, Data}, Destination}
    formats::Tuple{Printf.Format, Printf.Format}
    dir::String
    Outputs(redirect, time_format::String, data_format::String, dir::String = "") =
        new(isnothing(redirect) ? nothing : Stream(redirect),
            Dict{Union{String, Data}, Destination}(),
            (Printf.Format(time_format), Printf.Format(data_format)), dir)
end

"""
    open_destination(dir, name)

the destination a name stands for, see `Outputs`. A text file is created, emptied if it
exists.
"""
function open_destination(dir::AbstractString, name::AbstractString)
    if name == "stdout" || name == "-"
        return Stream(stdout)
    elseif name == ""
        return Stream(devnull)
    elseif name == "stderr"
        return Stream(stderr)
    end
    path = in_dir(dir, name)
    if last(splitext(name)) == ".json"
        # refused now, as `open` refuses a text file, rather than when it is written at the end
        # of the run
        parent = dirname(path)
        if !isempty(parent) && !isdir(parent)
            error("cannot write $name: there is no directory $parent")
        end
        return JsonFile(path, Series())
    end
    io = open(path, "w")
    return normpath(name) == "log" ? LogFile(io) : TextFile(io)
end

"""
    destination(::Outputs, name)
    destination(::Outputs, ::Data)

the destination of the given name, opened on first use, or the redirect stream if there is
one. A `Data` destination holds the series of that name in `sim.data`, created if needed.

A file is known by its normalized path, so that `"data"` and `"./data"` share one destination
rather than each emptying the other; the empty name is kept as it is.
"""
destination(o::Outputs, name::AbstractString) =
    isnothing(o.redirect) ?
        get!(() -> open_destination(o.dir, name), o.destinations,
             isempty(name) ? String(name) : normpath(name)) :
        o.redirect

destination(o::Outputs, d::Data) = get!(() -> DataStore(Series()), o.destinations, d)

"""
    data_series(::Outputs)

the series of the `Data` destinations by name, as `sim.data` gives them
"""
data_series(o::Outputs) =
    Dict{String, Series}(key.name => handle(d) for (key, d) in o.destinations if key isa Data)

"""
    output_marks(::Outputs)

how far every destination but the streams has got, by key, as a checkpoint records it, see
`reached`
"""
function output_marks(o::Outputs)
    marks = Dict{Union{String, Data}, Any}()
    for (key, d) in o.destinations
        m = reached(d)
        if !isnothing(m)
            marks[key] = m
        end
    end
    return marks
end

"""
    persist_outputs(::Outputs, marks)

what a checkpoint carries of the destinations, up to the `marks` of `output_marks`: one entry
per destination, its name and what `persist` gives.
"""
persist_outputs(o::Outputs, marks) =
    [ Dict("name" => key isa Data ? key.name : key, persist(o.destinations[key], m)...)
      for (key, m) in marks ]

"""
    restored_destination(dir, entry)

the key and the destination an entry of `persist_outputs` stands for, with whether a text
file is shorter than its recorded position: a text file is cut back to that position and
continued, the log continued as it is, keeping the history of every run, a series holds again
what it held.
"""
function restored_destination(dir::AbstractString, p::AbstractDict)
    name, kind = p["name"], p["kind"]
    if kind == "data"
        return Data(name), DataStore(restore_series(p["series"])), false
    end
    path = in_dir(dir, name)
    if kind == "json"
        return name, JsonFile(path, restore_series(p["series"])), false
    elseif kind == "log"
        return name, LogFile(open(path, "a")), false
    end
    len = isfile(path) ? filesize(path) : 0
    if len > p["position"]
        open(path, "a") do io
            Base.truncate(io, p["position"])
        end
    end
    return name, TextFile(open(path, "a")), len < p["position"]
end

"""
    restore_outputs!(::Outputs, persisted)

put back the destinations a checkpoint carries, see `restored_destination`, a destination the
checkpoint does not know being created on first use.

Return the names of the text files shorter than their recorded position, cut or removed by
hand, which are continued as they are, for the log to say so.
"""
function restore_outputs!(o::Outputs, persisted)
    shortened = String[]
    for p in persisted
        key, d, short = restored_destination(o.dir, p)
        o.destinations[key] = d
        if short
            push!(shortened, key)
        end
    end
    return shortened
end

"""
    close_outputs!(::Outputs)

close every destination, which writes the json files, each even when another fails, the
first failure being raised at the end
"""
function close_outputs!(o::Outputs)
    failure = nothing
    for d in values(o.destinations)
        try
            close!(d)
        catch e
            if isnothing(failure)
                failure = e
            end
        end
    end
    if !isnothing(failure)
        throw(failure)
    end
    return nothing
end
