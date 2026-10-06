# The destinations of measurements: text files, streams, json files and Data, kept in memory,
# the sinks they open to, how values are written to each, and how they are persisted and
# restored with a checkpoint.

export Data

"""
    Destination

what a destination of measurements is, as `output` and `get_sim_file` take it: a description,
compared by value, which the simulation opens to a sink the first time it is used. A name is
read as one by `Destination(name)`. The kinds are `TextFile`, `LogFile`, `Stream`, `JsonFile`
and `Data`.
"""
abstract type Destination end

"""
    Data(name)

a destination of measurements kept in memory: what is measured into `Data(name)` is gathered
in `sim.data[name]`, which `data_to_frame` turns into a table.

# Examples

    output(sim, Data("magnetization") => [X, Z(1)])
    data_to_frame(sim.data["magnetization"])
"""
struct Data <: Destination
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
    TextFile(name)

a text file of the simulation, written line by line as the measurements are made, known by
its normalized path, so that `"data"` and `"./data"` are one file rather than each emptying
the other
"""
struct TextFile <: Destination
    name::String
    TextFile(name::AbstractString) = new(normpath(name))
end

"""
    LogFile()

the log of the simulation, the file `log`: a text file that a resumed run continues without
cutting it back, keeping the history of every run
"""
struct LogFile <: Destination end

"""
    Stream(io)

a stream of the process, `stdout`, `stderr`, `devnull` or the stream `runTMS` redirects
everything to: it belongs to the process, so it is neither closed nor resumed
"""
struct Stream <: Destination
    io::IO
end

"""
    JsonFile(name)

a `.json` file, accumulated in memory and written when the files are closed, known by its
normalized path as a `TextFile` is
"""
struct JsonFile <: Destination
    name::String
    JsonFile(name::AbstractString) = new(normpath(name))
end

"""
    Destination(name)

the destination a name stands for: `"stdout"` (or `"-"`), `"stderr"` and `""` the streams of
the process, the last being `devnull`, `"log"` the log, a name ending in `.json` a json file,
any other a text file
"""
function Destination(name::AbstractString)
    if name == "stdout" || name == "-"
        return Stream(stdout)
    elseif name == ""
        return Stream(devnull)
    elseif name == "stderr"
        return Stream(stderr)
    elseif normpath(name) == "log"
        return LogFile()
    elseif last(splitext(name)) == ".json"
        return JsonFile(name)
    end
    return TextFile(name)
end

Destination(d::Destination) = d

############### sinks ###############

"""
    Sink

what a destination opens to, holding what is written to it: a `TextSink` for a text file or
the log, a `StreamSink` for a stream, a `JsonSink` or a `DataSink` for an accumulating
destination. Each kind has its methods of `emit!`, `new_event`, `reached`, `persist` and
`close!`; `emit_line!` and `stream` are for those of text.
"""
abstract type Sink end

"""
    TextSink(io)

the open file of a `TextFile` or of the `LogFile`, see `Sink`.
"""
struct TextSink <: Sink
    io::IO
end

"""
    StreamSink(io)

the stream of a `Stream`, see `Sink`.
"""
struct StreamSink <: Sink
    io::IO
end

"""
    JsonSink(path, series)

the series of a `JsonFile`, written to `path` when it is closed, see `Sink`.
"""
struct JsonSink <: Sink
    path::String
    series::Series
end

"""
    DataSink(series)

the series of a `Data`, the one `sim.data` holds, see `Sink`.
"""
struct DataSink <: Sink
    series::Series
end

"""
    emit!(::Sink, formats, time, values; event)

take the values of one call of `output`, pairs `header => value`: a text sink writes them as
rows and flushes, an accumulating sink appends them to its series as the event `event`, a
new one by default.
"""
function emit!(s::Union{TextSink, StreamSink}, formats, time, values; event = nothing)
    for (header, value) in values
        write_row(s.io, formats, time, header, value)
    end
    # so that a run can be followed, and a crash loses nothing already measured
    flush(s.io)
end

# `formats` named although unused: Julia 1.10 refuses `_` beside a keyword whose default is
# computed
function emit!(s::Union{JsonSink, DataSink}, formats, time, values;
               event = next_event(s.series))
    for (header, value) in values
        push_value!(s.series, header, time, value, event)
    end
end

"""
    emit_line!(::Sink, text)

write a line of the log to a text sink, flushed at once.
"""
function emit_line!(s::Union{TextSink, StreamSink}, text)
    println(s.io, text)
    flush(s.io)
end

"""
    new_event(::Sink)

the event the values of a call of `output` take in an accumulating sink, `nothing` for a
text sink, which does not number them
"""
new_event(::Union{TextSink, StreamSink}) = nothing
new_event(s::Union{JsonSink, DataSink}) = next_event(s.series)

"""
    reached(::Sink)

how far a sink has got, as a checkpoint records it: the position of a file, once flushed,
the lengths of the series of an accumulating sink, `nothing` for a stream, which is not
resumed.
"""
reached(s::TextSink) = (flush(s.io); position(s.io))
reached(::StreamSink) = nothing
reached(s::Union{JsonSink, DataSink}) = series_lengths(s.series)

"""
    persist(::Sink, reached)

what a checkpoint carries of a sink up to `reached`, as pairs: the position of a file,
already on the disk, the part of the series of an accumulating sink
"""
persist(::TextSink, pos::Int) = ("position" => pos,)
persist(s::Union{JsonSink, DataSink}, lengths) =
    ("series" => persist_series(s.series, lengths),)

"""
    close!(::Sink)

close a file, or write a json file; a stream or a `Data` is left as it is.
"""
close!(s::TextSink) = close(s.io)
close!(::Union{StreamSink, DataSink}) = nothing

function close!(s::JsonSink)
    # serialized before the file is opened, which empties it: a value json cannot hold leaves
    # the previous file
    text = JSON.json(json_value(s.series))
    open(s.path, "w") do io
        print(io, text)
    end
end

"""
    stream(::Sink)

the stream `get_sim_file` gives, that of a text sink: the series of an accumulating sink are
written by `output` alone, which numbers the events.
"""
stream(s::Union{TextSink, StreamSink}) = s.io
stream(::Union{JsonSink, DataSink}) =
    error("get_sim_file gives text files and streams: a json file or a Data is written by " *
          "output")

"""
    in_dir(dir, name)

the file `name` taken in the directory `dir` when it is relative, an empty `dir` standing
for the working directory
"""
in_dir(dir::AbstractString, name::AbstractString) =
    isempty(dir) || isabspath(name) ? name : joinpath(dir, name)

"""
    open_sink(::Destination, dir)

the sink a destination opens to, its file taken in the directory `dir`. A text file is
created, emptied if it exists.
"""
open_sink(d::TextFile, dir) = TextSink(open(in_dir(dir, d.name), "w"))
open_sink(::LogFile, dir) = TextSink(open(in_dir(dir, "log"), "w"))
open_sink(d::Stream, _) = StreamSink(d.io)
open_sink(::Data, _) = DataSink(Series())

function open_sink(d::JsonFile, dir)
    path = in_dir(dir, d.name)
    # refused now, as `open` refuses a text file, rather than when it is written at the end of
    # the run
    parent = dirname(path)
    if !isempty(parent) && !isdir(parent)
        error("cannot write $(d.name): there is no directory $parent")
    end
    return JsonSink(path, Series())
end

"""
    describe(::Destination)

how a checkpoint names a destination, as pairs, which `described` reads back
"""
describe(d::TextFile) = ("kind" => "text", "name" => d.name)
describe(::LogFile) = ("kind" => "log", "name" => "log")
describe(d::JsonFile) = ("kind" => "json", "name" => d.name)
describe(d::Data) = ("kind" => "data", "name" => d.name)

"""
    described(entry)

the destination an entry of a checkpoint names, see `describe`
"""
function described(p::AbstractDict)
    kind, name = p["kind"], p["name"]
    if kind == "text"
        return TextFile(name)
    elseif kind == "log"
        return LogFile()
    elseif kind == "json"
        return JsonFile(name)
    elseif kind == "data"
        return Data(name)
    end
    error("a checkpoint names a destination of unknown kind $kind")
end

"""
    restore_sink(::Destination, dir, entry)

the sink of a destination put back from its entry in a checkpoint, with whether a text file
is shorter than its recorded position: a text file is cut back to that position and
continued, the log continued as it is, a series holds again what it held.
"""
function restore_sink(d::TextFile, dir, p)
    path = in_dir(dir, d.name)
    len = isfile(path) ? filesize(path) : 0
    if len > p["position"]
        open(path, "a") do io
            Base.truncate(io, p["position"])
        end
    end
    return TextSink(open(path, "a")), len < p["position"]
end

restore_sink(::LogFile, dir, _) = TextSink(open(in_dir(dir, "log"), "a")), false
restore_sink(d::JsonFile, dir, p) =
    JsonSink(in_dir(dir, d.name), restore_series(p["series"])), false
restore_sink(::Data, _, p) = DataSink(restore_series(p["series"])), false

############### the destinations of a simulation ###############

"""
    Outputs(redirect, time_format, data_format[, dir])

the sinks of a simulation, by destination, and the formats of what is written. A destination
is opened the first time it is used and is the same afterwards, its file taken in the
directory `dir` when it is relative. With `redirect`, every destination is that stream, a
`Data` still being kept in memory.
"""
struct Outputs
    redirect::Union{Nothing, StreamSink}
    sinks::Dict{Destination, Sink}
    formats::Tuple{Printf.Format, Printf.Format}
    dir::String
    Outputs(redirect, time_format::String, data_format::String, dir::String = "") =
        new(isnothing(redirect) ? nothing : StreamSink(redirect), Dict{Destination, Sink}(),
            (Printf.Format(time_format), Printf.Format(data_format)), dir)
end

"""
    sink(::Outputs, ::Destination)

the sink of a destination, opened on first use, or the redirect stream if there is one
"""
sink(o::Outputs, d::Destination) =
    isnothing(o.redirect) ? get!(() -> open_sink(d, o.dir), o.sinks, d) : o.redirect

sink(o::Outputs, d::Data) = get!(() -> open_sink(d, o.dir), o.sinks, d)

"""
    data_series(::Outputs)

the series of the `Data` destinations by name, as `sim.data` gives them
"""
data_series(o::Outputs) =
    Dict{String, Series}(d.name => s.series for (d, s) in o.sinks if d isa Data)

"""
    output_marks(::Outputs)

how far every sink but the streams has got, by destination, as a checkpoint records it, see
`reached`
"""
function output_marks(o::Outputs)
    marks = Dict{Destination, Any}()
    for (d, s) in o.sinks
        m = reached(s)
        if !isnothing(m)
            marks[d] = m
        end
    end
    return marks
end

"""
    persist_outputs(::Outputs, marks)

what a checkpoint carries of the destinations, up to the `marks` of `output_marks`: one entry
per destination, what `describe` and `persist` give.
"""
persist_outputs(o::Outputs, marks) =
    [ Dict(describe(d)..., persist(o.sinks[d], m)...) for (d, m) in marks ]

"""
    restore_outputs!(::Outputs, persisted)

put back the destinations a checkpoint carries, see `restore_sink`, a destination the
checkpoint does not know being created on first use.

Return the names of the text files shorter than their recorded position, cut or removed by
hand, which are continued as they are, for the log to say so.
"""
function restore_outputs!(o::Outputs, persisted)
    shortened = String[]
    for p in persisted
        d = described(p)
        s, short = restore_sink(d, o.dir, p)
        o.sinks[d] = s
        if short
            push!(shortened, d.name)
        end
    end
    return shortened
end

"""
    close_outputs!(::Outputs)

close every sink, which writes the json files, each even when another fails, the first
failure being raised at the end
"""
function close_outputs!(o::Outputs)
    failure = nothing
    for s in values(o.sinks)
        try
            close!(s)
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
