# The destinations of measurements: text files, streams, json files and Data, kept in memory,
# how values are written to each, and how they are persisted and restored with a checkpoint.

export Data

"""
    Data(name)

a destination of measurements kept in memory rather than written to a file: what is measured
into `Data(name)` is gathered in `sim.data[name]`, which `data_to_frame` turns into a table.

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

`x` as a json destination writes it: a complex number as `{"re": …, "im": …}` at any depth,
rather than as JSON.jl would, which depends on its version; a matrix as the list of its rows,
as a text file writes it; a float that is not finite as the string Julia prints, `"Inf"`,
`"-Inf"` or `"NaN"`, json having no number for it.
"""
json_value(x::AbstractFloat) = isfinite(x) ? x : string(x)
json_value(x::Complex) = Dict("re" => json_value(real(x)), "im" => json_value(imag(x)))
json_value(x::AbstractMatrix) = [ json_value.(x[i, :]) for i in axes(x, 1) ]
json_value(x::AbstractArray) = map(json_value, x)
json_value(x::AbstractDict) = Dict(k => json_value(v) for (k, v) in x)
json_value(x) = x

"""
    checkpoint_value(x)

`x` encoded for a checkpoint, which `restored_value` decodes. Json holds neither complex
numbers, nor matrices (one comes back as the vector of its columns), nor floats that are not
finite, nor integers other than `Int` (one comes back as an `Int64` or a `BigInt`, and with
JSON 0.21 above `typemax(Int64)` as a wrong number), so these are wrapped in a dictionary
saying what they are: a resumed run then hands back the values an uninterrupted one would.
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

a value read back from a checkpoint, decoded from what `checkpoint_value` wrote. An array is
given back the element type of its values, which reading json loses.
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

what an accumulating destination, a json file or a `Data` one, holds: for each header, the
times, the values and the calls of `output` they came from, as
`Dict("times" => [...], "data" => [...], "events" => [...])`. It is what `sim.data` holds for
each `Data` destination, and what `data_to_frame` reads.

A series only grows, so a checkpoint records it by its lengths and can later write back
exactly the part it held then.
"""
const Series = Dict{String, Dict}

"""
    next_event(::Series)

the number of the next call of `output` on a series. Each value records the call it came
from, so that values measured together can be told from values that only share a time: the
time repeats over the sweeps of dmrg, in a circuit, or when a phase sets it back.
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

the number of values held under each header, which is how a checkpoint marks a series.
"""
series_lengths(s::Series) = Dict{String, Int}(h => length(d["times"]) for (h, d) in s)

"""
    persist_series(::Series, lengths)

the part of a series that `lengths` covers, encoded for a checkpoint. A header opened since is
left out, the uninterrupted run not having opened it at that point either.
"""
persist_series(s::Series, mark::Dict{String, Int}) =
    Dict(h => Dict(k => checkpoint_value(s[h][k][1:n]) for k in ("times", "data", "events"))
         for (h, n) in mark)

"""
    restore_series(d)

the series a checkpoint carried, decoded. The times and values stay vectors of element type
`Any`, whatever the values read back, since the measurements to come are pushed onto them.
"""
restore_series(d) =
    Series(h => Dict("times" => Any[ restored_value(x) for x in s["times"] ],
                     "data" => Any[ restored_value(x) for x in s["data"] ],
                     "events" => Int[ x for x in s["events"] ]) for (h, s) in d)

############### the rows of a text destination ###############

"""
    complex_columns(put, io, z, format)

write the complex number `z` as two columns of a text destination, its real part then its
imaginary part, each written by `put(io, part, format)`
"""
function complex_columns(put, file, z, format)
    put(file, real(z), format)
    print(file, "\t")
    put(file, imag(z), format)
end

"""
    output_one(io, x, format)

write a measured value in a row of a text destination: a float with `format`, a complex
number as two columns, real part then imaginary part, a vector or a matrix number by number,
a matrix row by row, rather than as the literal Julia prints, anything else as `print` does.
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

write the time column of a row. It always takes the time format, unlike a measured value,
which keeps its own printed form when it is not a float: `Linkdim` reads 8, not 8.000.
"""
function output_time(file, t::Number, format)
    Printf.format(file, format, t)
end

# `Printf` refuses a complex number outright, so a complex simulation time is written as
# the two columns a complex measurement takes, real then imaginary
output_time(file, t::Complex, format) = complex_columns(output_time, file, t, format)

"""
    write_row(io, formats, time, header, data)

write a value as rows of a text destination, `formats` being the time and data formats: the
header, the time and the numbers of the value, separated by tabs. A matrix takes a line with
the header alone, then one row per line of the matrix, headed `header:l`.
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

where the measurements of an `output` go, of four kinds:

- `TextFile`: a file of the simulation, written line by line as the measurements are made;
- `Stream`: `stdout`, `stderr`, `devnull` or the stream `runTMS` was told to redirect
  everything to; it belongs to the process, so it is neither closed nor resumed;
- `JsonFile`: a `.json` file, accumulated in memory and written when the files are closed;
- `DataStore`: a `Data(name)` destination, accumulated in memory and handed over in
  `sim.data`.

Each kind answers those of `emit!`, `emit_line!`, `reached`, `persist` and `close!` that
apply to it.
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
    emit!(::Destination, formats, time, values)

take the values of one call of `output`, pairs `header => value`: a text file or a stream
writes them as rows and flushes, an accumulating destination appends them to its series as
one new event.
"""
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

"""
    emit_line!(::Destination, text)

write a line of the log to a text file or a stream, flushed at once.
"""
function emit_line!(d::Union{TextFile, Stream}, text)
    println(d.io, text)
    flush(d.io)
end

"""
    reached(::Destination)

how far a text file or an accumulating destination has got: the position of the file, once
flushed, or the lengths of the series.
"""
reached(d::TextFile) = (flush(d.io); position(d.io))
reached(d::Union{JsonFile, DataStore}) = series_lengths(d.series)

"""
    persist(::Destination, reached)

what a checkpoint carries of a file up to `reached`. A text file is on the disk up to there
already, so only its position, which a resume cuts it back to; a json file is only written
at the end, so the part of its series.
"""
persist(::TextFile, pos::Int) = Dict("text" => pos)
persist(d::JsonFile, lengths) = Dict("json" => persist_series(d.series, lengths))

"""
    close!(::Destination)

close a text file, or write a json file; a stream or a `Data` destination is left as it is.
"""
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

"""
    handle(::Destination)

what `get_sim_file` returns for a destination: the stream of a text file or a stream, the
series of an accumulating destination.
"""
handle(d::Union{TextFile, Stream}) = d.io
handle(d::Union{JsonFile, DataStore}) = d.series

"""
    Outputs(redirect, time_format, data_format)

the destinations of a simulation and the formats of what is written. A destination is opened
the first time its name is asked for and is the same afterwards. This is the only place that
reads a destination name: `"stdout"` (or `"-"`), `"stderr"` and `""` (`devnull`) are the
streams of the process, a name ending in `.json` an accumulating file, any other a text file.
With `redirect`, every name is that stream, and a `Data` destination is still kept in memory.
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

"""
    open_destination(name)

the destination a name stands for, see `Outputs`. A text file is created, emptied if it
exists.
"""
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

"""
    destination(::Outputs, name)
    destination(::Outputs, ::Data)

the destination of the given name, opened on first use, or the redirect stream if there is
one. A `Data` destination holds the series of that name in `sim.data`, created if needed.
"""
destination(o::Outputs, name::AbstractString) =
    isnothing(o.redirect) ? get!(() -> open_destination(name), o.files, name) : Stream(o.redirect)

destination(o::Outputs, d::Data) = DataStore(get!(Series, o.data, d.name))

"""
    output_marks(::Outputs)

how far every destination has got, which a checkpoint records so as to write back later
exactly what they hold now: the position of each text file and the lengths of each series.
The streams of the process are left out, having no position to keep.
"""
output_marks(o::Outputs) =
    (files = Dict{String, Any}(name => reached(d) for (name, d) in o.files if !(d isa Stream)),
     data = Dict{String, Dict{String, Int}}(name => series_lengths(s) for (name, s) in o.data))

"""
    persist_outputs(::Outputs, marks)

what a checkpoint carries of the destinations, up to the `marks` of `output_marks`.
"""
persist_outputs(o::Outputs, marks) =
    Dict("files" => Dict(name => persist(o.files[name], m) for (name, m) in marks.files),
         "data" => Dict(name => persist_series(o.data[name], m) for (name, m) in marks.data))

"""
    restore_outputs!(::Outputs, persisted)

put back the destinations a checkpoint carries: a text file is cut back to its recorded
position and continued, a series holds again what it held. A file the checkpoint does not
know is left alone here, and is created, emptied, on first use, as in the uninterrupted run.
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

"""
    close_outputs!(::Outputs)

close every destination opened by name, which writes the json files, see `close!`.
"""
close_outputs!(o::Outputs) = foreach(close!, values(o.files))
