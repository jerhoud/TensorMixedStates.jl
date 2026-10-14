# The outputs of the measurements: Destination, the open type of where they are written, which
# other packages extend, the text tables, the JSON Lines files and the Data of TMS, and Outputs,
# which measures lists on a state and writes each to its destination.

export Destination, Data, data_to_frame

"""
    Destination

where measurements are written: a description, compared by value, opened once and given a row
at each call. TMS writes text tables, to files and to the streams of the process, JSON Lines
files, all named by a string, and `Data`. A package adds its own by a subtype and these methods
of TensorMixedStates, `o` being the destination opened, `o = open(d)`:

- `write_row(o, time, row)`: write `row`, the pairs `name => value` that `measure` returns,
  measured at `time`
- `close(o)`
- `check_list(o, names)`, optional: refuse a list of measurements, by their names, that it
  cannot take, before anything is measured
- `resume_mark(o)` and `open(d, mark)`, optional: a mark that a resumed simulation reopens the
  destination at, `nothing` by default, the destination then starting anew
"""
abstract type Destination end

"""
    write_row(o, time, row)

write to the destination opened `o` the values `row`, pairs `name => value`, measured at
`time`, see `Destination`
"""
function write_row end

"""
    check_list(o, names)

refuse a list of measurements, by their names, that the destination opened `o` cannot take,
see `Destination`; it takes any by default
"""
check_list(_, _) = nothing

"""
    resume_mark(o)

the mark that a resumed simulation reopens the destination opened `o` at, `open(d, mark)`, see
`Destination`; `nothing` by default, the destination then starting anew
"""
resume_mark(_) = nothing

open(d::Destination, _) = open(d)

"""
    prepare_list!(o, m)

check that the destination opened `o` can take the `Measure` `m`, see `check_list`, a text
table keeping what it needs to name its columns
"""
prepare_list!(o, m::Measure) = check_list(o, first.(m.entries))

"""
    put_row!(o, time, m, row)

write the values `row` of the `Measure` `m` to the destination opened `o`, see `write_row`
"""
put_row!(o, time, ::Measure, row) = write_row(o, time, row)

"""
    TextFile(path, format)
    Stream(io, format)

a text table written to the file `path` or to the stream `io` of the process, its numbers in
the Printf `format`
"""
struct TextFile <: Destination
    path::String
    format::String
end

struct Stream <: Destination
    io::IO
    format::String
end

show(io::IO, d::TextFile) = print(io, d.path)
show(io::IO, d::Stream) = print(io, d.io)

"""
    TextTable

a text table opened: its stream, the format of its numbers, the path of its file, `nothing`
for a stream of the process, the list of measurements a file holds, and the columns of the
times and of each measurement, those of its last header
"""
mutable struct TextTable
    io::IO
    format::Printf.Format
    path::Union{Nothing, String}
    names::Union{Nothing, Vector{String}}
    columns::Union{Nothing, Vector{Vector{String}}}
end

open(d::TextFile) = TextTable(open(d.path, "w"), Printf.Format(d.format), d.path, nothing, nothing)

"""
    reopened(path, position)

the file `path` opened at `position`, what was written after it removed; a file shorter than
that, cut or removed by hand, is continued as it is, with a warning
"""
function reopened(path::String, position::Integer)
    io = open(path, isfile(path) ? "r+" : "w")
    if filesize(path) < position
        @warn "$path is shorter than at the checkpoint, the lines it lost are not written again"
    else
        truncate(io, position)
    end
    seekend(io)
    return io
end

open(d::TextFile, mark::AbstractDict) =
    TextTable(reopened(d.path, mark["position"]), Printf.Format(d.format), d.path,
              Vector{String}(mark["names"]), [ Vector{String}(c) for c in mark["columns"] ])

open(d::Stream) = TextTable(d.io, Printf.Format(d.format), nothing, nothing, nothing)

function close(t::TextTable)
    if isnothing(t.path)
        flush(t.io)
    else
        close(t.io)
    end
    return nothing
end

resume_mark(t::TextTable) =
    if isnothing(t.path)
        nothing
    else
        flush(t.io)
        Dict("position" => position(t.io), "names" => t.names, "columns" => t.columns)
    end

# a file holds one list, a stream takes any, writing a header when it changes
function prepare_list!(t::TextTable, m::Measure)
    names = first.(m.entries)
    if isnothing(t.path) || t.names == names
        return nothing
    elseif !isnothing(t.names)
        error("$(t.path) holds the list $(join(t.names, ", ")), and another needs another file")
    end
    t.names = names
    return nothing
end

"""
    cells!(out, col, x, format)

the cells of the value `x` under the column `col`, pushed to `out` as pairs `column => text`:
a real number in `format`, an integer as it is, a complex number as `Re(col)` and `Im(col)`,
`missing` or `nothing` as an empty cell, a vector as `col[i]`, a matrix row by row as
`col[i,j]`, a named tuple as `col.field`, anything else as it prints, on one line
"""
cells!(out, col, x::Integer, _) = push!(out, col => string(x))
cells!(out, col, x::Real, f) = push!(out, col => Printf.format(f, float(x)))

function cells!(out, col, x::Complex, f)
    cells!(out, "Re($col)", real(x), f)
    cells!(out, "Im($col)", imag(x), f)
end

cells!(out, col, ::Union{Missing, Nothing}, _) = push!(out, col => "")

function cells!(out, col, x::AbstractVector, f)
    for (i, y) in enumerate(x)
        cells!(out, "$col[$i]", y, f)
    end
end

function cells!(out, col, x::AbstractMatrix, f)
    for i in axes(x, 1), j in axes(x, 2)
        cells!(out, "$col[$i,$j]", x[i, j], f)
    end
end

function cells!(out, col, x::NamedTuple, f)
    for (k, y) in pairs(x)
        cells!(out, "$col.$k", y, f)
    end
end

cells!(out, col, x, _) = push!(out, col => replace(string(x), r"[\t\r\n]" => " "))

"""
    entry_cells(name, m, x, format)

the cells of the value `x` of the measurement `m` named `name`, see `cells!`: an operator of one
site on site `i` under `name(i)`, a pair on sites `i` and `j` under `A(i)B(j)`, or `name(i,j)`
when it is named, and a measurement of the state named `φ` under the cells it has on its own,
prefixed by `φ.`
"""
function entry_cells(name::String, m, x, f)
    out = Pair{String, String}[]
    if m isa OnNamed
        prefix = "$(m.state)."
        inner = entry_cells(String(chopprefix(name, prefix)), m.entry, x, f)
        return [ prefix * c => v for (c, v) in inner ]
    elseif m isa SimpleOp && x isa AbstractVector
        for (i, y) in enumerate(x)
            cells!(out, "$name($i)", y, f)
        end
    elseif m isa Tuple{SimpleOp, SimpleOp} && x isa AbstractMatrix
        a, b = map(obs_name, m)
        col(i, j) = name == entry_name(m) ? "$a($i)$b($j)" : "$name($i,$j)"
        for i in axes(x, 1), j in axes(x, 2)
            cells!(out, col(i, j), x[i, j], f)
        end
    else
        cells!(out, name, x, f)
    end
    return out
end

function put_row!(t::TextTable, time, m::Measure, row)
    times = Pair{String, String}[]
    cells!(times, "time", time, t.format)
    parts = [ [times] ;
              [ entry_cells(name, e, x, t.format) for ((name, e), (_, x)) in zip(m.entries, row) ] ]
    columns = [ first.(p) for p in parts ]
    if columns ≠ t.columns
        if !isnothing(t.path) && !isnothing(t.columns)
            k = findfirst(k -> columns[k] ≠ t.columns[k], eachindex(columns))
            what = k == 1 ? "the time does" : "the values of $(first(row[k - 1])) do"
            error("$what not fit the columns $(t.path) took at its first line")
        end
        t.columns = columns
        println(t.io, join(reduce(vcat, columns), '\t'))
    end
    println(t.io, join([ last(c) for p in parts for c in p ], '\t'))
    flush(t.io)
    return nothing
end

"""
    JsonFile(path)

a JSON Lines file, one object per row
"""
struct JsonFile <: Destination
    path::String
end

show(io::IO, d::JsonFile) = print(io, d.path)

"""
    JsonLines(io)

a JSON Lines file opened
"""
struct JsonLines
    io::IO
end

open(d::JsonFile) = JsonLines(open(d.path, "w"))
open(d::JsonFile, mark::Integer) = JsonLines(reopened(d.path, mark))

close(o::JsonLines) = close(o.io)

resume_mark(o::JsonLines) = (flush(o.io); position(o.io))

"""
    json_text(x)

the value `x` in JSON: a vector as an array, a matrix as the array of its rows, a complex
number as `{"re": …, "im": …}`, `missing` as `null`, a float that is not finite as the string
`"Inf"`, `"-Inf"` or `"NaN"`, which JSON does not have, a named tuple as an object, anything
else that JSON has not as the string it prints
"""
json_text(x::AbstractFloat) = isfinite(x) ? JSON.json(x) : JSON.json(string(x))
json_text(x::Complex) = json_object(("re" => real(x), "im" => imag(x)))
json_text(::Missing) = "null"
json_text(x::AbstractVector) = "[" * join(map(json_text, x), ",") * "]"
json_text(x::AbstractMatrix) = json_text([ x[i, :] for i in axes(x, 1) ])
json_text(x::NamedTuple) = json_object(pairs(x))
json_text(x::Union{Real, AbstractString, Nothing}) = JSON.json(x)
json_text(x) = JSON.json(string(x))

"""
    json_object(ps)

the pairs `ps` as a JSON object, in their order
"""
json_object(ps) =
    "{" * join([ JSON.json(string(k)) * ":" * json_text(v) for (k, v) in ps ], ",") * "}"

function write_row(o::JsonLines, time, row)
    println(o.io, json_object([ "time" => time ; row ]))
    flush(o.io)
    return nothing
end

"""
    Data(name)

a destination kept in memory, a vector of rows, each the time then the pairs `name => value`
measured, which `data_to_frame` turns into a table
"""
struct Data <: Destination
    name::String
end

"""
    DataRows(rows)

a `Data` opened, its rows
"""
struct DataRows
    rows::Vector{Vector{Pair{String, Any}}}
end

"""
    RowsMark(rows, n)

the mark of a `Data`, its first `n` rows: the rows themselves, which only grow, rather than a
copy, which a run marking every step would make of all its rows each time
"""
struct RowsMark
    rows::Vector{Vector{Pair{String, Any}}}
    n::Int
end

"""
    kept(m)

the mark `m` as it is kept beyond the run, a `RowsMark` holding its rows alone
"""
kept(m::RowsMark) = RowsMark(m.rows[1:m.n], m.n)
kept(m) = m

open(::Data) = DataRows(Vector{Pair{String, Any}}[])
open(::Data, m::RowsMark) = DataRows(m.rows[1:m.n])
close(::DataRows) = nothing
resume_mark(o::DataRows) = RowsMark(o.rows, length(o.rows))

function write_row(o::DataRows, time, row)
    push!(o.rows, Pair{String, Any}[ "time" => time ; row ])
    return nothing
end

"""
    data_to_frame(rows)

a `DataFrame` of the rows of a `Data`: a column `time` and a column per measurement, in the
order they first come, each value as it was measured, `missing` in a row without it. The
`DataFrames` package must be loaded.

# Examples

    using DataFrames
    df = data_to_frame(sim.data["magnetization"])
"""
function data_to_frame end

"""
    Outputs(; dir = "", format = "%.8g")

where the lists of measurements of a run go: the destinations opened, by their description,
the files named relative to `dir`, the numbers of the text tables in the Printf `format`, the
time as the values
"""
struct Outputs
    dir::String
    format::String
    opened::Dict{Destination, Any}
    lock::ReentrantLock
end

function Outputs(; dir::AbstractString = "", format::AbstractString = "%.8g")
    ok = try
        Printf.format(Printf.Format(format), 1.5)
        true
    catch
        false
    end
    if !ok
        error("the format $(repr(format)) does not write one number, as \"%.8g\" does")
    end
    return Outputs(dir, format, Dict{Destination, Any}(), ReentrantLock())
end

"""
    destination(o, d)

the destination a name stands for in the outputs `o`: `"stdout"` or `"-"` and `"stderr"`
the streams of the process, `""` nowhere, a name ending in `.json` or `.jsonl` a JSON Lines
file, any other a text file, relative to the directory of `o`; a `Destination` as it is
"""
function destination(o::Outputs, name::AbstractString)
    if name == "stdout" || name == "-"
        return Stream(stdout, o.format)
    elseif name == "stderr"
        return Stream(stderr, o.format)
    elseif name == ""
        return Stream(devnull, o.format)
    end
    path = normpath(joinpath(o.dir, name))
    return last(splitext(name)) in (".json", ".jsonl") ? JsonFile(path) : TextFile(path, o.format)
end

destination(::Outputs, d::Destination) = d

destination(::Outputs, x) =
    error("$(repr(x)) is not a destination: a file name, \"stdout\", \"stderr\", \"\", " *
          "Data(name) or a Destination")

"""
    output!(o, state, destination => measurements; time, values...)
    output!(o, state, [destination => measurements, ...]; time, values...)

measure each list on `state`, see `measure`, and write it to its destination, see `destination`:
the destinations are opened and checked first, then everything is measured, then written
"""
output!(o::Outputs, state::AbstractState, p::Pair; kwargs...) = output!(o, state, [p]; kwargs...)

function output!(o::Outputs, state::AbstractState, list::AbstractVector; time = nothing, values...)
    lock(o.lock) do
        items = [ destination(o, first(p)) => as_measure(last(p)) for p in list ]
        for (d, m) in items
            if any(e -> first(e) == "time", m.entries)
                error("a measurement written to $d cannot be named time, the column of the times")
            end
            prepare_list!(get!(() -> open(d), o.opened, d), m)
        end
        rows = [ measure(state, m; time, values...) for (_, m) in items ]
        for ((d, m), row) in zip(items, rows)
            put_row!(o.opened[d], time, m, row)
        end
    end
    return nothing
end

as_measure(m::Measure) = m
as_measure(x) = Measure(x)

"""
    close!(o)

close every destination of the outputs `o`, the rows of a `Data` staying there
"""
function close!(o::Outputs)
    lock(o.lock) do
        foreach(close, values(o.opened))
    end
    return nothing
end

"""
    data_rows(o, name)

the rows the outputs `o` gathered in `Data(name)`
"""
data_rows(o::Outputs, name::String) = o.opened[Data(name)].rows

"""
    resume_marks(o)
    reopen!(o, marks)

the marks of the destinations of the outputs `o` that have one, see `resume_mark`, and the
outputs reopened at them, as a resumed run does
"""
resume_marks(o::Outputs) =
    lock(o.lock) do
        Dict{Destination, Any}(d => m for (d, m) in ((d, resume_mark(h)) for (d, h) in o.opened)
                               if !isnothing(m))
    end

function reopen!(o::Outputs, marks::AbstractDict)
    lock(o.lock) do
        for (d, m) in marks
            o.opened[d] = open(d, m)
        end
    end
    return o
end
