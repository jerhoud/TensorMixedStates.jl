# The outputs of the measurements: Destination, the open type of where they are written, which
# other packages extend, the text tables, the JSON files and the Data of TMS, and Outputs,
# which measures lists on a state and writes each to its destination.

export Data, data_to_frame

"""
    Destination

where measurements are written: a description, compared by value, opened once and given a row
at each call. TMS writes text tables, to files and to the streams of the process, JSON files,
all named by a string, and `Data`. A package adds its own by a subtype and these methods
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
    cells!(out, col, x)

`out`, the cells of the value `x` under the column `col` pushed to it as pairs
`column => value`: a complex number as `Re(col)` and `Im(col)`, a vector as `col[i]`, a matrix
row by row as `col[i,j]`, a named tuple as `col.field`, `nothing` as `missing`, anything else
as it is, which a text table writes by `text_cell`
"""
function cells!(out, col, x::Complex)
    cells!(out, "Re($col)", real(x))
    return cells!(out, "Im($col)", imag(x))
end

cells!(out, col, ::Nothing) = push!(out, col => missing)

function cells!(out, col, x::AbstractVector)
    for (i, y) in enumerate(x)
        cells!(out, "$col[$i]", y)
    end
    return out
end

function cells!(out, col, x::AbstractMatrix)
    for i in axes(x, 1), j in axes(x, 2)
        cells!(out, "$col[$i,$j]", x[i, j])
    end
    return out
end

function cells!(out, col, x::NamedTuple)
    for (k, y) in pairs(x)
        cells!(out, "$col.$k", y)
    end
    return out
end

cells!(out, col, x) = push!(out, col => x)

"""
    text_cell(x, format)

the text of the cell `x` in a text table, see `cells!`: a real number in `format`, an integer
as it is, `missing` as an empty cell, anything else as it prints, on one line
"""
text_cell(x::Integer, _) = string(x)
text_cell(x::Real, f) = Printf.format(f, float(x))
text_cell(::Missing, _) = ""
text_cell(x, _) = replace(string(x), r"[\t\r\n]" => " ")

"""
    entry_cells(name, m, x)

the cells of the value `x` of the measurement `m` named `name`, see `cells!`: an operator of one
site on site `i` under `name(i)`, a pair on sites `i` and `j` under `A(i)B(j)`, or `name(i,j)`
when it is named, and a measurement of the state named `φ` under the cells it has on its own,
prefixed by `φ.`. The cells are those `cells!` gives of `x`, in the same order, under other
names.
"""
function entry_cells(name::String, m, x)
    out = Pair{String, Any}[]
    if m isa OnNamed
        prefix = "$(m.state)."
        inner = entry_cells(String(chopprefix(name, prefix)), m.entry, x)
        return Pair{String, Any}[ prefix * c => v for (c, v) in inner ]
    elseif m isa SimpleOp && x isa AbstractVector
        for (i, y) in enumerate(x)
            cells!(out, "$name($i)", y)
        end
    elseif m isa Tuple{SimpleOp, SimpleOp} && x isa AbstractMatrix
        a, b = map(obs_name, m)
        col(i, j) = name == entry_name(m) ? "$a($i)$b($j)" : "$name($i,$j)"
        for i in axes(x, 1), j in axes(x, 2)
            cells!(out, col(i, j), x[i, j])
        end
    else
        cells!(out, name, x)
    end
    return out
end

function put_row!(t::TextTable, time, m::Measure, row)
    times = cells!(Pair{String, Any}[], "time", time)
    parts = [ [times] ;
              [ entry_cells(name, e, x) for ((name, e), (_, x)) in zip(m.entries, row) ] ]
    columns = [ first.(p) for p in parts ]
    if columns ≠ t.columns
        if !isnothing(t.path) && !isnothing(t.columns)
            k = findfirst(k -> columns[k] ≠ t.columns[k], eachindex(columns))
            what = k == 1 ? "the time does" : "the values of $(first(row[k - 1])) do"
            error("$what not fit the columns $(t.path) took at its first line: a value changing " *
                  "its cells, as one missing at first, needs a .json file")
        end
        names = reduce(vcat, columns)
        twice = findfirst(c -> count(==(c), names) > 1, names)
        if !isnothing(twice)
            error("two columns of $(something(t.path, "the table")) would be named " *
                  "$(names[twice]): give one of the measurements a name of its own, as " *
                  "\"name\" => op")
        end
        t.columns = columns
        println(t.io, join(names, '\t'))
    end
    println(t.io, join([ text_cell(last(c), t.format) for p in parts for c in p ], '\t'))
    flush(t.io)
    return nothing
end

"""
    JsonFile(path)

a JSON file of one object per row: a `.jsonl` file JSON Lines, an object per line, and a
`.json` file the array of the objects, an object per line too, which closing the file ends
"""
struct JsonFile <: Destination
    path::String
end

show(io::IO, d::JsonFile) = print(io, d.path)

"""
    is_array(d)

whether the JSON file `d` is an array, a `.json` file, rather than JSON Lines
"""
is_array(d::JsonFile) = last(splitext(d.path)) == ".json"

"""
    JsonOpened(io, array, empty)

a JSON file opened, an array or JSON Lines, and for an array whether it holds no object yet
"""
mutable struct JsonOpened
    io::IO
    array::Bool
    empty::Bool
end

function open(d::JsonFile)
    io = open(d.path, "w")
    if is_array(d)
        print(io, "[")
    end
    return JsonOpened(io, is_array(d), true)
end

# an array cut at a mark lost its closing bracket, and holds an object past its opening one
function open(d::JsonFile, mark::Integer)
    io = reopened(d.path, mark)
    if is_array(d) && position(io) == 0
        print(io, "[")
    end
    return JsonOpened(io, is_array(d), position(io) ≤ 1)
end

function close(o::JsonOpened)
    if o.array
        print(o.io, "\n]\n")
    end
    close(o.io)
end

resume_mark(o::JsonOpened) = (flush(o.io); position(o.io))

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

function write_row(o::JsonOpened, time, row)
    text = json_object([ "time" => time ; row ])
    if o.array
        print(o.io, o.empty ? "\n" : ",\n", text)
        o.empty = false
    else
        println(o.io, text)
    end
    flush(o.io)
    return nothing
end

"""
    Data(name)

a destination kept in memory, a vector of rows, each the time then the pairs `name => value`
measured, see `DataRow`, and a table of the cells of these values, see `DataTable`
"""
struct Data <: Destination
    name::String
end

"""
    DataRow

a row of a `Data`: the pairs `name => value` measured, the time first, each value as it was
measured, `row[name]` giving the value of `name`, and the columns of the cells of these values,
see `cells!`
"""
struct DataRow <: AbstractVector{Pair{String, Any}}
    measured::Vector{Pair{String, Any}}
    # the column of each cell, the time first, shared by the rows that have the same
    columns::Vector{String}
end

Base.size(r::DataRow) = size(r.measured)
Base.IndexStyle(::Type{DataRow}) = IndexLinear()
Base.getindex(r::DataRow, k::Int) = r.measured[k]

function Base.getindex(r::DataRow, name::AbstractString)
    k = findfirst(p -> first(p) == name, r.measured)
    if isnothing(k)
        throw(KeyError(name))
    end
    return last(r.measured[k])
end

"""
    measured_cells(measured)

the cells of the values of the pairs `measured`, see `cells!`, in their order, without their
columns
"""
measured_cells(measured) =
    Any[ last(c) for (name, x) in measured for c in cells!(Pair{String, Any}[], name, x) ]

# a row of the pairs `measured`, the columns of its cells named as `cells!` names those of
# values of no measurement of their own
DataRow(measured::AbstractVector) =
    DataRow(collect(Pair{String, Any}, measured),
            [ first(c) for (name, x) in measured for c in cells!(Pair{String, Any}[], name, x) ])

"""
    DataTable

the rows of a `Data`, see `DataRow`, as `sim.data[name]` gives them: a vector of rows, and a
table of `Tables.jl` whose columns are those a text table of the same measurements has, see
`cells!`, a cell that a row lacks being `missing`, so that `DataFrame(rows)` and
`CSV.write(file, rows)` take it
"""
struct DataTable <: AbstractVector{DataRow}
    rows::Vector{DataRow}
end

Base.size(t::DataTable) = size(t.rows)
Base.IndexStyle(::Type{DataTable}) = IndexLinear()
Base.getindex(t::DataTable, k::Int) = t.rows[k]

"""
    DataColumns

the columns of a `DataTable`, by `Tables.columns`
"""
struct DataColumns <: Tables.AbstractColumns
    names::Vector{Symbol}
    columns::Vector{AbstractVector}
    index::Dict{Symbol, Int}
end

Tables.columnnames(c::DataColumns) = getfield(c, :names)
Tables.getcolumn(c::DataColumns, k::Int) = getfield(c, :columns)[k]
Tables.getcolumn(c::DataColumns, name::Symbol) = getfield(c, :columns)[getfield(c, :index)[name]]

"""
    data_columns(rows)

the columns of the cells of `rows`, in the order they first come, a column only `time` when
there is no row, refused when two cells of a row share a column, as a text table refuses them
"""
function data_columns(rows::Vector{DataRow})
    names = String["time"]
    columns = Vector{Any}[Any[]]
    index = Dict("time" => 1)
    for (r, row) in enumerate(rows)
        cells = measured_cells(row.measured)
        if length(cells) ≠ length(row.columns)
            error("bug: a row of a Data has $(length(cells)) cells and $(length(row.columns)) " *
                  "columns")
        end
        for (name, x) in zip(row.columns, cells)
            k = get!(index, name) do
                push!(names, name)
                push!(columns, Any[ missing for _ in 1:r - 1 ])
                return length(columns)
            end
            if length(columns[k]) ≥ r
                error("two columns of the table would be named $name: give one of the " *
                      "measurements a name of its own, as \"name\" => op")
            end
            push!(columns[k], x)
        end
        for c in columns
            if length(c) < r
                push!(c, missing)
            end
        end
    end
    return DataColumns(Symbol.(names), AbstractVector[ map(identity, c) for c in columns ],
                       Dict(Symbol(n) => k for (k, n) in enumerate(names)))
end

Tables.istable(::Type{DataTable}) = true
Tables.columnaccess(::Type{DataTable}) = true
Tables.columns(t::DataTable) = data_columns(t.rows)

"""
    DataRows(rows)

a `Data` opened, its rows
"""
struct DataRows
    rows::Vector{DataRow}
end

"""
    RowsMark(rows, n)

the mark of a `Data`, its first `n` rows: the rows themselves, which only grow, rather than a
copy, which a run marking every step would make of all its rows each time
"""
struct RowsMark
    rows::Vector{DataRow}
    n::Int
end

"""
    kept(m)

the mark `m` as it is kept beyond the run, a `RowsMark` holding its rows alone
"""
kept(m::RowsMark) = RowsMark(m.rows[1:m.n], m.n)
kept(m) = m

open(::Data) = DataRows(DataRow[])
open(::Data, m::RowsMark) = DataRows(m.rows[1:m.n])
close(::DataRows) = nothing
resume_mark(o::DataRows) = RowsMark(o.rows, length(o.rows))

function write_row(o::DataRows, time, row)
    push!(o.rows, DataRow(Pair{String, Any}[ "time" => time ; row ]))
    return nothing
end

# the columns of the cells named as a text table of the list `m` names them, see `entry_cells`,
# those of the row before when they are the same
function put_row!(o::DataRows, time, m::Measure, row)
    measured = Pair{String, Any}[ "time" => time ; row ]
    cells = [ c for ((name, e), (_, x)) in zip(m.entries, row) for c in entry_cells(name, e, x) ]
    columns = String[ "time" ; first.(cells) ]
    if !isempty(o.rows) && last(o.rows).columns == columns
        columns = last(o.rows).columns
    end
    push!(o.rows, DataRow(measured, columns))
    return nothing
end

"""
    data_to_frame(rows)

the `DataFrame` of the rows of a `Data`, `DataFrame(rows)`, see `DataTable`: the column `time`,
then the columns a text table of the same measurements has, in the order they first come,
`missing` in a row without them. The `DataFrames` package must be loaded.

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
the streams of the process, `""` nowhere, a name ending in `.json` or `.jsonl` a JSON file,
see `JsonFile`, any other a text file, relative to the directory of `o`; a `Destination` as it
is
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

as_measure(m::Measure) = m
as_measure(x) = Measure(x)

"""
    by_destination(o, list)

the pairs `destination => measurements` of `list` as pairs `destination => Measure`, one for
each destination, in the order they first come, the lists going to one destination gathered
into one, so that a call writes one row to each
"""
function by_destination(o::Outputs, list::AbstractVector)
    ps = [ destination(o, first(p)) => last(p) for p in list ]
    ds = unique(first.(ps))
    ms = [ [ last(p) for p in ps if first(p) == d ] for d in ds ]
    return [ d => (length(m) == 1 ? as_measure(only(m)) : Measure(m)) for (d, m) in zip(ds, ms) ]
end

"""
    output!(o, state, destination => measurements; time, values...)
    output!(o, state, [destination => measurements, ...]; time, values...)

measure each list on `state`, see `measure`, and write it to its destination, see `destination`,
one row to each destination, the lists going to one destination measured as one: the
destinations are opened and checked first, then everything is measured, then written
"""
output!(o::Outputs, state::AbstractState, p::Pair; kwargs...) = output!(o, state, [p]; kwargs...)

function output!(o::Outputs, state::AbstractState, list::AbstractVector; time = nothing, values...)
    lock(o.lock) do
        items = by_destination(o, list)
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
data_rows(o::Outputs, name::String) = DataTable(o.opened[Data(name)].rows)

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
