# The outputs of the measurements: Destination, the open type of where they are written, which
# other packages extend, the text tables, the JSON files and the Data of TMS, the averages of the
# samples of a simulation and EachSample, and Outputs, which measures lists on a state and writes
# each to its destination.

export Data, EachSample, data_to_frame

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
    MeasuredRow(time, values, parts)

a row of measurements as the destinations take it: the time, the pairs `name => value`
measured, and the cells of each value, see `entry_cells`, which hold no operator, so that the
row is written again in another process, see `measured_row`
"""
struct MeasuredRow
    time::Any
    values::Vector{Pair{String, Any}}
    parts::Vector{Vector{Pair{String, Any}}}
end

# the cells of a row in their order, as a text table names its columns
row_cells(r::MeasuredRow) = Pair{String, Any}[ c for p in r.parts for c in p ]

"""
    prepare_list!(o, names)

check that the destination opened `o` can take the list of measurements named `names`, see
`check_list`, a text table keeping what it needs to name its columns
"""
prepare_list!(o, names::Vector{String}) = check_list(o, names)

"""
    put_row!(o, r)

write the `MeasuredRow` `r` to the destination opened `o`, see `write_row`
"""
put_row!(o, r::MeasuredRow) = write_row(o, r.time, r.values)

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
function prepare_list!(t::TextTable, names::Vector{String})
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

"""
    measured_row(m, time, row)

the `MeasuredRow` of the values `row` of the `Measure` `m` measured at `time`
"""
measured_row(m::Measure, time, row) =
    MeasuredRow(time, Pair{String, Any}[ p for p in row ],
                [ entry_cells(name, e, x) for ((name, e), (_, x)) in zip(m.entries, row) ])

function put_row!(t::TextTable, r::MeasuredRow)
    row = r.values
    parts = [ [cells!(Pair{String, Any}[], "time", r.time)] ; r.parts ]
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
measured, see `DataRow`, and a table of the cells of these values, see `DataTable`; in a
simulation with samples, the cells averaged, `N`, `N_eff`, `c`, `err(c)` and `std(c)`, see
`Averaged`
"""
struct Data <: Destination
    name::String
end

"""
    DataRow

a row of a `Data`: the pairs `name => value` measured, the time first, each value as it was
measured, `row[name]` giving the value of `name`, the cell `name` for a row averaged, and the
columns of the cells of these values, see `cells!`
"""
struct DataRow <: AbstractVector{Pair{String, Any}}
    measured::Vector{Pair{String, Any}}
    # the column of each cell, the time first, shared by the rows that have the same
    columns::Vector{String}
    # typed: the default constructor, of any two arguments, is ambiguous with those the oldest
    # versions of ITensors and NDTensors allowed define for every type of array
    DataRow(measured::Vector{Pair{String, Any}}, columns::Vector{String}) = new(measured, columns)
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

"""
    DataTable

the rows of a `Data`, see `DataRow`, as `sim.data[name]` gives them: a vector of rows, and a
table of `Tables.jl` whose columns are those a text table of the same measurements has, see
`cells!`, a cell that a row lacks being `missing`, so that `DataFrame(rows)` and
`CSV.write(file, rows)` take it; a slice or a filter of the rows is one as well
"""
struct DataTable <: AbstractVector{DataRow}
    rows::Vector{DataRow}
    # typed, as for DataRow
    DataTable(rows::Vector{DataRow}) = new(rows)
end

Base.size(t::DataTable) = size(t.rows)
Base.IndexStyle(::Type{DataTable}) = IndexLinear()
Base.getindex(t::DataTable, k::Int) = t.rows[k]
Base.setindex!(t::DataTable, r::DataRow, k::Int) = setindex!(t.rows, r, k)
# so that a slice or a filter of the rows is a table as well
Base.similar(::DataTable, ::Type{DataRow}, dims::Dims{1}) =
    DataTable(Vector{DataRow}(undef, dims...))

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

# the columns of the cells named as a text table of the list `m` names them, see `entry_cells`,
# those of the row before when they are the same
function put_row!(o::DataRows, r::MeasuredRow)
    measured = Pair{String, Any}[ "time" => r.time ; r.values ]
    cells = row_cells(r)
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
    EachSample(d)

the destination `d`, a file name or a destination, written a row for each sample in a simulation
with samples, after the time the column `sample` giving its index, 0 before the samples, rather
than the averages over the samples its other destinations hold, see `SimData`, a trajectory
writing the logarithm of its weight in the column `log_weight`; refused in a simulation without
samples

# Examples

    final_measurements = [ "energy" => :energy, EachSample("energies") => :energy ]
"""
struct EachSample <: Destination
    inner::Union{AbstractString, Destination}
end

show(io::IO, d::EachSample) = print(io, "EachSample(", d.inner, ")")

# the columns the rows of the samples write, or will: no measurement written there takes their name
const sample_columns = ("sample", "weight", "trajectory", "log_weight")

# opened as the destination it stands for, once a name is made a destination, see by_destination
open(d::EachSample) = open(d.inner)
open(d::EachSample, mark) = open(d.inner, mark)

"""
    Averaged(d)

the destination `d`, a text or JSON file or a `Data`, of a simulation with samples: for each row,
by its phase and step, the average of the samples, see `AveragedOpened`
"""
struct Averaged <: Destination
    inner::Union{TextFile, JsonFile, Data}
end

show(io::IO, d::Averaged) = show(io, d.inner)

"""
    RowKey

`(phase, step, n)`, the `n`-th row of a step in a sample; `typemax(Int)` as the step for the final
measurements of a phase, as the phase for those of the simulation
"""
const RowKey = NTuple{3, Int}

# the place of the row `key` in the messages
row_place((phase, step, _)::RowKey) =
    phase == typemax(Int) ? "the final measurements" :
    step == typemax(Int) ? "the end of phase $phase" : "step $step of phase $phase"

"""
    RowSum(time, cells)

the sums over the samples of a row of an averaged destination, by the algorithm of West (1979),
each cell shifted by its value in the first sample (Chan, Golub and LeVeque 1983)
"""
struct RowSum
    time::Any
    columns::Vector{String}
    missing::Vector{Bool}
    shift::Vector{Float64}
    n::Int
    W::Float64
    W2::Float64
    mean::Vector{Float64}
    M2::Vector{Float64}
end

function RowSum(time, cells::Vector{Pair{String, Any}})
    k = length(cells)
    return RowSum(time, first.(cells), [ ismissing(x) for (_, x) in cells ],
                  [ ismissing(x) ? 0. : Float64(x) for (_, x) in cells ], 0, 0., 0., zeros(k),
                  zeros(k))
end

# whether two times of a row are the same but for rounding, see `rounding_tol`: a step computed
# from another duration differs by an ulp
same_time(a, b) =
    isequal(a, b) ||
    (a isa Number && b isa Number && abs(a - b) ≤ rounding_tol * max(abs(a), abs(b)))

"""
    added(s, cells, time, key, d)

the sums `s` of the row `key` of the averaged destination `d` with a sample, refused unless it has
the same time, columns and missing cells; a trajectory weighed by the exponential of its cell
`log_weight`, which is shifted by the largest so far, the sums rescaled when it grows, so that no
weight overflows
"""
function added(s::RowSum, cells::Vector{Pair{String, Any}}, time, key::RowKey, d)
    if !same_time(time, s.time)
        error("the samples of $d have the times $(s.time) and $time at $(row_place(key)): " *
              "averaged rows need the same times, which EachSample(…) does not")
    elseif first.(cells) ≠ s.columns
        error("the samples of $d have other columns at $(row_place(key)): averaged rows need the " *
              "same, which EachSample(…) does not")
    end
    shift, mean, M2 = copy(s.shift), copy(s.mean), copy(s.M2)
    W, W2 = s.W, s.W2
    iw = findfirst(==("log_weight"), s.columns)
    w = 1.
    if !isnothing(iw)
        lw = Float64(last(cells[iw]))
        if lw > shift[iw]
            f = exp(shift[iw] - lw)
            W, W2, M2 = f * W, f^2 * W2, f * M2
            shift[iw] = lw
        end
        w = exp(lw - shift[iw])
    end
    W += w
    for (i, (column, x)) in enumerate(cells)
        if i == iw
            continue
        elseif ismissing(x) ≠ s.missing[i]
            error("$column of $d is missing in some samples only at $(row_place(key)): it cannot " *
                  "be averaged, which EachSample(…) does not need")
        elseif !ismissing(x)
            # exact where the deviations are small against the value, which keeps their digits
            y = x - shift[i]
            delta = y - mean[i]
            mean[i] += (w / W) * delta
            M2[i] += w * delta * (y - mean[i])
        end
    end
    return RowSum(s.time, s.columns, s.missing, shift, s.n + 1, W, W2 + w^2, mean, M2)
end

# the standard deviation of the samples in the cell `i` of the row of sums `s` and the error on
# their mean, NaN for one sample
function spread(s::RowSum, i::Int)
    unbiased = s.W - s.W2 / s.W
    sd = unbiased > 0 ? sqrt(s.M2[i] / unbiased) : NaN
    return sd, sd / sqrt(s.W^2 / s.W2)
end

"""
    averaged_cells(s, samples = s.n)

the cells of the row of sums `s`: `N`, `N_eff`, and for each cell its mean, the error on it and
the standard deviation of the samples, `NaN` for one sample; the mean weight over the `samples`
merged, those a weight of zero stopped before the row included
"""
function averaged_cells(s::RowSum, samples::Int = s.n)
    out = Pair{String, Any}[ "N" => s.n, "N_eff" => s.W^2 / s.W2 ]
    for (i, c) in enumerate(s.columns)
        if c == "log_weight"
            # the mean weight, its error relative to it from the spread of the weights
            n = max(samples, s.n)
            m = s.W / n
            err = n > 1 ? sqrt(max(s.W2 / n - m^2, 0.) / (n - 1)) / m : NaN
            append!(out, [ c => s.shift[i] + log(m), "err($c)" => err ])
        elseif s.missing[i]
            append!(out, [ c => missing, "err($c)" => missing, "std($c)" => missing ])
        else
            sd, err = spread(s, i)
            append!(out, [ c => s.shift[i] + s.mean[i], "err($c)" => err, "std($c)" => sd ])
        end
    end
    return out
end

# the cells of a row written before the samples, which all share it: no error, no deviation
function fixed_cells(cells::Vector{Pair{String, Any}}, samples::Int)
    out = Pair{String, Any}[ "N" => samples, "N_eff" => Float64(samples) ]
    for (c, x) in cells
        z = ismissing(x) ? missing : 0.
        append!(out, [ c => x, "err($c)" => z, "std($c)" => z ])
    end
    return out
end

"""
    AveragedOpened

an averaged destination opened: the sums of the samples ended, the rows before the samples, the
samples ended, the rows of the sample being run and the list of a text file
"""
mutable struct AveragedOpened
    destination::Averaged
    sums::Dict{RowKey, RowSum}
    fixed::Dict{RowKey, Tuple{Any, Vector{Pair{String, Any}}}}
    samples::Int
    rows::Vector{Tuple{RowKey, Any, Vector{Pair{String, Any}}}}
    names::Union{Nothing, Vector{String}}
end

"""
    AveragedMark

the mark of an averaged destination, its first `n` rows of the sample being run as `RowsMark`;
the sums and the rows before the samples, replaced rather than changed, are held without a copy
"""
struct AveragedMark
    sums::Dict{RowKey, RowSum}
    fixed::Dict{RowKey, Tuple{Any, Vector{Pair{String, Any}}}}
    samples::Int
    rows::Vector{Tuple{RowKey, Any, Vector{Pair{String, Any}}}}
    n::Int
    names::Union{Nothing, Vector{String}}
end

kept(m::AveragedMark) = AveragedMark(m.sums, m.fixed, m.samples, m.rows[1:m.n], m.n, m.names)

open(d::Averaged) =
    AveragedOpened(d, Dict{RowKey, RowSum}(), Dict{RowKey, Tuple{Any, Vector{Pair{String, Any}}}}(),
                   0, Tuple{RowKey, Any, Vector{Pair{String, Any}}}[], nothing)

open(d::Averaged, m::AveragedMark) =
    AveragedOpened(d, m.sums, m.fixed, m.samples, m.rows[1:m.n], m.names)

resume_mark(o::AveragedOpened) =
    AveragedMark(o.sums, o.fixed, o.samples, o.rows, length(o.rows), o.names)

function prepare_list!(o::AveragedOpened, names::Vector{String})
    inner = o.destination.inner
    if !(inner isa TextFile) || o.names == names
        return nothing
    elseif !isnothing(o.names)
        error("$(inner.path) holds the list $(join(o.names, ", ")), and another needs another file")
    end
    o.names = names
    return nothing
end

"""
    put_keyed_row!(o, at, r)

`put_row!(o, r)`, the row being written at `at`, `(sample, phase, step)`, which an averaged
destination keeps it by, see `RowKey`
"""
put_keyed_row!(o, _, r::MeasuredRow) = put_row!(o, r)

"""
    check_averaged(d, cells, phase, step)

refuse the cells of a row of the averaged destination `d` at step `step` of phase `phase` that
cannot be averaged: two columns of one name, a column named as the averages or the samples
name theirs, see `sample_columns`, a value other than a number
"""
function check_averaged(d::Averaged, cells::Vector{Pair{String, Any}}, phase::Int, step::Int)
    columns = first.(cells)
    twice = findfirst(c -> count(==(c), columns) > 1, columns)
    taken = findfirst(c -> c in ("N", "N_eff") || c in sample_columns || startswith(c, "err(") ||
                           startswith(c, "std("), columns)
    if !isnothing(twice)
        error("two columns of $d would be named $(columns[twice]): give one of the measurements " *
              "a name of its own, as \"name\" => op")
    elseif !isnothing(taken)
        error("the column $(columns[taken]) of $d is one the averages or the samples write: " *
              "give the measurement another name, as \"name\" => op")
    end
    bad = findfirst(c -> !(last(c) isa Union{Real, Missing}), cells)
    if !isnothing(bad)
        error("$(first(cells[bad])) of $d is $(repr(last(cells[bad]))) at " *
              "$(row_place((phase, step, 0))), which cannot be averaged: write it with " *
              "EachSample(…)")
    end
    return nothing
end

function put_keyed_row!(o::AveragedOpened, (sample, phase, step), r::MeasuredRow)
    cells = row_cells(r)
    weighted = !isempty(cells) && first(last(cells)) == "log_weight"
    check_averaged(o.destination, weighted ? cells[1:end - 1] : cells, phase, step)
    if sample == 0
        n = 1 + count(k -> k[1] == phase && k[2] == step, keys(o.fixed))
    else
        # the rows of a sample come in their order: those of the same step are the last ones
        n = 1
        while n ≤ length(o.rows) && o.rows[end - n + 1][1][1:2] == (phase, step)
            n += 1
        end
    end
    if sample == 0
        o.fixed = merge(o.fixed, Dict((phase, step, n) => (r.time, cells)))
    else
        push!(o.rows, ((phase, step, n), r.time, cells))
    end
    return nothing
end

"""
    fold!(o)

the rows of the sample being run added to the sums of the averaged destination `o`
"""
function fold!(o::AveragedOpened)
    sums = copy(o.sums)
    for (key, time, cells) in o.rows
        s = get(() -> RowSum(time, cells), sums, key)
        sums[key] = added(s, cells, time, key, o.destination)
    end
    o.sums = sums
    o.samples += 1
    o.rows = Tuple{RowKey, Any, Vector{Pair{String, Any}}}[]
    return o
end

# the rows of an averaged destination, `time => cells`, in the order of their keys
function averaged_rows(sums::Dict{RowKey, RowSum}, fixed::AbstractDict, samples::Int)
    ks = sort!(collect(union(keys(sums), keys(fixed))))
    return [ haskey(fixed, k) ? first(fixed[k]) => fixed_cells(last(fixed[k]), samples) :
                                sums[k].time => averaged_cells(sums[k], samples) for k in ks ]
end

"""
    write_averaged(d, sums, fixed, samples)

write the averaged destination `d` whole to a file renamed into place, so that it is never seen
half written; nothing for a `Data`
"""
function write_averaged(d::Averaged, sums, fixed, samples::Int)
    inner = d.inner
    if inner isa Data
        return nothing
    end
    rows = averaged_rows(sums, fixed, samples)
    tmp = inner.path * ".tmp"
    if inner isa TextFile
        open(tmp, "w") do io
            format = Printf.Format(inner.format)
            header = nothing
            for (time, cells) in rows
                line = [ cells!(Pair{String, Any}[], "time", time) ; cells ]
                # rows of other columns, written by phases of other lists
                if first.(line) ≠ header
                    header = first.(line)
                    println(io, join(header, '\t'))
                end
                println(io, join([ text_cell(last(c), format) for c in line ], '\t'))
            end
        end
    else
        array = is_array(inner)
        io = open(tmp, "w")
        print(io, array ? "[" : "")
        o = JsonOpened(io, array, true)
        for (time, cells) in rows
            write_row(o, time, cells)
        end
        close(o)
    end
    Base.Filesystem.rename(tmp, inner.path)
    return nothing
end

close(o::AveragedOpened) = write_averaged(o.destination, o.sums, o.fixed, o.samples)

# the rows of an averaged `Data`, as `sim.data` gives them
function averaged_data(o::AveragedOpened)
    rows = averaged_rows(o.sums, o.fixed, o.samples)
    return DataTable(DataRow[ DataRow(Pair{String, Any}[ "time" => t ; cells ],
                                      String[ "time" ; first.(cells) ]) for (t, cells) in rows ])
end

"""
    Recorder(destination)

a file or `Data` of a simulation with samples as a sample opens it: its rows kept, each with
where it is written and the names of its list, for `replay!`; checked as `check_averaged` does
"""
mutable struct Recorder
    destination::Destination
    names::Union{Nothing, Vector{String}}
    list::Vector{String}
    rows::Vector{Tuple{NTuple{3, Int}, Vector{String}, MeasuredRow}}
end

Recorder(d::Destination) =
    Recorder(d, nothing, String[], Tuple{NTuple{3, Int}, Vector{String}, MeasuredRow}[])

# whether a part of a simulation with samples records what goes to `d`, see `Recorder`
recorded(d::Destination) =
    d isa Averaged || (d isa EachSample && d.inner isa Union{TextFile, JsonFile, Data})

function prepare_list!(r::Recorder, names::Vector{String})
    inner = r.destination.inner
    if inner isa TextFile
        if !isnothing(r.names) && r.names ≠ names
            error("$(inner.path) holds the list $(join(r.names, ", ")), and another needs " *
                  "another file")
        end
        r.names = names
    end
    r.list = names
    return nothing
end

function put_keyed_row!(r::Recorder, at::NTuple{3, Int}, row::MeasuredRow)
    if r.destination isa Averaged
        cells = row_cells(row)
        weighted = !isempty(cells) && first(last(cells)) == "log_weight"
        check_averaged(r.destination, weighted ? cells[1:end - 1] : cells, at[2], at[3])
    end
    push!(r.rows, (at, r.list, row))
    return nothing
end

close(::Recorder) = nothing

"""
    RecordMark(names, rows, n)

the mark of a `Recorder`, its first `n` rows, held without a copy as `RowsMark` holds them
"""
struct RecordMark
    names::Union{Nothing, Vector{String}}
    rows::Vector{Tuple{NTuple{3, Int}, Vector{String}, MeasuredRow}}
    n::Int
end

resume_mark(r::Recorder) = RecordMark(r.names, r.rows, length(r.rows))
kept(m::RecordMark) = RecordMark(m.names, m.rows[1:m.n], m.n)
open(d::Averaged, m::RecordMark) = Recorder(d, m.names, String[], m.rows[1:m.n])
open(d::EachSample, m::RecordMark) = Recorder(d, m.names, String[], m.rows[1:m.n])

"""
    Outputs(; dir = "", format = "%.13g", sampled = false, recording = false)

where the lists of measurements of a run go: the destinations opened, by their description,
the files named relative to `dir`, the numbers of the text tables in the Printf `format`, the
time as the values; for a simulation `sampled`, its files and `Data` averaged over the samples,
see `Averaged`, and where the rows are written, see `put_keyed_row!`; `recording` for a part of
a simulation with samples, which records its rows instead of writing them, see `Recorder`
"""
struct Outputs
    dir::String
    format::String
    opened::Dict{Destination, Any}
    lock::ReentrantLock
    sampled::Bool
    recording::Bool
    at::Base.RefValue{NTuple{3, Int}}
end

function Outputs(; dir::AbstractString = "", format::AbstractString = "%.13g",
                 sampled::Bool = false, recording::Bool = false)
    ok = try
        Printf.format(Printf.Format(format), 1.5)
        true
    catch
        false
    end
    if !ok
        error("the format $(repr(format)) does not write one number, as \"%.13g\" does")
    end
    return Outputs(dir, format, Dict{Destination, Any}(), ReentrantLock(), sampled, recording,
                   Ref((0, 0, 0)))
end

"""
    destination(o, d)

the destination a name stands for in the outputs `o`: `"stdout"` or `"-"` and `"stderr"`
the streams of the process, `""` nowhere, a name ending in `.json` or `.jsonl` a JSON file,
see `JsonFile`, any other a text file, relative to the directory of `o`; a `Destination` as it
is; a file or a `Data` averaged when `o` is sampled, see `Averaged`
"""
function destination(o::Outputs, d)
    x = plain_destination(o, d)
    return o.sampled && x isa Union{TextFile, JsonFile, Data} ? Averaged(x) : x
end

function plain_destination(o::Outputs, name::AbstractString)
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

plain_destination(::Outputs, d::Destination) = d

plain_destination(::Outputs, x) =
    error("$(repr(x)) is not a destination: a file name, \"stdout\", \"stderr\", \"\", " *
          "Data(name) or a Destination")

as_measure(m::Measure) = m
as_measure(x) = Measure(x)

"""
    each_sample(d, m)

the measurements `m` of the destination `EachSample(d)`, the lists going there gathered,
measuring `:sample` first unless they start with it; another of the `sample_columns` refused
"""
function each_sample(d::EachSample, m::Measure)
    sampled = !isempty(m.entries) && first(m.entries) == ("sample" => :sample)
    skip = sampled ? 1 : 0
    taken = findfirst(e -> first(e) in sample_columns, m.entries[skip + 1:end])
    if !isnothing(taken)
        e = m.entries[taken + skip]
        error("the column $(first(e)) of $d is one the samples write: " *
              (e == ("sample" => :sample) ? "put :sample first, or leave it out" :
               "give the measurement another name, as \"name\" => op"))
    end
    return sampled ? m : Measure([:sample, m])
end

"""
    by_destination(o, list)

the pairs `destination => measurements` of `list` as pairs `destination => Measure`, one for
each destination, in the order they first come, the lists going to one destination gathered
into one, so that a call writes one row to each; `EachSample(d)` is opened as `d`, not averaged,
measuring `:sample` first unless the list starts with it, so that a list made this way is made
again as it is, as `output` does with the lists of `phase_output`
"""
function by_destination(o::Outputs, list::AbstractVector)
    ps = [ first(p) isa EachSample ? EachSample(plain_destination(o, first(p).inner)) => last(p) :
           destination(o, first(p)) => last(p) for p in list ]
    ds = unique(first.(ps))
    ms = [ [ last(p) for p in ps if first(p) == d ] for d in ds ]
    gathered(m) = length(m) == 1 ? as_measure(only(m)) : Measure(m)
    return [ d => (d isa EachSample ? each_sample(d, gathered(m)) : gathered(m))
             for (d, m) in zip(ds, ms) ]
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
            elseif d isa EachSample && !o.sampled
                error("EachSample writes a row for each sample, and the simulation has no samples")
            end
            twin = d isa EachSample && d.inner isa Union{TextFile, JsonFile, Data} ?
                   Averaged(d.inner) : d isa Averaged ? EachSample(d.inner) : nothing
            if !haskey(o.opened, d) && !isnothing(twin) && haskey(o.opened, twin)
                error("$(d.inner) is written both averaged and for each sample: give " *
                      "EachSample(…) another destination")
            end
            opened = get!(() -> o.recording && recorded(d) ? Recorder(d) : open(d), o.opened, d)
            prepare_list!(opened, first.(m.entries))
        end
        rows = [ measure(state, m; time, values...) for (_, m) in items ]
        weight = get(values, :log_weight, nothing)
        for ((d, m), row) in zip(items, rows)
            r = measured_row(m, time, row)
            if !isnothing(weight) && d isa Union{Averaged, EachSample}
                r = weighted_row(r, weight)
            end
            put_keyed_row!(o.opened[d], o.at[], r)
        end
    end
    return nothing
end

# the row `r` of a trajectory with the logarithm of its weight, `log_weight`, which the averages
# weigh it by and the rows of each sample write
weighted_row(r::MeasuredRow, log_weight::Real) =
    MeasuredRow(r.time, Pair{String, Any}[ r.values ; "log_weight" => Float64(log_weight) ],
                [ r.parts ; [Pair{String, Any}["log_weight" => Float64(log_weight)]] ])

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

the rows the outputs `o` gathered in `Data(name)`, its averages when it is averaged
"""
function data_rows(o::Outputs, name::String)
    d = Data(name)
    if o.recording
        error("the rows of Data(\"$name\") are gathered once the samples are merged: read them " *
              "in the Simulation runTMS returns")
    end
    return haskey(o.opened, Averaged(d)) ? averaged_data(o.opened[Averaged(d)]) :
           DataTable(o.opened[haskey(o.opened, d) ? d : EachSample(d)].rows)
end

# the name of the `Data` a destination opened writes to, `nothing` for another
data_name(d::Data) = d.name
data_name(d::Union{Averaged, EachSample}) = data_name(d.inner)
data_name(_) = nothing

"""
    end_sample!(o)
    rewrite_averaged!(o, marks)

the rows of the sample ended added to the averages of the outputs `o`, see `fold!`; and their
averaged files written whole as the `marks` of a commit have them, as a checkpoint is written
"""
function end_sample!(o::Outputs)
    lock(o.lock) do
        foreach(fold!, (h for h in values(o.opened) if h isa AveragedOpened))
    end
    return nothing
end

function rewrite_averaged!(o::Outputs, marks::AbstractDict)
    for (d, m) in marks
        if m isa AveragedMark
            write_averaged(d, m.sums, m.fixed, m.samples)
        end
    end
    return nothing
end

"""
    take_records!(o)
    clear_records!(o)

the rows the recorders of the outputs `o` of a part hold, `destination => rows` for each that
has some, see `Recorder`, taken away from them; and those rows dropped
"""
function take_records!(o::Outputs)
    lock(o.lock) do
        records = [ d => h.rows for (d, h) in o.opened if h isa Recorder && !isempty(h.rows) ]
        clear_records!(o)
        return records
    end
end

function clear_records!(o::Outputs)
    lock(o.lock) do
        # anew rather than emptied: the marks of the commits hold them
        for h in values(o.opened)
            if h isa Recorder
                h.rows = Tuple{NTuple{3, Int}, Vector{String}, MeasuredRow}[]
            end
        end
    end
    return nothing
end

"""
    replay!(o, records)

the rows of a sample, as `take_records!` took them in a part, written to the outputs `o` of the
directory of the simulation, then added to its averages, see `end_sample!`, as the sample does
when it is run in the process
"""
function replay!(o::Outputs, records::AbstractVector)
    lock(o.lock) do
        for (d, rows) in records
            opened = get!(() -> open(d), o.opened, d)
            for (at, names, row) in rows
                prepare_list!(opened, names)
                put_keyed_row!(opened, at, row)
            end
        end
    end
    end_sample!(o)
    return o
end

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

"""
    data_outputs(dir, format, marks)

outputs in `dir` holding the `Data` of `marks` alone, averaged or for each sample, reopened at
their marks, so that `sim.data` reads them: the files of `marks` are left as they are
"""
function data_outputs(dir::AbstractString, format::AbstractString, marks::AbstractDict)
    o = Outputs(; dir, format, sampled = true)
    for (d, m) in marks
        if !isnothing(data_name(d))
            o.opened[d] = open(d, m)
        end
    end
    return o
end

"""
    within(sums, samples, (atol, rtol), min_samples)

whether the averages `sums` of `samples` samples, at least `min_samples`, have every error within
`atol + rtol * |mean|`, cell by cell, every row holding at least `min_samples` samples: a row of
one sample, whose error is not known, is not
"""
function within(sums::AbstractDict{RowKey, RowSum}, samples::Int, (atol, rtol), min_samples::Int)
    if samples < max(min_samples, 2)
        return false
    end
    for s in values(sums)
        if s.n < min_samples
            return false
        end
        for i in eachindex(s.columns)
            if !s.missing[i]
                _, err = spread(s, i)
                # NaN fails it
                if !(err ≤ atol + rtol * abs(s.shift[i] + s.mean[i]))
                    return false
                end
            end
        end
    end
    return true
end

"""
    within(opened, dir, format, tolerance, min_samples)

whether the averaged destinations the `tolerance` names, `destination => (atol, rtol)`, have
their errors within it, see `within`, `opened` holding them opened or their marks, their files
named relative to `dir`; one not written yet is not
"""
function within(opened::AbstractDict, dir::AbstractString, format::AbstractString,
                tolerance::AbstractVector, min_samples::Int)
    o = Outputs(; dir, format, sampled = true)
    for (d, tol) in tolerance
        h = get(opened, destination(o, d), nothing)
        if !(h isa Union{AveragedOpened, AveragedMark} &&
             within(h.sums, h.samples, (tol.atol, tol.rtol), min_samples))
            return false
        end
    end
    return true
end

