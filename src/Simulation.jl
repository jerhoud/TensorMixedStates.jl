export Simulation, get_sim_file, Data, DataToFrame, data_to_frame

"""
    Data(name)

represent a storage with the given name where to put measurement data
"""
struct Data
    name::String
end

"""
    data_to_frame(data)

return a `DataFrame` object corresponding to the data. The `DataFrames` package must be imported
before using this function.
"""
function data_to_frame end

# the docstring goes through `@doc` rather than sitting above the call, because the macro
# expands to a toplevel block and a docstring cannot be attached to one
Base.@deprecate DataToFrame(data) data_to_frame(data) false

@doc """
    DataToFrame(data)

deprecated, use [`data_to_frame`](@ref) instead, which spells it the way the other
functions of the package are spelled.
""" DataToFrame

"""
    default_time_format
    default_data_format

the C like formats used to write simulation times and measured values. `Simulation` and
`SimData` both default to them, and have to agree: a `Simulation` built by `runTMS` is
given the formats of the `SimData`, one built directly falls back to these.
"""
const default_time_format = "%8.4g"
const default_data_format = "%14.8g"

"""
    Simulation(state[; time = 0.])
    Simulation(sim, state[, time = sim.time])

A type to represent simulation data and store time and file data. It is used and returned by runTMS.
The first form creates a simulation object. The second updates the state in the simulation object. (see also `get_sim_file`)

Most functions applicable to States can be applied to Simulations

# Fields
- `state`       : the state of the system
- `time`        : the simulation time
- `output`      : if not nothing an io where to redirect output
- `files`       : a dictionary holding io or dict where to write data
- `data`        : a dictionary holding data collected for the `Data` objects
- `formats`     : format info for the output
- `checkpoint`  : the checkpointing machinery, see `Checkpointer`

A `Simulation` is immutable, and the state is threaded through a run by building a new one
at each step rather than by assigning to a field. Three of the fields are shared rather than
copied, on purpose: the second form above hands the new object the very `files`, `data` and
`checkpoint` of the old one. They are the parts that must not fork — the open output files,
the accumulated data, and the bookkeeping that says where the run has got to. A copy made
while a phase is running therefore sees, and can advance, the same checkpoint as the
simulation it was made from.
"""
struct Simulation
    state::Union{Nothing, State}
    time::Number
    output::Union{Nothing, IO}
    files::Dict{String, Union{IO, Dict}}
    data::Dict{String, Dict}
    formats::Tuple{Printf.Format, Printf.Format}
    checkpoint::Checkpointer
    Simulation(state::Union{Nothing, State}; time::Number = 0., output = nothing,
               time_format::String = default_time_format, data_format::String = default_data_format,
               checkpoint::Checkpointer = Checkpointer()) =
        new(state, time, output, Dict(), Dict(), (Printf.Format(time_format), Printf.Format(data_format)), checkpoint)
    Simulation(s::Simulation, st::Union{Nothing, State}, t::Number = s.time) =
        new(st, t, s.output, s.files, s.data, s.formats, s.checkpoint)
end

show(io::IO, s::Simulation) = print(io, "Simulation($(s.state), $(s.time), ...)")

"""
    get_sim_file(::Simulation, filename)

return the corresponding file of the given simulation "stdout" (or "-"), "stderr" and "" respectively
redirect to stdout, stderr and devnull, other names are interpreted as file names.

Filename finishing by ".json" will return a Dict
where to store data and this data will be output in JSON format in the file by `runTMS` at the end.

Special filenames of the form `Data(name)` return a Dict where to store Data.
Those Dict are gathered as a Dict in the `data` field of the Simulation
"""
get_sim_file(sim::Simulation, filename::AbstractString) =
    if !isnothing(sim.output)
        sim.output
    else
        get!(sim.files, filename) do
            if filename == "stdout" || filename == "-"
                stdout
            elseif filename == ""
                devnull
            elseif filename == "stderr"
                stderr
            elseif last(splitext(filename)) == ".json"
                Dict()
            else
                open(filename, sim.checkpoint.appending ? "a" : "w")
            end
        end
    end

get_sim_file(sim::Simulation, data::Data) =
    get!(sim.data, data.name) do
        Dict()
    end

"""
    close_sim_files(::Simulation)

write the dictionaries collected for the json destinations and close the files opened for
the simulation. The standard streams are destinations like any other but they belong to
the process, so they are left alone, the same way `save_checkpoint` leaves them alone.
"""
function close_sim_files(sim::Simulation)
    for (filename, data) in sim.files
        if data isa Dict
            open(filename, "w") do io
                JSON.print(io, data)
            end
        elseif data ∉ (stdout, stderr, devnull)
            close(data)
        end
    end
    return nothing
end

length(sim::Simulation) = length(sim.state)

maxlinkdim(sim::Simulation) = maxlinkdim(sim.state)

truncate(sim::Simulation; kwargs...) = Simulation(sim, truncate(sim.state; kwargs...))

mix(sim::Simulation) = Simulation(sim, mix(sim.state))

weaken(sim::Simulation) = Simulation(sim, weaken(sim.state))

weaken(sim::Simulation, spec) = Simulation(sim, weaken(sim.state, spec))

apply(op, sim::Simulation; kwargs...) = Simulation(sim, apply(op, sim.state; kwargs...))

PreMPO(sim::Simulation, args...) = PreMPO(sim.state, args...)

tdvp(op, t::Number, sim::Simulation; kwargs...) =
    Simulation(sim, tdvp(op, t, sim.state; time_start = sim.time, kwargs...), sim.time + t)


function dmrg(op, sim::Simulation; kwargs...)
    e, st = dmrg(op, sim.state; kwargs...)
    return (e, Simulation(sim, st))
end

approx_W(op, t::Number, sim::Simulation; kwargs...) =
    Simulation(sim, approx_W(op, t, sim.state; time_start = sim.time, kwargs...), sim.time + t)

function steady_state(op, sim::Simulation; kwargs...)
    e, st = steady_state(op, sim.state; kwargs...)
    return (e, Simulation(sim, st))
end

partial_trace(sim::Simulation, pos; kwargs...) =
    Simulation(sim, partial_trace(sim.state, pos; kwargs...))