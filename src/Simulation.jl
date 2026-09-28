export Simulation, get_sim_file, DataToFrame, data_to_frame

"""
    data_to_frame(data)

return a `DataFrame` object corresponding to the data, with a row for each set of values
measured together, in the order they were measured. The `DataFrames` package must be imported
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
- `outputs`     : the destinations of the measurements and the formats they are written in,
                  see `Outputs`
- `checkpoint`  : the checkpointing machinery, see `Checkpointer`

`sim.data` is the dictionary of the `Data` destinations, each a `Dict` of the series measured
into it, as `data_to_frame` reads them.

A `Simulation` is immutable, and the state is threaded through a run by building a new one
at each step rather than by assigning to a field. The other two fields are shared rather than
copied, on purpose: the second form above hands the new object the very `outputs` and
`checkpoint` of the old one. They are the parts that must not fork — the destinations and
what they hold, and the bookkeeping that says where the run has got to. A copy made while a
phase is running therefore sees, and can advance, the same checkpoint as the simulation it
was made from.
"""
struct Simulation
    state::Union{Nothing, State}
    time::Number
    outputs::Outputs
    checkpoint::Checkpointer
    Simulation(state::Union{Nothing, State}; time::Number = 0., output = nothing,
               time_format::String = default_time_format, data_format::String = default_data_format,
               checkpoint::Checkpointer = Checkpointer()) =
        new(state, time, Outputs(output, time_format, data_format), checkpoint)
    Simulation(s::Simulation, st::Union{Nothing, State}, t::Number = s.time) =
        new(st, t, s.outputs, s.checkpoint)
end

function Base.getproperty(s::Simulation, f::Symbol)
    if f === :data
        return getfield(s, :outputs).data
    end
    return getfield(s, f)
end

Base.propertynames(::Simulation) = (fieldnames(Simulation)..., :data)

show(io::IO, s::Simulation) = print(io, "Simulation($(s.state), $(s.time), ...)")

"""
    simulation_files

the files `runTMS`, the log and the checkpoints write in the directory of a simulation, which
no destination may be named after: a destination called `stop` stopped the simulation at its
first sweep and was erased by the next run, and one called `checkpoint.json` overwrote the
checkpoint.
"""
const simulation_files = Set(["log", "stop", "error", "running", "stamp", "description",
    "prog.jl", "checkpoint.json", "checkpoint.json.tmp", "checkpoint-1.h5", "checkpoint-2.h5"])

# a simulation with a directory, the one `runTMS` writes in, keeps its files for itself
function check_destination(sim::Simulation, name::AbstractString)
    if !isempty(sim.checkpoint.dir) && normpath(name) in simulation_files
        error("cannot write to $name, a file of the simulation directory: choose another name")
    end
end
check_destination(::Simulation, ::Data) = nothing

"""
    get_sim_file(::Simulation, filename)

return the corresponding file of the given simulation "stdout" (or "-"), "stderr" and "" respectively
redirect to stdout, stderr and devnull, other names are interpreted as file names.

Filename finishing by ".json" will return a Dict
where to store data and this data will be output in JSON format in the file by `runTMS` at the end.

Special filenames of the form `Data(name)` return a Dict where to store Data.
Those Dict are gathered as a Dict in the `data` field of the Simulation

In a simulation run by `runTMS` in its directory, a file of that directory, the log, the
checkpoint and the markers, see `simulation_files`, cannot be asked for.
"""
function get_sim_file(sim::Simulation, name::Union{AbstractString, Data})
    check_destination(sim, name)
    return handle(destination(sim.outputs, name))
end

"""
    close_sim_files(::Simulation)

write the dictionaries collected for the json destinations and close the files opened for
the simulation. The standard streams are destinations like any other but they belong to
the process, so they are left alone.
"""
close_sim_files(sim::Simulation) = close_outputs!(sim.outputs)

length(sim::Simulation) = length(sim.state)

maxlinkdim(sim::Simulation) = maxlinkdim(sim.state)

truncate(sim::Simulation; kwargs...) = Simulation(sim, truncate(sim.state; kwargs...))

mix(sim::Simulation) = Simulation(sim, mix(sim.state))

weaken(sim::Simulation) = Simulation(sim, weaken(sim.state))

weaken(sim::Simulation, spec) = Simulation(sim, weaken(sim.state, spec))

apply(op, sim::Simulation; kwargs...) = Simulation(sim, apply(op, sim.state; kwargs...))

PreMPO(sim::Simulation, args...) = PreMPO(sim.state, args...)

partial_trace(sim::Simulation, pos; kwargs...) =
    Simulation(sim, partial_trace(sim.state, pos; kwargs...))