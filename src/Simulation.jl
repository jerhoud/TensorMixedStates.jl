# Simulation, a state with its simulation time and the destinations of its measurements, as
# runTMS runs it, and data_to_frame, which turns the values gathered in a Data into a table.

export Simulation, get_sim_file, close_sim_files, DataToFrame, data_to_frame

"""
    data_to_frame(data)

a `DataFrame` of the values gathered in a `Data` destination, `sim.data[name]`: a `time`
column and a column for each measurement, in the order of their names, with a row for each
call of `output`, in the order of the calls. The `DataFrames` package must be loaded.

# Examples

    using DataFrames
    sim = runTMS(sim_data)
    df = data_to_frame(sim.data["magnetization"])
"""
function data_to_frame end

# data_to_frame lives in the extension DataFramesExt, loaded with DataFrames: without it, a call
# is a MethodError, which says so
function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, e, _, _
        if e.f === data_to_frame && isnothing(Base.get_extension(@__MODULE__, :DataFramesExt))
            print(io, "\ndata_to_frame needs the DataFrames package: run using DataFrames")
        end
    end
end

# the docstring goes through `@doc` rather than sitting above the call, because the macro
# expands to a toplevel block and a docstring cannot be attached to one
Base.@deprecate DataToFrame(data) data_to_frame(data) false

@doc """
    DataToFrame(data)

deprecated, use [`data_to_frame`](@ref) instead.
""" DataToFrame

"""
    default_time_format

the C like format simulation times are written in by default. `Simulation` and `SimData` both
default to it and to `default_data_format`, and have to agree: a `Simulation` built by
`runTMS` is given the formats of the `SimData`, one built directly falls back to these.
"""
const default_time_format = "%8.4g"

"""
    default_data_format

the C like format measured values are written in by default, see `default_time_format`.
"""
const default_data_format = "%14.8g"

"""
    Simulation(state; time = 0., output = nothing, time_format, data_format)
    Simulation(sim, state[, time = sim.time])

a state with its simulation time and the destinations of its measurements, which `runTMS`
returns. The first form builds one, `output` being a stream every destination is redirected
to, as for `runTMS`. The second gives `sim` another state, and possibly another time.

`length`, `maxlinkdim`, `truncate`, `mix`, `weaken`, `apply`, `partial_trace`, `PreMPO`,
`tdvp`, `approx_W`, `dmrg`, `steady_state` and `set_threading` take a `Simulation` as they take
a `State`. The measurement functions do not: measure a simulation with `output`, or its state
with `measure(sim.state, measurements, sim.time)`.

# Fields

- `state`: the state of the system
- `time`: the simulation time
- `outputs`: the destinations of the measurements and the formats they are written in
- `checkpoint`: the checkpointing machinery

`sim.data` is the dictionary of the `Data` destinations, see `Data`.

A `Simulation` is immutable: the functions acting on it return a new one. The second form
shares the destinations and the checkpoint of `sim` rather than copying them, so that the new
simulation writes to the same destinations and advances the same checkpoint.

# Examples

    sim = Simulation(state)
    sim = tdvp(-im * H, 1., sim; nsweeps = 10)
    output(sim, "data.dat" => [X, Z(1)])
"""
struct Simulation
    state::Union{Nothing, AbstractState}
    time::Number
    outputs::Outputs
    checkpoint::Checkpointer
    Simulation(state::Union{Nothing, AbstractState}; time::Number = 0., output = nothing,
               time_format::String = default_time_format,
               data_format::String = default_data_format,
               checkpoint::Checkpointer = Checkpointer()) =
        new(state, time, Outputs(output, time_format, data_format), checkpoint)
    Simulation(s::Simulation, st::Union{Nothing, AbstractState}, t::Number = s.time) =
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
no destination may be named after: a destination called `stop` would stop the simulation, and
one called `checkpoint.json` would overwrite the checkpoint.
"""
const simulation_files = Set(["log", "stop", "error", "running", "stamp", "description",
    "prog.jl", basename(checkpoint_json("")), basename(checkpoint_json("")) * ".tmp",
    state_file(1), state_file(2)])

"""
    check_destination(::Simulation, name)

refuse a destination named after one of the `simulation_files`, in a simulation with a
directory, the one `runTMS` writes in, which keeps its files for itself.
"""
function check_destination(sim::Simulation, name::AbstractString)
    if !isempty(sim.checkpoint.dir) && normpath(name) in simulation_files
        error("cannot write to $name, a file of the simulation directory: choose another name")
    end
end
check_destination(::Simulation, ::Data) = nothing

"""
    get_sim_file(::Simulation, name)

the destination `output` writes to under this name, to write to it directly. `"stdout"` (or
`"-"`), `"stderr"` and `""` give `stdout`, `stderr` and `devnull`, any other name the stream
of a file of that name. A name ending in `.json` gives instead a `Dict` gathering the data,
written to the file as json when the files of the simulation are closed, see
`close_sim_files`, and `Data(name)` the `Dict` of
`sim.data[name]`. When the output of the simulation is redirected, every name but a `Data`
one gives that stream.

In a simulation run by `runTMS` in its directory, the files `runTMS` writes there itself, the
log, the checkpoint and the markers, cannot be asked for.

It is meant for a phase of your own, see `TensorMixedStates.run_phase`, to write what is not a
measurement in a file of the simulation: that file is cut back on a resume as the others are,
where one opened with `open` would get the lines written since the last checkpoint twice.
Once `runTMS` has returned, the files of the simulation are closed.

# Examples

    function TensorMixedStates.run_phase(sim::Simulation, p::MyPhase)
        println(get_sim_file(sim, "notes.txt"), "starting at time ", sim.time)
        return sim
    end
"""
function get_sim_file(sim::Simulation, name::Union{AbstractString, Data})
    check_destination(sim, name)
    return handle(destination(sim.outputs, name))
end

"""
    close_sim_files(::Simulation)

write the json destinations and close the files opened for the simulation. The standard
streams belong to the process and are left open. `runTMS` calls it when it ends; a
`Simulation` built by hand calls it once its measurements are made, its json files being
written only then.

# Examples

    sim = Simulation(state)
    output(sim, "magnetization.json" => Z)
    close_sim_files(sim)
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