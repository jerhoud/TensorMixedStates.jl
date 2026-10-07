# Simulation, a state with its simulation time and the destinations of its measurements, as
# runTMS runs it, and data_to_frame, which turns the values gathered in a Data into a table.

export Simulation, get_sim_file, close_sim_files, data_to_frame

"""
    data_to_frame(data)

a `DataFrame` of the values gathered in a `Data` destination, `sim.data[name]`: a `time`
column and a column per measurement, in the order of their names, and a row per call of
`output`, in the order of the calls. The `DataFrames` package must be loaded.

# Examples

    using DataFrames
    sim = runTMS(sim_data)
    df = data_to_frame(sim.data["magnetization"])
"""
function data_to_frame end

# data_to_frame lives in DataFramesExt: without DataFrames, its MethodError says so
function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, e, _, _
        if e.f === data_to_frame && isnothing(Base.get_extension(@__MODULE__, :DataFramesExt))
            print(io, "\ndata_to_frame needs the DataFrames package: run using DataFrames")
        end
    end
end

"""
    default_time_format

the C like format simulation times are written in by default, on which `Simulation` and
`SimData` have to agree, as on `default_data_format`
"""
const default_time_format = "%8.4g"

"""
    default_data_format

the C like format measured values are written in by default, see `default_time_format`
"""
const default_data_format = "%14.8g"

"""
    Simulation(state; time = 0., output = nothing, time_format, data_format)
    Simulation(sim, state[, time = sim.time])

a state with its simulation time and the destinations of its measurements, which `runTMS`
returns. The first form builds one, `output` being a stream every destination is redirected
to, as for `runTMS`. The second gives `sim` another state, and possibly another time, sharing
its destinations and its checkpoint.

A `Simulation` is immutable: the functions acting on it return a new one. `length`,
`maxlinkdim`, `truncate`, `mix`, `weaken`, `apply`, `partial_trace`, `collapse`, `PreMPO`,
`tdvp`, `approx_W`, `dmrg`, `steady_state`, `thermal_state` and `set_threading` take a
`Simulation` as they take a `State`. The measurement functions do not: measure a simulation with `output`,
or its state with `measure(sim.state, measurements, sim.time)`.

# Fields

- `state`: the state of the system
- `time`: the simulation time

The other fields are internal. `sim.data` is the dictionary of the `Data` destinations, see
`Data`.

# Examples

    sim = Simulation(state)
    sim = tdvp(-im * H, 1., sim; nsteps = 10)
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
        new(state, time, Outputs(output, time_format, data_format, checkpoint.dir), checkpoint)
    Simulation(s::Simulation, st::Union{Nothing, AbstractState}, t::Number = s.time) =
        new(st, t, s.outputs, s.checkpoint)
end

function Base.getproperty(s::Simulation, f::Symbol)
    if f === :data
        return data_series(getfield(s, :outputs))
    end
    return getfield(s, f)
end

Base.propertynames(::Simulation) = (fieldnames(Simulation)..., :data)

show(io::IO, s::Simulation) = print(io, "Simulation($(s.state), $(s.time), ...)")

"""
    simulation_files

the files `runTMS`, the log and the checkpoints write in the directory of a simulation, which
no destination may be named after
"""
const simulation_files = Set(["log", "stop", "error", "running", "stamp", "description",
    "prog.jl", "prog_args.json", basename(checkpoint_json("")), basename(checkpoint_json("")) * ".tmp",
    state_file(1), state_file(2)])

"""
    check_destination(::Simulation, name)

refuse a destination named after one of the `simulation_files` in a simulation with a
directory
"""
function check_destination(sim::Simulation, name::AbstractString)
    if !isempty(sim.checkpoint.dir) && normpath(name) in simulation_files
        error("cannot write to $name, a file of the simulation directory: choose another name")
    end
end
check_destination(sim::Simulation, d::Union{TextFile, JsonFile}) = check_destination(sim, d.name)
check_destination(sim::Simulation, ::LogFile) = check_destination(sim, "log")
check_destination(::Simulation, ::Destination) = nothing

"""
    save_state(filename, statename, ::Simulation)

save the state of the simulation as `save_state` saves a state, refusing a `filename` that
names a file of the simulation directory, as its checkpoint.

# Examples

    save_state("ground.h5", "gs", sim)
"""
function save_state(filename::String, statename::String, sim::Simulation)
    check_destination(sim, filename)
    return save_state(in_dir(sim.checkpoint.dir, filename), statename, sim.state)
end

"""
    get_sim_file(::Simulation, name)

the stream `output` writes to under this name, to write to it directly. `"stdout"` (or
`"-"`), `"stderr"` and `""` give `stdout`, `stderr` and `devnull`, any other name the stream
of a file. When the output of the simulation is redirected, every name gives that stream. A
name ending in `.json` and a `Data`, whose values `output` alone writes, are refused, as are
the files `runTMS` writes itself in its directory, the log, the checkpoint and the markers.

It is meant for a phase of your own, see `TensorMixedStates.run_phase`, to write what is not a
measurement: unlike a file opened with `open`, such a file is cut back on a resume as the
others are. Once `runTMS` has returned, the files of the simulation are closed.

# Examples

    function TensorMixedStates.run_phase(sim::Simulation, p::MyPhase)
        println(get_sim_file(sim, "notes.txt"), "starting at time ", sim.time)
        return sim
    end
"""
function get_sim_file(sim::Simulation, name::Union{AbstractString, Destination})
    d = Destination(name)
    check_destination(sim, d)
    return stream(sink(sim.outputs, d))
end

"""
    close_sim_files(::Simulation)

write the json destinations and close the files opened for the simulation, the standard
streams being left open. `runTMS` calls it when it ends; a `Simulation` built by hand calls
it once its measurements are made, its json files being written only then.

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

function collapse(sim::Simulation, args...; kwargs...)
    x, st = collapse(sim.state, args...; kwargs...)
    return (x, Simulation(sim, st))
end