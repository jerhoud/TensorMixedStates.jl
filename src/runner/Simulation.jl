# Simulation, a state, its time and its values as runTMS runs it, output, which measures it into
# the destinations of the run, log_message, which writes to its log, and the logger sending the
# warnings of the package there.

export Simulation, output, log_message, stopped

"""
    simulation_files

the files `runTMS` writes in the directory of a simulation, which no destination may be named
after
"""
const simulation_files = Set(["log", "stop", "restart", "error", "running", "stamp", "description",
                              "prog.jl", "prog_args.json", checkpoint_files...])

"""
    reserved_values

the names no value of a simulation may take: those the measurements read otherwise, `c.state`
and `c.time`, those the run gives, `c.step` and `c.rng`, and `sim.data`
"""
const reserved_values = (:state, :time, :step, :rng, :data)

"""
    RunData(dir, outputs, log, redirect, checkpointer, provided)

what the simulations of a run share: its directory, empty for a run whose output is redirected
to the stream `redirect`, its outputs, its log, its checkpointing, the values its phases
provide, the phase being run and the time it started from, and whether it stopped
"""
mutable struct RunData
    dir::String
    outputs::Outputs
    log::IO
    redirect::Union{Nothing, IO}
    checkpointer::Checkpointer
    provided::Vector{Symbol}
    phase::Int
    phase_time::Number
    stopped::Bool
end

RunData(dir, outputs, log, redirect, checkpointer, provided) =
    RunData(dir, outputs, log, redirect, checkpointer, provided, 0, 0., false)

"""
    Simulation(sim, state[, time = sim.time]; values...)

the simulation `runTMS` runs and returns: a state, its time and its values, which the
measurements read, `sim.energy` as `c.energy`. The form above gives `sim` another state, and
possibly another time and other values, sharing its run. A value is a number, an array, a
state, anything a phase leaves for the measurements and the phases after it; a checkpoint keeps
it, so that it must not be changed in place.

- `sim.state`, `sim.time`: the state and the time
- `sim.step`: the steps the phase being run has done
- `sim.rng`: the random number generator of the run, for every draw of a phase, see `SimData`
- `sim.name`: the value `name`, a state named by `CopyState` for instance
- `sim.data`: the rows gathered in each `Data(name)`, under `name`

# Examples

    sim = runTMS(sim_data)
    rows = sim.data["magnetization"]
"""
struct Simulation{S <: Union{Nothing, AbstractState}}
    state::S
    time::Number
    values::Dict{Symbol, Any}
    run::RunData
end

function Simulation(sim::Simulation, state::Union{Nothing, AbstractState}, time::Number = sim.time;
                    values...)
    for name in keys(values)
        if name in reserved_values
            error("a simulation cannot be given a value named $name, a name the run reserves")
        end
    end
    return Simulation(state, time, merge(getfield(sim, :values), values), getfield(sim, :run))
end

"""
    with_values(sim; values...)

the simulation with the values given, the reserved ones included, which the run gives
"""
with_values(sim::Simulation; values...) =
    Simulation(sim.state, sim.time, merge(getfield(sim, :values), values), getfield(sim, :run))

"""
    data_rows(sim)

the rows gathered in each `Data` of the run of `sim`, by its name
"""
function data_rows(sim::Simulation)
    o = getfield(sim, :run).outputs
    return lock(o.lock) do
        Dict{String, Vector{Vector{Pair{String, Any}}}}(d.name => h.rows for (d, h) in o.opened
                                                         if d isa Data)
    end
end

function getproperty(sim::Simulation, name::Symbol)
    if name === :state || name === :time
        return getfield(sim, name)
    elseif name === :data
        return data_rows(sim)
    end
    values = getfield(sim, :values)
    if !haskey(values, name)
        error("the simulation has no value $name")
    end
    return values[name]
end

propertynames(sim::Simulation) = (:state, :time, :data, keys(getfield(sim, :values))...)

show(io::IO, sim::Simulation) = print(io, "Simulation($(sim.state), $(sim.time), ...)")

"""
    Editing(sim)
    Simulation(c)

the simulation `sim` open to changes, as `Update` hands it to its function, read as a
`Simulation` and changed by assignments: `c.state = st`, `c.time = t`, `c.x = v` for a value,
`c.x = nothing` removing it; and the simulation it makes once changed
"""
mutable struct Editing
    state::Union{Nothing, AbstractState}
    time::Number
    values::Dict{Symbol, Any}
    run::RunData
end

Editing(sim::Simulation) =
    Editing(sim.state, sim.time, copy(getfield(sim, :values)), getfield(sim, :run))

Simulation(c::Editing) =
    Simulation(getfield(c, :state), getfield(c, :time), getfield(c, :values), getfield(c, :run))

getproperty(c::Editing, name::Symbol) = getproperty(Simulation(c), name)

propertynames(c::Editing) = propertynames(Simulation(c))

function setproperty!(c::Editing, name::Symbol, x)
    if name === :state || name === :time
        setfield!(c, name, x)
    elseif name in reserved_values
        error("$name cannot be changed, a name the run reserves")
    elseif isnothing(x)
        delete!(getfield(c, :values), name)
    else
        getfield(c, :values)[name] = x
    end
    return x
end

"""
    stopped(sim)

whether the run `runTMS` returned stopped before the end of its phases, at `max_time`, on the
file `stop` or on an interrupt: running it again resumes it

# Examples

    sim = runTMS(sim_data)
    if stopped(sim)
        println("stopped at time ", sim.time, ", to be resumed")
    end
"""
stopped(sim::Simulation) = getfield(sim, :run).stopped

"""
    run_destination(run, d)

where the run writes what goes to `d`: the stream of a redirected run for a name, `d` otherwise,
a name of a file of the simulation refused
"""
function run_destination(run::RunData, name::AbstractString)
    if !isnothing(run.redirect)
        return Stream(run.redirect, run.outputs.format)
    elseif !isempty(run.dir) && normpath(name) in simulation_files
        error("cannot write to $name, a file of the simulation directory: choose another name")
    end
    return name
end

run_destination(::RunData, d) = d

"""
    measured_values(sim)

the values the measurements of `sim` read, `missing` for those its phases provide that are not
there yet
"""
function measured_values(sim::Simulation)
    values = getfield(sim, :values)
    return merge(Dict{Symbol, Any}(k => missing for k in getfield(sim, :run).provided), values)
end

"""
    output(sim, destination => measurements)
    output(sim, [destination1 => measurements1, ...])

measure each list on the state of the simulation, at its time and with its values, and write it
to its destination: a file name, `"stdout"`, `"stderr"`, `""` or `Data(name)`, see `Destination`
and `measure`. The files are taken in the directory of the simulation, whose own files are
refused.

# Examples

    output(sim, "magnetization" => [X, Z(1)])
    output(sim, [ "energy" => :energy, Data("purity") => Purity ])
"""
function output(sim::Simulation, list::AbstractVector)
    if isempty(list)
        return nothing
    end
    run = getfield(sim, :run)
    items = [ run_destination(run, first(p)) => last(p) for p in list ]
    # called for the state as it is: inferred for any state, the measurements would take
    # minutes to compile
    invokelatest(output!, run.outputs, sim.state, items; time = sim.time, measured_values(sim)...)
    return nothing
end

output(c::Editing, list::AbstractVector) = output(Simulation(c), list)
output(sim::Union{Simulation, Editing}, p::Pair) = output(sim, [p])

"""
    log_message(sim, text)

write the line `text` to the log of the simulation, or to the stream its output is redirected
to, at once
"""
function log_message(sim::Union{Simulation, Editing}, text)
    log_line(getfield(sim, :run), text)
    return nothing
end

function log_line(run::RunData, text)
    println(run.log, text)
    flush(run.log)
end

"""
    SimLogger(run, parent)

the logger of a run: the warnings of this package go to its log, everything else to `parent`,
the logger in place
"""
struct SimLogger <: Logging.AbstractLogger
    run::RunData
    parent::Logging.AbstractLogger
end

"""
    simulation_warning(level, _module)

whether a log message is a warning of this package, which `SimLogger` sends to the log of the
run
"""
simulation_warning(level, _module) =
    level >= Logging.Warn && !isnothing(_module) && Base.moduleroot(_module) === @__MODULE__

Logging.min_enabled_level(l::SimLogger) = min(Logging.Warn, Logging.min_enabled_level(l.parent))

Logging.shouldlog(l::SimLogger, level, _module, group, id) =
    simulation_warning(level, _module) || Logging.shouldlog(l.parent, level, _module, group, id)

Logging.catch_exceptions(l::SimLogger) = Logging.catch_exceptions(l.parent)

function Logging.handle_message(l::SimLogger, level, message, _module, group, id, file, line;
                                kwargs...)
    if simulation_warning(level, _module)
        log_line(l.run, "WARNING: $message")
    elseif level >= Logging.min_enabled_level(l.parent)
        # the level of this logger is the lower of the two: the parent's is checked again
        Logging.handle_message(l.parent, level, message, _module, group, id, file, line; kwargs...)
    end
    return nothing
end
