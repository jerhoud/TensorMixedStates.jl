# output, which measures a simulation and writes the values to their destinations, log_message,
# and the logger sending the warnings of the package to the log of the simulation.

export output, log_message

"""
    output(::Simulation, destination => measurements; kwargs...)
    output(::Simulation, [ destination1 => measurements1, ... ]; kwargs...)

compute the given measurements on a simulation at its time, in a single call of `measure`,
and write them to their destinations. The keyword arguments give their values to the `Symbol`
measurements, as for `measure`. A destination is a name, a text file but for `"stdout"` (or
`"-"`), `"stderr"`, `""` for `devnull` and a name ending in `.json` for a json file, or a
`Data(name)`; the measurements are anything `measure` takes.

A text file takes one line per measurement: its name, the time and its values, separated by
tabs. A matrix takes a line with its name alone, then one line per row, named `name:l` for row
`l`. A complex value, see `RealValue`, takes two columns, real then imaginary part, and is
written `{"re": …, "im": …}` in a json file.

A json file is written when the files of the simulation are closed, by `runTMS` when it ends,
or by `close_sim_files` for a `Simulation` built by hand.

# Examples

    output(sim, "file" => [X, X(1)Y(2), (X, Y)])
    output(sim, [ "file1" => [X, Y(2)], "file2" => Trace])
    output(sim, Data("magnetization") => Z)
"""
output(sim::Simulation, m::Pair; kwargs...) =
    output(sim, [m]; kwargs...)

function output(sim::Simulation, measurements::Vector; kwargs...)
    if isempty(measurements)
        return
    end
    vals = Logging.with_logger(SimLogger(sim, Logging.current_logger())) do
        measure(sim.state, Measure.(last.(measurements)), sim.time; kwargs...)
    end
    # one event per destination for the whole call, so that data_to_frame makes one row of
    # the pairs going to one Data
    events = IdDict()
    for (v, name) in zip(vals, first.(measurements))
        d = Destination(name)
        check_destination(sim, d)
        s = sink(sim.outputs, d)
        event = get!(() -> new_event(s), events, s)
        emit!(s, sim.outputs.formats, sim.time, v; event)
    end
end

"""
    measure_sets(measurements)

the measurements given as `output` takes them, each destination with its `Measure`, for an
observer to simplify and compact the operators once for its phase
"""
measure_sets(m::Pair) = measure_sets([m])
measure_sets(ms::Vector) = Pair[ first(m) => Measure(last(m)) for m in ms ]

"""
    log_message(::Simulation, text)

write the line `text` to the `log` file of the simulation, or to the stream its output is
redirected to, flushed at once. A `Simulation` built by hand without `output` has its log in a
file `log` of the current directory, which its first line empties.
"""
function log_message(sim::Simulation, text)
    # not through `output`, where `dest => "text"` is a measurement
    emit_line!(sink(sim.outputs, LogFile()), text)
end

"""
    SimLogger(sim, parent)

the logger `output` measures under: the warnings of this package go to the log of the
simulation, everything else to `parent`, the logger in place. `measure` itself only warns, so
that a direct call shows its warnings as any other.
"""
struct SimLogger <: Logging.AbstractLogger
    sim::Simulation
    parent::Logging.AbstractLogger
end

"""
    simulation_warning(level, _module)

whether a log message is a warning of this package, which `SimLogger` sends to the log of the
simulation
"""
simulation_warning(level, _module) = level >= Logging.Warn && _module === @__MODULE__

Logging.min_enabled_level(l::SimLogger) = min(Logging.Warn, Logging.min_enabled_level(l.parent))

Logging.shouldlog(l::SimLogger, level, _module, group, id) =
    simulation_warning(level, _module) || Logging.shouldlog(l.parent, level, _module, group, id)

Logging.catch_exceptions(l::SimLogger) = Logging.catch_exceptions(l.parent)

function Logging.handle_message(l::SimLogger, level, message, _module, group, id, file, line; kwargs...)
    if simulation_warning(level, _module)
        log_message(l.sim, "WARNING: $message")
    elseif level >= Logging.min_enabled_level(l.parent)
        # the level of this logger is the lower of the two: the parent's is checked again
        Logging.handle_message(l.parent, level, message, _module, group, id, file, line; kwargs...)
    end
    return nothing
end