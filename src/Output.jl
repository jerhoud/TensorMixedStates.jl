# output, which measures a simulation and writes the values to their destinations, log_msg,
# which writes to its log, and the logger sending the warnings of the package to that log.

export output, log_msg

"""
    output(::Simulation, destination => measurements)
    output(::Simulation, [ destination1 => measurements1, ... ])

compute the given measurements on a simulation, at its time and in a single call of
`measure`, and write them to their destinations. A destination is a name, read as
`get_sim_file` reads it, or a `Data(name)`; the measurements are anything `measure` takes.

A complex value takes two columns, its real part then its imaginary part, and a json file
writes it as `{"re": …, "im": …}`: see `RealValue` for which values are complex.

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
    for (v, name) in zip(vals, first.(measurements))
        check_destination(sim, name)
        emit!(destination(sim.outputs, name), sim.outputs.formats, sim.time, v)
    end
end

"""
    log_msg(::Simulation, text)

write the given line to the `log` file of the simulation, or to the stream its output is
redirected to, flushed at once.
"""
function log_msg(sim::Simulation, text)
    # written here rather than through `output`, where `dest => "text"` is a measurement. A
    # comment between a docstring and what it documents detaches it, so this one is inside
    emit_line!(destination(sim.outputs, "log"), text)
end

"""
    SimLogger(sim, parent)

the logger `output` measures under. The warnings of this package, such as `measure`
dropping a part of a value that is more than rounding, go to the log of the simulation, and
everything else goes on to `parent`, the logger in place. `measure` itself only warns, so that
a direct call shows its warnings as any other.
"""
struct SimLogger <: Logging.AbstractLogger
    sim::Simulation
    parent::Logging.AbstractLogger
end

"""
    simulation_warning(level, _module)

whether a log message is a warning of this package, which `SimLogger` sends to the log of the
simulation.
"""
simulation_warning(level, _module) = level >= Logging.Warn && _module === @__MODULE__

Logging.min_enabled_level(l::SimLogger) = min(Logging.Warn, Logging.min_enabled_level(l.parent))

Logging.shouldlog(l::SimLogger, level, _module, group, id) =
    simulation_warning(level, _module) || Logging.shouldlog(l.parent, level, _module, group, id)

Logging.catch_exceptions(l::SimLogger) = Logging.catch_exceptions(l.parent)

function Logging.handle_message(l::SimLogger, level, message, _module, group, id, file, line; kwargs...)
    if simulation_warning(level, _module)
        log_msg(l.sim, "WARNING: $message")
    elseif level >= Logging.min_enabled_level(l.parent)
        # the level of this logger is the lower of the two, so the parent's has to be
        # checked again here
        Logging.handle_message(l.parent, level, message, _module, group, id, file, line; kwargs...)
    end
    return nothing
end