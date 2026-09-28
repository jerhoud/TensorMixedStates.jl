export output, log_msg

# a row of a text file, written the way `output` writes a measurement
output(sim::Simulation, file::IO, header, data) =
    write_row(file, sim.outputs.formats, sim.time, header, data)

"""
    output(::Simulation, [ filename => measure1, ... ])

compute the given measurements on a simulation and output them to the associated file or dict

filenames are interpreted by get\\_sim\\_file (see there for special values)

A complex value takes two columns, its real part then its imaginary part, and a json file
writes it as `{"re": …, "im": …}`: see `RealValue` for which values are complex.

# Examples

    output(sim, "file" => [X, X(1)Y(2), (X, Y)])
    output(sim, [ "file1" => [X, Y(2)], "file2" => Trace])
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
        emit!(destination(sim.outputs, name), sim.outputs.formats, sim.time, v)
    end
end

"""
    log_msg(::Simulation, text)

log the given message on the "log" file of the simulation
"""
# written here rather than through `output`, where `dest => "text"` is a measurement
log_msg(sim::Simulation, text) = emit_line!(destination(sim.outputs, "log"), text)

"""
    SimLogger(sim, parent)

the logger `output` measures under. The warnings of this package, `measure` dropping a part
of a value that is more than rounding, go to the log of the simulation with the rest of what
it reports, and everything else goes on to `parent`, the logger in place, as it would outside
a measurement. `measure` itself only warns, so that a direct call shows its warnings as any
other would.
"""
struct SimLogger <: Logging.AbstractLogger
    sim::Simulation
    parent::Logging.AbstractLogger
end

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