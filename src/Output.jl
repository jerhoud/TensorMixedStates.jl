export output, log_msg

function output_one(file, x::AbstractFloat, format)
    Printf.format(file, format, x)
end

function output_one(file, x::Complex, format)
    output_one(file, real(x), format)
    print(file, "\t")
    output_one(file, imag(x), format)
end

# a value holding several, the part of a `Check` made on a vector observable for instance,
# is written number by number, a matrix row by row as the lines of a matrix are, rather than
# as the literal Julia prints
function output_one(file, x::AbstractArray, format)
    for (k, y) in enumerate(row_major(x))
        if k > 1
            print(file, "\t")
        end
        output_one(file, y, format)
    end
end

function output_one(file, x, _)
    print(file, x)
end

# the time column always takes the time format, unlike a measured value, which keeps its
# own printed form when it is not a float: `Linkdim` is meant to read as 8, not as 8.000
function output_time(file, t::Number, format)
    Printf.format(file, format, t)
end

# `Printf` refuses a complex number outright, so a complex simulation time is written as
# the two columns a complex measurement takes, real then imaginary
function output_time(file, t::Complex, format)
    output_time(file, real(t), format)
    print(file, "\t")
    output_time(file, imag(t), format)
end

output(sim::Simulation, file::IO, header, data) =
    output(sim, file, header, [data])

function output(sim::Simulation, file::IO, header, data::Vector)
    print(file, header, "\t")
    output_time(file, sim.time, first(sim.formats))
    for x in data
        print(file, "\t")
        output_one(file, x, last(sim.formats))
    end
    println(file)
end

function output(sim::Simulation, file::IO, header, data::Matrix)
    println(file, header)
    for l in 1:size(data, 1)
        output(sim, file, "$header:$l", data[l,:])
    end
end

function output(sim::Simulation, dict::Dict, header, data)
    t = sim.time
    d = get!(dict, header, Dict("times"=>[], "data"=>[]))
    push!(d["times"], t)
    push!(d["data"], data)
end

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
    files = [ get_sim_file(sim, filename) for filename in first.(measurements) ]
    for (v, f) in zip(vals, files)
        for x in v
            output(sim, f, first(x), last(x))
        end
        if f isa IO
            flush(f)
        end
    end
end

function output(sim::Simulation, text::Pair{<:Any, <:AbstractString})
    file = get_sim_file(sim, first(text))
    println(file, last(text))
    flush(file)
end

"""
    log_msg(::Simulation, text)

log the given message on the "log" file of the simulation
"""
log_msg(sim::Simulation, text) = output(sim, "log" => text)

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