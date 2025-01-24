export output, log_msg

struct LogWarn <: Exception
    msg
end

function output_one(file, x::AbstractFloat, format)
    Printf.format(file, format, x)
end

function output_one(file, x::Complex, format)
    output_one(file, real(x), format)
    if abs(imag(x)) > 1e-12
        throw(LogWarn(@sprintf("WARNING: large imaginary part %8.1e (rel %8.1e)", imag(x), imag(x) / real(x))))
    end
end

function output_one(file, x, _)
    print(file, x)
end

output(sim, file, header, data) =
    output(sim, file, header, [data])

function output(sim::Simulation, file, header, data::Vector)
    print(file, header, "\t")
    Printf.format(file, first(sim.formats), sim.time)
    for (i, x) in enumerate(data)
        try
            print(file, "\t")
            output_one(file, x, last(sim.formats))
        catch e
            if e isa LogWarn
                log_msg(sim, "$(e.msg): $header $(sim.time) ($i)")
            else
                rethrow(e)
            end
        end
    end 
    println(file)
end

function output(sim::Simulation, file, header, data::Matrix)
    println(file, header)
    for l in 1:size(data, 1)
        output(sim, file, "$header:$l", data[l,:])
    end
end

"""
    output(sim, measurments)

compute and output the given measurments on simultaion data sim

# Examples
    output(sim, [X, X(1)Y(2), (X, Y)])
"""
output(sim::Simulation, m::Pair; kwargs...) =
    output(sim, [m]; kwargs...)

function output(sim::Simulation, measurements::Vector; kwargs...)
    if isempty(measurements)
        return
    end
    vals = measure(sim.state, Measure.(last.(measurements)), sim.time; kwargs...)
    files = [ get_sim_file(sim, filename) for filename in first.(measurements) ]
    for (v, f) in zip(vals, files)
        for x in v
            output(sim, f, first(x), last(x))
        end
        flush(f)
    end
end

function output(sim::Simulation, text::Pair{<:Any, <:AbstractString})
    file = get_sim_file(sim, first(text))
    println(file, last(text))
    flush(file)
end

"""
    log_msg(sim, text)

log the given message on the "log" file of the simulation
"""
log_msg(sim::Simulation, text) = output(sim, "log" => text)