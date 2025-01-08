export output, log_msg

function output_one(file, x::AbstractFloat, format)
    Printf.format(file, format, x)
end

function output_one(file, x, _)
    print(file, x)
end

output(file, header, t, data, formats) =
    output(file, header, t, [data], formats)

function output(file, header, t, data::Vector, formats)
    print(file, header, "\t")
    Printf.format(file, first(formats), t)
    for x in data
        print(file, "\t")
        output_one(file, x, last(formats))
    end 
    println(file)
end

function output(file, header, t, data::Matrix, formats)
    println(file, header)
    for l in 1:size(data, 1)
        output(file, "$header:$l", t, data[l,:], formats)
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

function output(sim::Simulation, measurements::Vector; debug = false, kwargs...)
    if isempty(measurements)
        return
    end
    vals = measure(sim.state, Measure.(last.(measurements)), sim.time; kwargs...)
    files = [ get_sim_file(sim, filename) for filename in first.(measurements) ]
    for (v, f) in zip(vals, files)
        for x in v
            output(f, first(x), sim.time, last(x), sim.formats)
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