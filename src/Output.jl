export output, log_msg

const thresh_warn_imag = 1e-6

make_real(::Simulation, _, data) = data

warn_imag(sim::Simulation, header, x) =
    log_msg(sim, "WARNING: large imaginary part: time $(sim.time), $header " * @sprintf("%8.1e (rel %8.1e)", imag(x), imag(x) / real(x)))

function make_real(sim::Simulation, header, data::Number)
    if abs(imag(data)) > thresh_warn_imag
        warn_imag(sim, header, data)
    end
    return real(data)
end

function make_real(sim::Simulation, header, data::Union{Vector, Matrix})
    map(keys(data)) do ij
        x = data[ij]
        if !(x isa Number)
            return X
        end
        if abs(imag(x)) > thresh_warn_imag
            warn_imag(sim, "$header $(Tuple(ij))", x)
        end
        return real(x)
    end
end

function output_one(file, x::AbstractFloat, format)
    Printf.format(file, format, x)
end

function output_one(file, x::Complex, format)
    output_one(file, real(x), format)
    print(file, "\t")
    output_one(file, imag(x), format)
end

function output_one(file, x, _)
    print(file, x)
end

output(sim::Simulation, file::IO, header, data) =
    output(sim, file, header, [data])

function output(sim::Simulation, file::IO, header, data::Vector)
    print(file, header, "\t")
    Printf.format(file, first(sim.formats), sim.time)
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

compute the given measurements on a simultation and output them to the associated file or dict

filenames are interpreted by get\\_sim\\_file (see there for special values)

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
    vals = measure(sim.state, Measure.(last.(measurements)), sim.time; kwargs...)
    files = [ get_sim_file(sim, filename) for filename in first.(measurements) ]
    for (v, f) in zip(vals, files)
        for x in v
            header = first(x)
            data = make_real(sim, header, last(x))
            output(sim, f, header, data)
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