export Simulation, output

using Printf

struct Simulation
    state::State
    time::Number
    files::Dict{String, IO}
    formats::Tuple{Printf.Format, Printf.Format}
    Simulation(st::State, t::Number = 0, tfmt::String = "%8.4g", dfmt::String = "%12.6g") =
        new(st, t, Dict(), (Printf.Format(tfmt), Printf.Format(dfmt)))
    Simulation(s::Simulation, st::State, t::Number) =
        new(st, t, s.files, s.formats)
end

get_sim_file(sim::Simulation, filename) =
    get!(sim.files, filename) do
        if filename == "stdout" || filename == "-"
            stdout
        elseif filename == ""
            devnull
        elseif filename == "stderr"
            stderr
        else
            open(filename, "a")
        end
    end

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

function output(sim::Simulation, filename, m::Measure)
    file = get_sim_file(sim, filename)
    v = measure(sim.state, m, sim.time)
    foreach(v) do x
        output(file, first(x), sim.time, last(x), sim.formats)
    end
    flush(file)
end