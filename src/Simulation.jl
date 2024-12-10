export Simulation

struct Simulation
    state::State
    time::Number
    files::Dict{String, IO}
    formats::Tuple{Printf.Format, Printf.Format}
    Simulation(st::State, t::Number = 0, tfmt::String = "%8.4g", dfmt::String = "%12.6g") =
        new(st, t, Dict(), (Printf.Format(tfmt), Printf.Format(dfmt)))
    Simulation(s::Simulation, st::State, t::Number = s.time) =
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
