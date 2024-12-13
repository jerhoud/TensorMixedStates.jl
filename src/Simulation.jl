export Simulation, get_sim_file

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

length(sim::Simulation) = length(sim.state)

maxlinkdim(sim::Simulation) = maxlinkdim(sim.state)

truncate(sim::Simulation; kwargs...) = Simulation(sim, truncate(sim.state; kwargs...))

mix(sim::Simulation) = Simulation(sim, mix(sim.state))

apply(op, sim::Simulation; kwargs...) = Simulation(sim, apply(op, sim.state; kwargs...))

tdvp(s::Simulation, mpo::MPO, t::Number, n::Int; kwargs...) =
    Simulation(s, tdvp(s.state, mpo, t, n; kwargs...), s.time + t)

function dmrg(s::Simulation, mpo::MPO, n::Int; kwargs...)
    e, st = dmrg(s.state, mpo, n; kwargs...)
    return (e, Simulation(s, st))
end
        
PreMPO(sim::Simulation, tp) = PreMPO(sim.state, tp)