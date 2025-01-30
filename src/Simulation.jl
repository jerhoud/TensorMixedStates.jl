export Simulation, get_sim_file

"""
    Simulation(state[, time = 0])
    Simulation(sim, state[, time = sim.time])

A type to represent simulation data ans store time and file data. It is used and returned by runTMS.
The first form creates a simulation object. The second updates the state in the simulation object. (see also `get_sim_file`)

# Fields
- `state`
- `time`
"""
struct Simulation
    state
    time::Number
    output::Union{Nothing, IO}
    files::Dict{String, IO}
    formats::Tuple{Printf.Format, Printf.Format}
    Simulation(state; t::Number = 0, output = nothing, time_format::String = "%8.4g", data_format::String = "%12.6g") =
        new(state, t, output, Dict(), (Printf.Format(time_format), Printf.Format(data_format)))
    Simulation(s::Simulation, st, t::Number = s.time) =
        new(st, t, s.output, s.files, s.formats)
end

"""
    get_sim_file(::Simulation, filename)

return the corresponding file of the given simulation
"""
get_sim_file(sim::Simulation, filename) =
    get!(sim.files, filename) do
        if !isnothing(sim.output)
            sim.output
        elseif filename == "stdout" || filename == "-"
            stdout
        elseif filename == ""
            devnull
        elseif filename == "stderr"
            stderr
        else
            open(filename, "w")
        end
    end

length(sim::Simulation) = length(sim.state)

maxlinkdim(sim::Simulation) = maxlinkdim(sim.state)

truncate(sim::Simulation; kwargs...) = Simulation(sim, truncate(sim.state; kwargs...))

mix(sim::Simulation) = Simulation(sim, mix(sim.state))

apply(op, sim::Simulation; kwargs...) = Simulation(sim, apply(op, sim.state; kwargs...))

PreMPO(sim::Simulation, tp) = PreMPO(sim.state, tp)

tdvp(op, t::Number, sim::Simulation; kwargs...) =
    Simulation(sim, tdvp(op, t, sim.state; time_start = sim.time, kwargs...), sim.time + t)


function dmrg(op, sim::Simulation; kwargs...)
    e, st = dmrg(op, sim.state; kwargs...)
    return (e, Simulation(sim, st))
end

approx_W(op, t::Number, sim::Simulation; kwargs...) =
    Simulation(sim, approx_W(op, t, sim.state; time_start = sim.time, kwargs...), sim.time + t)
