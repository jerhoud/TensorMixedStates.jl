export runTMS, SimData

@kwdef struct SimData
    description::String=""
    name::String=""
    time_start = 0.
    final_measures = []
    time_format::String = "%8.4g"
    data_format::String = "%12.6g"
    phases
end

function runTMS(sim_data::SimData; restart::Bool=false, clean::Bool=false, debug::Bool=false)
    live = !debug
    if !debug && (restart || clean)
        rm(sim_data.name; recursive = true, force = true)
    end
    if clean 
        return
    end
    start_dir = pwd()
    try
        if !debug
            mkpath(sim_data.name);
            cd(sim_data.name);
            touch("running")
            if sim_data.description ≠ ""
                write("description", sim_data.description)
            end
        end
        sim = Simulation(nothing; debug, sim_data.time_format, sim_data.data_format)
        log_phase(sim, sim_data)
        if !debug
            rm("running")
            cd(start_dir)
        end
    catch
        if !debug
            touch("error")
            rm("running"; force = true)
            cd(start_dir)
        end
        rethrow()
    end
end

function log_phase(sim::Simulation, phase)
    log_msg(sim, "\n***** Starting phase \"$(phase.name)\" *****")
    if !isnothing(phase.time_start)
        sim = Simulation(sim, sim.state, phase.time_start)
    end
    _, elapsed, bytes = @timed begin
        sim = run_phase(sim, phase)
        output(sim, phase.final_measures)
    end
    elapsed = round(elapsed; digits=3)
    log_msg(sim, "***** Ending phase \"$(phase.name)\" after $elapsed seconds, $(Base.format_bytes(bytes)) allocated *****")
    return sim
end
