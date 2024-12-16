export runTMS, SimData

@kwdef struct SimData
    description::String=""
    name::String=""
    end_measures = []
    tfmt::String = "%8.4g"
    dfmt::String = "%12.6g"
    phases
end

function runTMS(sim_data::SimData; restart::Bool=false, clean::Bool=false, debug::Bool=false)
    live = !debug
    if !debug && (restart || clean)
        rm(sim_data.name; recursive = true)
    end
    if clean 
        return
    end
    try
        start_dir = pwd()
        if !debug
            mkpath(sim_data.name);
            cd(sim_data.name);
            touch("running")
            if sim_data.description ≠ ""
                write("description", sim_data.description)
            end
        end
        sim = Simulation(nothing; debug, sim_data.tfmt, sim_data.dfmt)
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
    _, elapsed, bytes = @timed begin
        sim = run_phase(sim, phase)
        output(sim, phase.end_measures)
    end
    elapsed = round(elapsed; digits=3)
    log_msg(sim, "***** Ending phase \"$(phase.name)\" after $elapsed seconds, $(Base.format_bytes(bytes)) allocated *****")
    return sim
end
