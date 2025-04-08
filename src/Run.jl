export runTMS, SimData

"""
    SimData(name = "my_simulation", phases::Vector{Phases} = [phase1, phase2...])

A type for describing a simulation to use with `runTMS`

# Fields

- `name`:            the name of the simulation used as the name of the directory to store the results
- `phases`:          the list of phases of the simulation (see Phases for a list of possible values)
- `descritpion`:     text put in the description file of the simulation (default "")
- `time_start`:      initial simulation time (default 0.)
- `final_measures`:  measures to make at the end of simulation (default []) see `measure` and `output`
- `time_format`:     C like format for output of simulation time (default "%8.4g")
- `data_format`:     C like format for output of simulation data (default "%12.6g")
"""
@kwdef struct SimData
    description::String=""
    name::String
    time_start = 0.
    final_measures = []
    time_format::String = "%8.4g"
    data_format::String = "%14.8g"
    phases
end

show(io::IO, s::SimData) =
    print(io,
    """
    SimData(
        description = $(repr(s.description)),
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        time_format = $(repr(s.time_format)),
        data_format = $(repr(s.data_format)),
        phases =
    $(s.phases))"""
    )


"""
    runTMS(::SimData)
    runTMS(::SimData; clean = true)
    runTMS(::SimData; restart = true)
    runTMS(::SimData; output = myoutput)

run the given simulation (see SimData for details) and return a Simulation object containing the result.
clean remove the simulation directory and exit,
restart remove the simulation directory and run the simulation,
output redirect all output to the given IO channel (no output directory created), usefull values are stdout or devnull (to suppress all output).
"""
function runTMS(sim_data::SimData; restart::Bool=false, clean::Bool=false, output::Union{Nothing, IO} = nothing)
    live = isnothing(output)
    if live && (restart || clean)
        rm(sim_data.name; recursive = true, force = true)
    end
    if clean 
        return
    end
    start_dir = pwd()
    try
        if live
            mkpath(sim_data.name);
            cd(sim_data.name);
            touch("running")
            if sim_data.description ≠ ""
                write("description", sim_data.description)
            end
            write("stamp", """
                    Julia $VERSION
                    TensorMixedStates $(pkgversion(TensorMixedStates))
                    Date $(now())
                    """)
            src_path = Base.source_path()
            if src_path ≠ ""
                cp(src_path, "prog.jl"; force = true)
            else
                write_prog("prog.jl", sim_data)
            end
        end
        sim = Simulation(nothing; output, sim_data.time_format, sim_data.data_format)
        sim = log_phase(sim, sim_data)
        if live
            rm("running")
            cd(start_dir)
        end
        return sim
    catch
        if live
            touch("error")
            rm("running"; force = true)
            cd(start_dir)
        end
        rethrow()
    end
end

function log_phase(sim::Simulation, phases::Vector)
    for phase in phases
        sim = log_phase(sim, phase)
    end
    return sim
end

function log_phase(sim::Simulation, phase)
    log_msg(sim, "\n***** Starting phase \"$(phase.name)\" *****")
    if !isnothing(phase.time_start)
        sim = Simulation(sim, sim.state, phase.time_start)
    end
    td = @timed begin
        sim = run_phase(sim, phase)
        output(sim, phase.final_measures)
    end
    elapsed = round(td.time; digits=3)
    comp =
        if haskey(td, :compile_time)
            round(td.compile_time + td.recompile_time; digits=3)
        else
            nothing
        end
    log_msg(sim, "***** Ending phase \"$(phase.name)\" after $elapsed seconds, $(Base.format_bytes(td.bytes)) allocated *****")
    if !isnothing(comp)
        log_msg(sim, "compilation time was $comp seconds")
    end
    return sim
end

write_prog(filename, s::SimData) =
    write(filename,
        """
        using TensorMixedStates

        simdata = $s

        runTMS(simdata)
        """)
