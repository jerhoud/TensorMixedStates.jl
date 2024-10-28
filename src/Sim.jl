struct Die <: Exception end

sim_time::Float64 = 0
sim_type::Union{Type{Pure}, Type{Mixed}} = Pure
sim_state::MPS = MPS()
files::Dict{String, IO} = Dict()
start_dir::String = ""

export runNL

function runNL(sim_data::SimData)
  live = !sim_data.debug
  global start_dir = pwd()
  global files = Dict()
  global sim_type = Pure
  global sim_time = 0
  global sim_state = MPS()
  try
    if live && (sim_data.restart || sim_data.clean)
      rm(sim_data.name; recursive = true)
    end
    if sim_data.clean 
      return
    end
    if live
      mkpath(sim_data.name);
      cd(sim_data.name);
      touch("running")
      files["log"] = open("log", "a")
      files["data"] = open("data", "a")
      files["expect"] = open("expect", "a")
      files["correl"] = open("correl", "a")
    else
      println("start run")
      files["log"] = stdout
      files["data"] = stdout
      files["expect"] = stdout
      files["correl"] = stdout
    end

    start_phase("Run $(sim_data.name)")

    time_data = @timed begin
      if (sim_data.description ≠ "")
        if live
          write("description", sim_data.description)
        else
          println("description:\n", sim_data.description)
        end
      end

      run_phases(sim_data.phases)
    end

    end_phase("Run $(sim_data.name)", time_data)
  catch e
    if live
      touch("error")
      rm("running")
      cd(start_dir)
    else
      println("error")
    end
    if !isa(e, Die)
      rethrow()
    end
  else
    if live
      log("Closing data files and removing empty data files")
      foreach(close, values(files))
      foreach(keys(files)) do filename
        if filesize(filename) == 0
          rm(filename)
        end
      end
      rm("running")
      cd(start_dir)
    else
      println("end run")
    end
  end
end

function die(msg::String)
  log_file = files["log"]
  println(log_file, msg, "\nAborting...")
  close(log_file)
  throw(Die)
end

function log(msg::String)
  log_file = files["log"]
  println(log_file, msg)
  flush(log_file)
end

function start_phase(name::String)
  log_file = files["log"]
  println(log_file, "\n\n***** Starting phase \"$name\" *****")
  flush(log_file)
end

function end_phase(name::String, time_data)
  _, elapsed, bytes = time_data
  elapsed = round(elapsed; digits=3)
  log_file = files["log"]
  println(log_file, "\n***** Ending phase \"$name\" after $elapsed seconds, $(Base.format_bytes(bytes)) allocated *****")
  flush(log_file)
end
