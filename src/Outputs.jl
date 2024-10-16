using Printf

function write_output(output::Output)
    tfmt = Printf.Format(output.time_format)
    dfmt = Printf.Format(output.data_format)
    if output.file ==""
        flog = files["log"]
        fdata = files["data"]
        fexpect = files["expect"]
        fcorrel = files["correl"]
    else
        flog = fdata = fexpect = fcorrel = open(output.file, "a")
    end
    prep = preprocess(sim_type, sim_state)
    if (output.state_info)
        write_state_info(flog, tfmt, dfmt)
    end
    write_observables(fdata, output.observables, prep, tfmt, dfmt)
    write_expect(fexpect, output.expect, prep, tfmt, dfmt)
    write_correl(fcorrel, output.correl, prep, tfmt, dfmt)
    if output.file ≠ ""
        close(flog)
    end
end

function write_state_info(file::IO, tfmt::Printf.Format, dfmt::Printf.Format)
    p = sim_state
    println(file)
    print(file, "simulation time: ")
    Printf.format(file, tfmt, sim_time)
    println(file)
    if sim_type == Pure
        println(file, "Pure state with $(length(p)) sites")
    else
        println(file, "Mixed state with $(length(p)) sites")
        print(file, "Trace is ")
        Printf.format(file, dfmt, trace(sim_type, sim_state))
        print(file, " and squared trace is ")
        Printf.format(file, dfmt, trace2(sim_type, sim_state))
        println(file)
    end
    println(file, "maxdim is $(maxlinkdim(p)) and memory usage $(Base.format_bytes(Base.summarysize(p)))")
end

function write_data(file::IO, header::Union{String, ProdLit}, data, tfmt::Printf.Format, dfmt::Printf.Format)
    print(file, header, "\t")
    Printf.format(file, tfmt, sim_time)
    for x in data
        print(file, "\t")
        Printf.format(file, dfmt, real(x))
        if imag(x) ≠ 0
            print(file, "+i*")
            Printf.format(file, dfmt, imag(x))
        end
    end 
    println(file)
end

function write_data(file::IO, header, data, tfmt, dfmt)
    for (h, d) in zip(header, data)
        write_data(file, h, d, tfmt, dfmt)
    end
end

function write_expect(file::IO, op::Vector{String}, prep, tfmt::Printf.Format, dfmt::Printf.Format)
    if op ≠ []
        write_data(file, op, expect(sim_type, sim_state, op, prep), tfmt, dfmt)
        flush(file)
    end
end

function write_correl(file::IO, op::Vector{Tuple{String, String}}, prep, tfmt::Printf.Format, dfmt::Printf.Format)
    if op ≠ []
        write_data(file, map(t->t[1]*t[2],op), correlations(sim_type, sim_state, op, prep), tfmt, dfmt)
        flush(file)
    end
end

function write_observables(file::IO, obs::Vector{ProdLit}, prep, tfmt::Printf.Format, dfmt::Printf.Format)
    if op ≠ []
        write_data(file, obs, expect(sim_type, sim_state, obs, prep), tfmt, dfmt)
        flush(file)
    end
end

