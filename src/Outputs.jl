using Printf

function write_output(output::Output)
    tfmt = Printf.Format(output.time_format)
    dfmt = Printf.Format(output.data_format)
    if output.file ==""
        flog = files["log"]
        fdata = files["data"]
        fexpect = files["expect"]
        fcorrel = files["correl"]
        fee = files["entanglement_entropy"]
    else
        flog = fdata = fexpect = fcorrel = fee = open(output.file, "a")
    end
    prep = preprocess(sim_type, sim_state)
    if (output.state_info)
        write_state_info(flog, tfmt, dfmt)
    end
    write_observables(fdata, output.observables, prep, tfmt, dfmt)
    write_expect(fexpect, output.expect, prep, tfmt, dfmt)
    write_correl(fcorrel, output.correl, prep, tfmt, dfmt)
    write_ee(fee, output.entanglement_entropy, output.show_ee_spectrum, tfmt, dfmt)
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
        println(file, "Pure state with $(length(p)) sites and norm $(norm(p))")
    else
        println(file, "Mixed state with $(length(p)) sites")
        print(file, "Trace deviation to 1 is ")
        t = abs(trace(sim_type, sim_state) - 1)
        Printf.format(file, tfmt, t)
        print(file, " and effective number of states is ")
        Printf.format(file, tfmt, 1 / trace2(sim_type, sim_state))
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
        if abs(imag(x) > 1e-8)
            log_msg("imaginary part of data $header at time $sim_time is large ($(imag(x)))")
        end
    end 
    println(file)
end

function write_data(file::IO, header::Union{String, ProdLit}, data::Matrix, tfmt::Printf.Format, dfmt::Printf.Format)
    print(file, header, "\t")
    Printf.format(file, tfmt, sim_time)
    println(file)
    for l in 1:size(data, 1)
        write_data(file, "$header:$l", data[l,:], tfmt, dfmt)
    end
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

function write_ee(file::IO, ee, ss, tfmt, dfmt)
    if ee ≠ []
        for p in ee
            e, s = entanglement_entropy(sim_state, p)
            write_data(file, "EE-$p", e, tfmt, dfmt)
            if ss > 0
                if length(s) > ss
                    s = s[1:ss]
                end
                write_data(file, "SP-$p", s, tfmt, dfmt)
            end
        end
        flush(file)
    end
end
