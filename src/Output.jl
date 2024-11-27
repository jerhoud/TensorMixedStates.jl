struct Trace end
struct Trace2 end
struct Norm end
struct Purity end

@kwdef struct EE
    pos::Int
    spectre::Int = 0
end

@kwdef struct Check
    obs
    value
    tol = nothing
end

@kwdef struct Output
    outputs::Dict{IOStream, Vector{String}} = Dict()
    prodlits::Dict{String, ProdLit} = Dict()
    obs1::Dict{String, Any} = Dict()
    obs2::Dict{String, Any} = Dict()
    other::Dict{String, Any} = Dict()
    files::Dict{String, IOStream} = Dict()
end

prep_obs(obs::ProdLit) = string(obs) => insertFfactors(obs)
prep_obs(obs::Tuple) = string((obs[1](), obs[2]())) => obs
prep_obs(obs::Vector) = map(prep_obs, obs)
prep_obs(obs) = string(obs()) => obs

function Output(out::Pair{IOStream, Any}, output::Output = Output())

end

function Output(out::Pair{String, Any}, output::Output = Output())
    io = get!(output.files, out[1]) do filename
        return open(filename, "a")
    end
    Output(io => out[2], output)
end

Output(out::Vector, output::Output = Output()) = 
    Output(stdout => out, output)

Output(outs...) = Output([outs...])

function Output(outs::Vector{<:Pair}, output::Output = Output())
    foreach(outs) do out
        Output(o, output)
    end
    return output
end

