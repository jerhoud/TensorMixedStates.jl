export save_state, load_state

const state_file_version = 1

site_params(site::AbstractSite) =
    Float64[ getfield(site, f) for f in fieldnames(typeof(site)) ]

"""
    save_state(filename, statename, state)

save the state to disk in a hdf5 file,
several states with different names can be saved in the same file,
saving under a name already present in the file replaces it

# Examples

    save_state("myfile.h5", "ground_state", state)
"""
function save_state(filename::String, statename::String, state::State{R}) where R
    sites = state.system.sites
    h5open(filename, "cw") do f
        if haskey(f, statename)
            delete_object(f, statename)
        end
        g = create_group(f, statename)
        attributes(g)["version"] = state_file_version
        attributes(g)["type"] = string(nameof(R))
        g["modules"] = [ string(nameof(parentmodule(typeof(s)))) for s in sites ]
        g["types"] = [ string(nameof(typeof(s))) for s in sites ]
        g["nparams"] = [ length(fieldnames(typeof(s))) for s in sites ]
        g["params"] = reduce(vcat, map(site_params, sites); init = Float64[])
        g["state"] = state.state
    end
    return nothing
end

function site_module(name::String)
    if name == string(nameof(@__MODULE__))
        return @__MODULE__
    end
    for m in values(Base.loaded_modules)
        if string(nameof(m)) == name
            return m
        end
    end
    error("cannot find module $name needed to rebuild sites, is it loaded ?")
end

function build_site(modname::String, typename::String, params::Vector{Float64})
    t = getfield(site_module(modname), Symbol(typename))
    if !(t isa Type && t <: AbstractSite)
        error("$modname.$typename is not a site type")
    end
    return t(params...)
end

"""
    load_state(filename, statename)

load a state previously saved by `save_state`

The site types of the state are rebuilt by name, so the modules defining them
must be loaded (this is automatic for the site types of this package)

# Examples

    load_state("myfile.h5", "ground_state")
"""
function load_state(filename::String, statename::String)
    h5open(filename, "r") do f
        g = open_group(f, statename)
        version = read(attributes(g)["version"])
        if version ≠ state_file_version
            error("state \"$statename\" has file version $version, expected $state_file_version")
        end
        type = read(attributes(g)["type"])
        modules = read(g, "modules")
        types = read(g, "types")
        nparams = read(g, "nparams")
        params = read(g, "params")
        st = read(g, "state", MPS)
        if length(types) ≠ length(st)
            error("state \"$statename\" has $(length(types)) sites but a state of length $(length(st))")
        end
        sites = AbstractSite[]
        j = 0
        for (m, t, n) in zip(modules, types, nparams)
            push!(sites, build_site(m, t, params[j+1:j+n]))
            j += n
        end
        sites = identity.(sites)
        idx = Index[ siteind(st, i) for i in 1:length(st) ]
        if type == "Pure"
            return State{Pure}(System(sites, idx, map(mix, idx)), st)
        elseif type == "Mixed"
            return State{Mixed}(System(sites, map(Index, sites), idx), st)
        else
            error("state \"$statename\" has unknown type \"$type\"")
        end
    end
end
