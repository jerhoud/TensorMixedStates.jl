export save_state, load_state

const state_file_version = 2
const readable_state_file_versions = (1, 2)

"""
    param_kind(x)

the name a state file gives to the kind of value a site field holds, or `nothing` when a
state file cannot carry it.

Version 1 of the format wrote every field as a `Float64`, so a site carrying a `Symbol`, a
name or a flag could not be saved at all: `save_state` failed on a `convert` raised deep
inside HDF5, naming neither the site nor the field. Since a site has to be able to declare
what it conserves, the fields now travel as a string each, together with the name of what
they are, which also means reading them back does not depend on the field types of the site
being declared concretely.

`Bool` comes first on purpose, being an `Integer` as far as dispatch is concerned.
"""
param_kind(::Bool) = "Bool"
param_kind(::Integer) = "Int"
param_kind(::AbstractFloat) = "Float"
param_kind(::Symbol) = "Symbol"
param_kind(::AbstractString) = "String"
param_kind(::Nothing) = "Nothing"
param_kind(_) = nothing

"""
    param_value(kind, s)

the value a site field had, read back from the kind and the string `save_state` wrote
"""
param_value(kind::AbstractString, s::AbstractString) =
    if kind == "Bool"
        parse(Bool, s)
    elseif kind == "Int"
        parse(Int, s)
    elseif kind == "Float"
        parse(Float64, s)
    elseif kind == "Symbol"
        Symbol(s)
    elseif kind == "String"
        String(s)
    elseif kind == "Nothing"
        nothing
    else
        error("state file describes a site field as \"$kind\", which this version does not know")
    end

"""
    site_params(site)

the fields of a site, as the kinds of value they hold and those values written as strings,
which is what a state file carries. See `param_kind`.
"""
function site_params(site::AbstractSite)
    kinds = String[]
    values = String[]
    for f in fieldnames(typeof(site))
        x = getfield(site, f)
        k = param_kind(x)
        if isnothing(k)
            error("cannot save a state on site $(typeof(site)): its field $f is a " *
                  "$(typeof(x)), which a state file cannot carry. A site field must be a " *
                  "number, a boolean, a symbol, a string or nothing")
        end
        push!(kinds, k)
        push!(values, string(x))
    end
    return (kinds, values)
end

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
    # read before the file is opened. A site field a state file cannot carry must not leave a
    # half written group behind, and above all must not reach `delete_object` first: saving a
    # state that cannot be written over a name already in the file would then destroy what was
    # there and put nothing in its place
    ps = map(site_params, sites)
    h5open(filename, "cw") do f
        if haskey(f, statename)
            delete_object(f, statename)
        end
        g = create_group(f, statename)
        attributes(g)["version"] = state_file_version
        attributes(g)["type"] = string(nameof(R))
        # the whole path, so that a site type of a module inside another one is found again
        g["modules"] = [ join(fullname(parentmodule(typeof(s))), ".") for s in sites ]
        g["types"] = [ string(nameof(typeof(s))) for s in sites ]
        g["nparams"] = [ length(fieldnames(typeof(s))) for s in sites ]
        g["pkinds"] = reduce(vcat, first.(ps); init = String[])
        g["params"] = reduce(vcat, last.(ps); init = String[])
        g["state"] = state.state
    end
    return nothing
end

"""
    site_module(name)

the module a state file names for a site type: the path from a root module down, as
`save_state` writes it, `Main.MySites` for a module defined in a script. Files written before
carried the last name only, which is the whole path of a root module, so they read as they did.
"""
function site_module(name::String)
    root, path... = split(name, '.')
    m = nothing
    if root == string(nameof(@__MODULE__))
        m = @__MODULE__
    else
        for r in values(Base.loaded_modules)
            if string(nameof(r)) == root
                m = r
                break
            end
        end
    end
    for p in path
        if !(m isa Module && isdefined(m, Symbol(p)))
            m = nothing
            break
        end
        m = getfield(m, Symbol(p))
    end
    if !(m isa Module)
        error("cannot find module $name needed to rebuild sites, is it loaded ?")
    end
    return m
end

function build_site(modname::String, typename::String, params::Vector)
    t = getfield(site_module(modname), Symbol(typename))
    if !(t isa Type && t <: AbstractSite)
        error("$modname.$typename is not a site type")
    end
    ft = fieldtypes(t)
    if length(params) > length(ft)
        error("state file gives $(length(params)) parameters for site $typename, which has " *
              "$(length(ft)) fields")
    end
    # version 1 wrote every field as a `Float64`, so the reader converts back to what the
    # site declares. Widening the constructors to take a `Real` dimension would put the
    # conversion in the wrong place and state something looser than the truth
    return t([ convert(ft[i], p) for (i, p) in enumerate(params) ]...)
end

"""
    load_state(filename, statename[; system])

load a state previously saved by `save_state`

The site types of the state are rebuilt by name, so the modules defining them
must be loaded (this is automatic for the site types of this package)

`system` reads the state onto an existing `System` rather than onto one built from the
file. The sites must match, and this is what makes the state comparable with one already
in hand, `inner` and the fidelities requiring their arguments to share a system.

# Examples

    load_state("myfile.h5", "ground_state")
    load_state("myfile.h5", "ground_state"; system = sim.state.system)
"""
function load_state(filename::String, statename::String;
                    system::Union{Nothing, System} = nothing)
    st = h5open(filename, "r") do f
        g = open_group(f, statename)
        version = read(attributes(g)["version"])
        if !(version in readable_state_file_versions)
            error("state \"$statename\" has file version $version, expected one of " *
                  join(readable_state_file_versions, ", "))
        end
        type = read(attributes(g)["type"])
        modules = read(g, "modules")
        types = read(g, "types")
        nparams = read(g, "nparams")
        # version 1 wrote every field as a Float64, which is what a checkpoint or a state
        # saved by an earlier version still holds
        if version == 1
            params = collect(read(g, "params"))
        else
            params = map(param_value, read(g, "pkinds"), read(g, "params"))
        end
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
            return State{Pure}(System(sites,
                idx, [ mixed_index(idx[k], sites[k]) for k in eachindex(sites) ]), st)
        elseif type == "Mixed"
            # the pure indices are rebuilt rather than read, only the mixed ones being in the
            # file, so they need the mode of the system the sites make up
            charged = is_charged(sites)
            return State{Mixed}(
                System(sites, [ site_index(s, charged) for s in sites ], idx), st)
        else
            error("state \"$statename\" has unknown type \"$type\"")
        end
    end
    # without a system the state comes back on one built from the file, whose indices are
    # its own. `system` puts it on an existing one instead, which is what comparing it with
    # a state already in hand requires
    if isnothing(system)
        return st
    else
        return State(system, st)
    end
end
