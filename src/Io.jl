# save_state and load_state, which write a state to a file and read it back with its system, the
# sites being rebuilt from their recorded parameters.

export save_state, load_state

"""
    state_file_version

the version of the format `save_state` writes.
"""
const state_file_version = 2

"""
    readable_state_file_versions

the versions of the format `load_state` reads.
"""
const readable_state_file_versions = (1, 2)

"""
    param_kind(x)

the name a state file gives to the kind of the value `x` of a site field, or `nothing` when a
state file cannot carry it. Each field is saved as a string along with this name, so that a
site may hold a `Symbol`, a name or a flag, which version 1 of the format, writing every field
as a `Float64`, could not, and so that reading it back does not depend on the field types of
the site being concrete. `Bool`, an `Integer` for dispatch, has a kind of its own.
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

the value of a site field, read back from its kind and the string `save_state` wrote. An
integer too large for an `Int` is read as a `BigInt`.
"""
param_value(kind::AbstractString, s::AbstractString) =
    if kind == "Bool"
        parse(Bool, s)
    elseif kind == "Int"
        # an integer beyond `Int` is read as the integer it is, and converted to the type of
        # its field by whoever rebuilds the site
        let n = tryparse(Int, s)
            isnothing(n) ? parse(BigInt, s) : n
        end
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
    param_string(x)

the value `x` of a site field as a state file writes it. A float is written as the `Float64`
it widens to exactly, and converted back to the type of its field on reading: `string(0.1f0)`
is `"0.1f0"`, which `parse(Float64, s)` refuses.
"""
param_string(x::AbstractFloat) = string(Float64(x))
param_string(x) = string(x)

"""
    site_params(site)

the fields of `site` as a state file carries them: the kinds of their values, see
`param_kind`, and the values written as strings. A field no kind fits is refused.
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
        push!(values, param_string(x))
    end
    return (kinds, values)
end

"""
    save_state(filename, statename, state)

save the state in the HDF5 file `filename` under the name `statename`. A file can hold several
states under different names, and saving under a name already present replaces that state.
Every field of every site must be an integer, a float, a boolean, a symbol, a string or
`nothing`.

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
`save_state` writes it, `Main.MySites` for a module defined in a script. Older files hold the
last name only, which is the whole path of a root module.
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

"""
    build_site(modname, typename, params)

the site of type `typename` of the module `modname`, built from the values `params` of its
fields, each converted to the type of its field, which a file of version 1, holding every field
as a `Float64`, requires.
"""
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

the state saved under the name `statename` in the file `filename` by `save_state`.

The site types are rebuilt by name, so the modules defining them must be loaded, which is
automatic for those of this package. The state comes back on a system built from the file,
or, when `system` is given, on that one, whose sites must match: this is what makes it
comparable with a state already in hand, `inner` and the fidelities requiring their arguments
to share a system.

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
            # file, so they take the mode of the stored ones: a partial trace of a charged
            # system keeps charged indices on sites that conserve nothing
            charged = hasqns(first(idx))
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
