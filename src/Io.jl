# save_state and load_state, which write a state to a file and read it back with its system, the
# sites being rebuilt from their recorded fields. A checkpoint reads the sites alone with
# `saved_sites`.

export save_state, load_state

"""
    state_file_version

the version of the format `save_state` writes
"""
const state_file_version = 3

"""
    readable_state_file_versions

the versions of the format `load_state` reads
"""
const readable_state_file_versions = (1, 2, 3)

"""
    param_kind(x)

the name a state file gives to the kind of the value `x` of a site field, or `nothing` when a
state file cannot carry it. What a site conserves is written apart, see `save_state`.
"""
param_kind(::Bool) = "Bool"
param_kind(::Integer) = "Int"
param_kind(::AbstractFloat) = "Float"
param_kind(::Symbol) = "Symbol"
param_kind(::AbstractString) = "String"
param_kind(::Nothing) = "Nothing"
param_kind(::Conservation) = "Conservation"
param_kind(_) = nothing

"""
    param_value(kind, s)

the value of a site field, read back from its kind and the string `save_state` wrote, an
integer too large for an `Int` as a `BigInt`; a conservation is read apart
"""
param_value(kind::AbstractString, s::AbstractString) =
    if kind == "Bool"
        parse(Bool, s)
    elseif kind == "Int"
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
    elseif kind == "Conservation"
        Conservation()
    else
        error("state file describes a site field as \"$kind\", which this version does not know")
    end

"""
    param_string(x)

the value `x` of a site field as a state file writes it, a float as the `Float64` it widens to
"""
param_string(x::AbstractFloat) = string(Float64(x))
param_string(::Conservation) = ""
param_string(x) = string(x)

"""
    site_params(site)

the fields of `site` as a state file carries them: the kinds of their values and the values
written as strings. A field no kind fits is refused.
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
                  "number, a boolean, a symbol, a string, nothing or a Conservation")
        end
        push!(kinds, k)
        push!(values, param_string(x))
    end
    return (kinds, values)
end

"""
    write_conservations(g, sites)

write in the group `g` what each site conserves, one row per quantity: its site, name,
modulus, strength, and its charges, those of every quantity following each other
"""
function write_conservations(g, sites)
    rows = [ (k, q) for (k, s) in enumerate(sites) for q in conserved(s).quantities ]
    if isempty(rows)
        return nothing
    end
    g["cons_site"] = [ k for (k, _) in rows ]
    g["cons_name"] = [ q.name for (_, q) in rows ]
    g["cons_modulus"] = [ q.modulus for (_, q) in rows ]
    g["cons_strong"] = [ Int(q.strong) for (_, q) in rows ]
    g["cons_ncharges"] = [ length(q.charges) for (_, q) in rows ]
    g["cons_charges"] = reduce(vcat, [ q.charges for (_, q) in rows ])
    return nothing
end

"""
    read_conservations(g, n)

what each of the `n` sites conserves, as `write_conservations` wrote it in the group `g`
"""
function read_conservations(g, n::Int)
    qs = [ Quantity[] for _ in 1:n ]
    if haskey(g, "cons_site")
        charges = read(g, "cons_charges")
        j = 0
        for (k, name, modulus, st, m) in zip(read(g, "cons_site"), read(g, "cons_name"),
                                             read(g, "cons_modulus"), read(g, "cons_strong"),
                                             read(g, "cons_ncharges"))
            push!(qs[k], Quantity(name, modulus, charges[j+1:j+m], st == 1))
            j += m
        end
    end
    return [ Conservation(q) for q in qs ]
end

"""
    decode_conserve(s)

what a site conserves, read from the string versions 1 and 2 of the format recorded it in:
quantities separated by `;`, each its name, `%` and its modulus when not 1, `!` when strong,
then `:` and the charges separated by `,`
"""
function decode_conserve(s::AbstractString)
    if isempty(s)
        return Conservation()
    end
    qs = map(split(s, ';')) do part
        i = findfirst(==(':'), part)
        if isnothing(i)
            error("a site records \"$part\" as a conserved quantity, which has no charges")
        end
        head, tail = part[1:i-1], part[i+1:end]
        st = endswith(head, '!')
        if st
            head = head[1:end-1]
        end
        j = findlast(==('%'), head)
        name, modulus = isnothing(j) ? (String(head), 1) :
                        (String(head[1:j-1]), parse(Int, head[j+1:end]))
        Quantity(name, modulus, parse.(Int, split(tail, ',')), st)
    end
    return Conservation(Vector{Quantity}(qs))
end

"""
    write_state(group, state)

write the tensors of `state` in the HDF5 group `group`, where `save_state` writes its type and
its sites: the MPS of a `State`. The state of an extension needs a method of its own, and one
of `read_state`.
"""
function write_state(g, state::State)
    g["state"] = state.state
    return nothing
end

"""
    state_type_name(state)

the type of `state` as a state file records it: `"Pure"` or `"Mixed"` for a `State`, the full
path of its type for the state of an extension
"""
state_type_name(::State{R}) where R = string(nameof(R))
state_type_name(state::AbstractState) =
    join((fullname(parentmodule(typeof(state)))..., nameof(typeof(state))), ".")

"""
    save_state(filename, statename, state)

save the state in the HDF5 file `filename` under the name `statename`, replacing a state saved
under that name. Every field of every site must be a number, a boolean, a symbol, a string,
`nothing` or a `Conservation`.

# Examples

    save_state("myfile.h5", "ground_state", state)
"""
function save_state(filename::String, statename::String, state::AbstractState)
    sites = state.system.sites
    # checked before the file is opened: a state that cannot be written must neither leave a
    # half written group nor delete the state already saved under its name
    ps = map(site_params, sites)
    if !hasmethod(write_state, Tuple{HDF5.Group, typeof(state)})
        error("cannot save a $(typeof(state)), its type has no method of " *
              "TensorMixedStates.write_state")
    end
    h5open(filename, "cw") do f
        if haskey(f, statename)
            delete_object(f, statename)
        end
        g = create_group(f, statename)
        attributes(g)["version"] = state_file_version
        attributes(g)["type"] = state_type_name(state)
        # the whole path, so that a site type of a module inside another one is found again
        g["modules"] = [ join(fullname(parentmodule(typeof(s))), ".") for s in sites ]
        g["types"] = [ string(nameof(typeof(s))) for s in sites ]
        g["nparams"] = [ length(fieldnames(typeof(s))) for s in sites ]
        g["pkinds"] = reduce(vcat, first.(ps); init = String[])
        g["params"] = reduce(vcat, last.(ps); init = String[])
        write_conservations(g, sites)
        write_state(g, state)
    end
    return nothing
end

"""
    site_module(name)

the module a state file names for a site type, or for the type of the state of an extension:
the path from a root module down, `Main.MySites` for a module defined in a script
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
        error("cannot find module $name needed to read the state, is it loaded ?")
    end
    return m
end

"""
    field_value(t, x)

the value `x` read for a field of type `t`: a `Conservation` from the string versions 1 and 2
wrote, any other value converted to the type of the field
"""
field_value(::Type{Conservation}, x::AbstractString) = decode_conserve(x)
field_value(t::Type, x) = convert(t, x)

"""
    build_site(modname, typename, params)

the site of type `typename` of the module `modname`, built from the values `params` of its
fields
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
    return rebuild_site(t, [ field_value(ft[i], p) for (i, p) in enumerate(params) ],
                        "a state file to be read")
end

"""
    read_sites(g, version)

the sites a state file records in the group `g`, as version `version` of the format writes them
"""
function read_sites(g, version)
    modules = read(g, "modules")
    types = read(g, "types")
    nparams = read(g, "nparams")
    # version 1 wrote every field as a Float64
    params = version == 1 ? collect(read(g, "params")) :
                            map(param_value, read(g, "pkinds"), read(g, "params"))
    cons = version ≥ 3 ? read_conservations(g, length(types)) : nothing
    sites = AbstractSite[]
    j = 0
    for (k, (m, t, n)) in enumerate(zip(modules, types, nparams))
        ps = params[j+1:j+n]
        if !isnothing(cons)
            ps = [ p isa Conservation ? cons[k] : p for p in ps ]
        end
        push!(sites, build_site(m, t, ps))
        j += n
    end
    return identity.(sites)
end

"""
    read_mps(g, system, R)

the MPS a state file holds in the group `g`, on the indices of `system` in the representation
`R`. Version 2 kept charged indices on a partial trace whose sites conserve nothing, of a single
sector each, which densifying puts in the order of the plain ones.
"""
function read_mps(g, system::System, ::Type{R}) where R
    st = read(g, "state", MPS)
    if length(system) ≠ length(st)
        error("state \"$(lstrip(HDF5.name(g), '/'))\" has $(length(system)) sites but a state " *
              "of length $(length(st))")
    end
    if hasqns(st[1]) && !is_charged(system)
        st = dense(st)
    end
    return replace_siteinds(st, siteinds(R, system))
end

"""
    read_state(::Type{S}, group, sites)

the state of type `S` that `save_state` wrote in the HDF5 group `group`, on the sites `sites`.
The state of an extension needs a method of its own, reading what its `write_state` wrote.
"""
function read_state(::Type{State{R}}, g, sites) where R
    system = System(sites)
    return State{R}(system, read_mps(g, system, R))
end

read_state(t::Type{<:AbstractState}, _, _) =
    error("cannot read a $t, its type has no method of TensorMixedStates.read_state")

"""
    state_type(name, statename)

the type of the state `statename`, which its file names `name`
"""
function state_type(name::String, statename::String)
    if name == "Pure"
        return State{Pure}
    elseif name == "Mixed"
        return State{Mixed}
    end
    i = findlast(==('.'), name)
    t = nothing
    if !isnothing(i)
        m = site_module(name[1:prevind(name, i)])
        s = Symbol(name[nextind(name, i):end])
        t = isdefined(m, s) ? getfield(m, s) : nothing
    end
    if !(t isa Type && t <: AbstractState)
        error("state \"$statename\" has unknown type \"$name\"")
    end
    return t
end

"""
    file_version(g, statename)

the version of the format of the state `statename`, refused when it is not readable
"""
function file_version(g, statename)
    version = read(attributes(g)["version"])
    if !(version in readable_state_file_versions)
        error("state \"$statename\" has file version $version, expected one of " *
              join(readable_state_file_versions, ", "))
    end
    return version
end

"""
    saved_sites(filename, statename)

the sites of the state saved under the name `statename` in the file `filename`
"""
saved_sites(filename::String, statename::String) =
    h5open(filename, "r") do f
        g = open_group(f, statename)
        return read_sites(g, file_version(g, statename))
    end

"""
    load_state(filename, statename)

the state saved under the name `statename` in the file `filename` by `save_state`, on the system
of its sites. The modules defining its site types, and the state type of an extension, must be
loaded.

# Examples

    load_state("myfile.h5", "ground_state")
"""
load_state(filename::String, statename::String) =
    h5open(filename, "r") do f
        g = open_group(f, statename)
        version = file_version(g, statename)
        type = state_type(read(attributes(g)["type"]), statename)
        return read_state(type, g, read_sites(g, version))
    end
