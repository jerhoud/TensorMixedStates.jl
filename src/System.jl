import Base: length, show
export Pure, Mixed, TPure, TMixed, System, site


"""
    struct Pure end
    Pure()

    Correspond to pure quantum representation
"""
struct Pure end

"""
    Struct Mixed end
    Mixed()

    Correspond to mixed quantum representation
"""
struct Mixed end

TPure = Type{Pure}
TMixed = Type{Mixed}

site(::TPure, type::String; kwargs...) =
    siteind(type; kwargs...)

function site(::TMixed, type::String; kwargs...)
    i = siteind(type; kwargs...)
    return addtags(Index(2dim(i), tags(i)), "Mixed" * type)
end

"""
    site(::String; qns_info...)

Represent a quantum site of the given type

# Examples
    site("Qubit")
"""
site(type::String; kwargs...) =
    tp -> site(tp, type; kwargs...)

"""
    struct System

Represent a quantum system

# Fields
- `pure_sites::Vector{Index}`: Site indices for pure representation
- `mixed_sites::Vector{Index}`: Site indices for mixed representation

# Examples
    System(10, "Qubit")
    System(["Qubit", "SpinOne", "Qubit"])
"""
@kwdef struct System
    pure_sites::Vector{Index}
    mixed_sites::Vector{Index}
    create_string::String
end
    
System(sitetypes::Vector; create_string::String="") =
    System(map(s->s(Pure), sitetypes), map(s->s(Mixed), sitetypes), create_string)

System(size::Int, sitetype; create_string::String="") =
    System(fill(sitetype, size); create_string)

System(sitenames::Vector{String}; create_string::String="System($sitenames)") =
    System(map(site, sitenames); create_string)

System(size::Int, sitename::String; create_string::String="System($size, \"$sitename\")") =
    System(fill(sitename, size); create_string)

function show(io::IO, s::System)
    if s.create_string ≠ ""
        print(io, s.create_string)
    else
        print(io, "System($(s.pure_sites), $(s.mixed_sites))")
    end
end

"""
    length(::System)

Return number of sites in the system
"""
length(system::System) = length(system.pure_sites)

"""
    sites(Pure|Mixed, ::System)

Return the corresponding site index vector from the system
"""
sites(::TPure, system::System) = system.pure_sites
sites(::Any, system::System) = system.mixed_sites