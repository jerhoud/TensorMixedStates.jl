import Base: length, show
export Pure, Mixed, TPure, TMixed, Site, System, sites


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

"""
    struct Site

Represent a qunatum site. Applied on Pure or Mixed return a corresponding ITensor index

# Examples
    Site("Qubit")
    Site("Boson"; dim = 5)
"""
struct Site
    type::String
    kwargs
    Site(type::String; kwargs...) = new(type, kwargs)
end

show(io::IO, s::Site) =
    show_func(io, "Site", [s.type], s.kwargs)

(s::Site)(::TPure) = siteind(s.type; s.kwargs...)
function (s::Site)(::TMixed)
    i = s(Pure)
    c = combiner(i, i'; tags = tags(i))
    return addtags(combinedind(c), "Mixed" * s.type)
end    

"""
    struct System

Represent a quantum system

# Fields
- `pure_sites::Vector{Index}`: Site indices for pure representation
- `mixed_sites::Vector{Index}`: Site indices for mixed representation

# Examples
    System(10, "Qubit")
    System(["Qubit", "SpinOne", "Qubit"])
    System(10, Site("Boson"; dim = 5))
    System([Site(Qubit), Site("Boson"; dim = 5), Site("Qubit")])
"""
@kwdef struct System
    pure_sites::Vector{Index}
    mixed_sites::Vector{Index}
    show_data::String
end
    
System(sitetypes::Vector; show_data = "System($sitetypes)") =
    System(map(s->s(Pure), sitetypes), map(s->s(Mixed), sitetypes), show_data)

System(size::Int, sitetype) =
    System(fill(sitetype, size); show_data = "System($size, $sitetype)")

System(sitenames::Vector{String}) =
    System(map(Site, sitenames))

System(size::Int, sitename::String; kwargs...) =
    System(size, Site(sitename; kwargs...))

show(io::IO, s::System) =
    print(io, s.show_data)

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