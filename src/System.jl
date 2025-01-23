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

function show_kwargs(kwargs)
    s = ""
    for (sym, val) in pairs(kwargs)
        s *= " $sym = $val,"
    end
    return s
end

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
    if isempty(s.kwargs)
        print(io, "Site($(repr(s.type)))")
    else
        print(io, "Site($(repr(s.type));$(show_kwargs(s.kwargs)))")
    end

(s::Site)(::TPure) = siteind(type; kwargs...)
function (s::Site)(::TMixed)
    i = s(Pure)
    return addtags(Index(2dim(i), tags(i)), "Mixed" * type)
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
    show_data
end
    
System(sitetypes::Vector; show_data = sitetypes) =
    System(map(s->s(Pure), sitetypes), map(s->s(Mixed), sitetypes), show_data)

System(size::Int, sitetype) =
    System(fill(sitetype, size); show_data = size => sitetype)

System(sitenames::Vector{String}) =
    System(map(Site, sitenames))

System(size::Int, sitename::String; kwargs...) =
    System(size, Site(sitename; kwargs...))

function show(io::IO, s::System)
    if s.show_data isa Pair
        print(io, "System($(first(s.show_data)), $(repr(last(s.show_data))))")
    else
        print(io, "System($(s.show_data))")
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