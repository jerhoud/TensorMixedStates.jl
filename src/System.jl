import Base: length
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
    site(::String)

    Represent a quantum site of the given type
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
end
    
System(sitetypes::Vector) =
    System(map(s->s(Pure), sitetypes), map(s->s(Mixed), sitetypes))

System(size::Int, sitetype) =
    System(fill(sitetype, size))

System(sitenames::Vector{String}) =
    System(map(site, sitenames))

System(size::Int, sitename::String) =
    System(fill(sitename, size))

length(system::System) = length(system.pure_sites)
