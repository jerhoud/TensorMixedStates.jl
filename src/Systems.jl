export System, Pure, Mixed

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

PM = Union{Pure, Mixed}

"""
    struct System

Represent a quantum system

# Fields
- `sites::Vector{<:AbstractSite}`: sites of the system
- `pure_indices::Vector{Index}`: Indices for pure representations
- `mixed_indices::Vector{Index}`: Indices for mixed representations

# Examples
    System(10, Qubit())
    System([Qubit(), SpinOne(), Qubit(), Boson(5)])

# Indexation
    system[i]          gives site i
    system[Pure(), i]  gives pure index i
    system[Mixed(), i] gives mixed index i
"""
struct System
    sites::Vector{<:AbstractSite}
    pure_indices::Vector{Index}
    mixed_indices::Vector{Index}
end

function System(sites::Vector{<:AbstractSite})
    pi = map(Index, sites)
    mi = map(mix, pi)
    return System(sites, pi, mi)
end

System(size::Int, a::AbstractSite) = System(fill(a, size))

getindex(s::System, i...) = s.sites[i...]
getindex(s::System, ::Pure, i...) = s.pure_indices[i...]
getindex(s::System, ::Mixed, i...) = s.mixed_indices[i...]


"""
    length(::System)

Return number of sites in the system
"""
length(system::System) = length(system.sites)

function tensor(type::PM, system::System, a::Indexed)
    s = map(i->system[i], a.index)
    t = tensor(a.op, s...)
    is = map(i->system[type, i], a.index)
    j = tensor_index(t)
    c = combinerto(j, reverse(is)...)
    return t * c * c'
end

