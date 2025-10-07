export System, sim

"""
    type System

represent a quantum system

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

show(io::IO, s::System) = print(io, "System($(s.sites))")

struct SysIndex{R <: PM}
    SysIndex{Pure}(s::System, i::Int...) = s.pure_indices[i...]
    SysIndex{Mixed}(s::System, i::Int...) = s.mixed_indices[i...]
    SysIndex{R}(s::System, is) where R = map(i->SysIndex{R}(system, i), is)
 end

"""
    length(::System)

return the number of sites in the system
"""
length(system::System) = length(system.sites)

"""
    sim(::System)

create a clone of the system: identical but with different indices
"""
sim(system::System) =
    System(system.sites, sim.(system.pure_indices), sim.(system.mixed_indices))

"""
    ::System ⊗ ::System
    tensor(::System, ::System)

create the tensorial product of two systems
"""
function (sys1::System ⊗ sys2::System)
    if sys2 == sys1
        sys2 = sim(sys1)
    end
    return System(
        [sys1.sites ; sys2.sites],
        [sys1.pure_indices ; sys2.pure_indices],
        [sys1.mixed_indices ; sys2.mixed_indices])
end

tensor(sys1::System, sys2::System) = sys1 ⊗ sys2


"""
    tensor(::System, ::Indexed)

returns a tensor representing the given base indexed operator acting on this system
"""
function tensor(system::System, a::AtIndex{R}) where R
    s = map(i->system[i], a.index)
    t = tensor(a.op, s...)
    is = SysIndex{R}(system, a.index)
    j = tensor_index(t)
    c = combinerto(j, reverse(is)...)
    return t * c * c'
end

