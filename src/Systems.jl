export System, SysIndex

"""
    type System

represent a quantum system

# Fields

- `sites::Vector{<:AbstractSite}`: sites of the system
- `pure_indices::Vector{Index}`: Indices for pure representations
- `mixed_indices::Vector{Index}`: Indices for mixed representations

# Examples

    System(10, Qubit())
    System([Qubit(), Spin(1), Qubit(), Boson(5)])

# Indexation

    system[i]                  # gives site i
    SysIndex{Pure}(system, i)  # gives pure index i
    SysIndex{Mixed}(system, i) # gives mixed index i
"""
struct System
    sites::Vector{<:AbstractSite}
    pure_indices::Vector{Index}
    mixed_indices::Vector{Index}
end

"""
    is_charged(sites)

whether a system of those sites carries quantum numbers, that is whether any of them
declares something conserved. It is a property of the whole list: one site declaring a charge
makes every index of the system a charged one, the others taking a trivial charge, since an
MPS cannot mix the two kinds.
"""
is_charged(sites) = any(s -> !isempty(conserved(s)), sites)

# read off an index rather than walked over the sites again: it is asked once per operator
# placed on the system, which an MPO does for every factor of every term
is_charged(system::System) = hasqns(first(system.pure_indices))

function System(sites::Vector{<:AbstractSite})
    check_charges(sites)
    charged = is_charged(sites)
    pidx = [ site_index(s, charged) for s in sites ]
    midx = [ mixed_index(pidx[k], sites[k]) for k in eachindex(sites) ]
    return System(sites, pidx, midx)
end

System(size::Int, a::AbstractSite) = System(fill(a, size))

"""
    strong_names(::System)

the names of every quantity the sites of the system conserve strongly. See `strong`.
"""
strong_names(system::System) =
    unique(reduce(vcat, map(strong_names, system.sites); init = String[]))

symmetries(system::System) = symmetries(system.sites[1])

"""
    weaken(::System)
    weaken(::System, target)

the system a state lands on when `weaken(::State)` is given the same target: the same sites,
conserving less. Without a target, one level down.
"""
function weaken(system::System, target::Conserved)
    check_target(symmetries(system), target, "this system")
    return System([ weaken(s, target) for s in system.sites ])
end

weaken(system::System, spec) = weaken(system, Conserved(spec_names(spec)))

weaken(system::System) = weaken(system, one_step_down(symmetries(system)))

getindex(s::System, i...) = s.sites[i...]

show(io::IO, s::System) = print(io, "System($(s.sites))")

"""
    SysIndex{Pure|Mixed}(system, i)

returns the pure or mixed ITensor.Index for site i
"""
struct SysIndex{R <: PM}
    SysIndex{Pure}(s::System, i::Int...) = s.pure_indices[i...]
    SysIndex{Mixed}(s::System, i::Int...) = s.mixed_indices[i...]
    SysIndex{R}(s::System, is) where R = map(i->SysIndex{R}(s, i), is)
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
    if sys2 === sys1
        sys2 = sim(sys1)
    end
    return System(
        [sys1.sites ; sys2.sites],
        [sys1.pure_indices ; sys2.pure_indices],
        [sys1.mixed_indices ; sys2.mixed_indices])
end

tensor(sys1::System, sys2::System) = sys1 ⊗ sys2


"""
    check_indices(system, op)

check that every site an indexed operator acts on is a site of the system.

Nothing between writing `X(10)` and contracting its tensor compares that number with the
size of the system, and the three paths an indexed operator can take reach a different
array first: `expect_norm` indexes the mps, `PreMPO` its own link dimensions and `apply`
the sites. Each used to report a `BoundsError` on an internal vector the caller has no
reason to know. The check is made at those three entries instead, and names the factor at
fault rather than the array.
"""
function check_index(system::System, i::Int, a)
    n = length(system)
    if i < 1 || i > n
        error("$a acts on site $i, which the system does not have: it has $n sites")
    end
    return nothing
end

function check_indices(system::System, a::AtIndex)
    for i in a.index
        check_index(system, i, a)
    end
    return nothing
end

check_indices(system::System, a::Multi_F) =
    (check_index(system, a.start, a); check_index(system, a.stop, a))
check_indices(system::System, a::Union{SumOp, ProdOp}) =
    foreach(x -> check_indices(system, x), a.subs)
check_indices(system::System, a::ScalarOp) = check_indices(system, a.arg)
check_indices(system::System, a::Evolver) = check_indices(system, a.arg)
# a vector is a time dependent evolver, one term per coefficient
check_indices(system::System, a::Vector) = foreach(x -> check_indices(system, x), a)
# a generic operator carries no index, and neither does anything else that may be passed
check_indices(::System, _) = nothing

"""
    tensor(::System, ::AtIndex)

returns a tensor representing the given simple indexed operator acting on this system
"""
function tensor(system::System, a::AtIndex{R}) where R
    s = map(i->system[i], a.index)
    t = tensor(a.op, s...; charged = is_charged(system))
    is = SysIndex{R}(system, a.index)
    j = tensor_index(t)
    c = combinerto(j, reverse(is)...)
    # the primed combiner is daggered so that the two sides of the operator carry opposite
    # directions, which is what a charged index requires and what a dense one ignores
    return t * c * dag(c')
end

