# Systems, the sites a state lives on, with the indices of each site in the pure and the mixed
# representations, and the weakening of what a system conserves.

export System, SysIndex

"""
    System(sites::Vector{<:AbstractSite})
    System(n::Int, site::AbstractSite)

a quantum system: its sites, and the ITensor index of each site for pure and for mixed
representations. `System(n, site)` is made of `n` copies of `site`.

When any site conserves something, every index carries charges, a site conserving nothing
taking a trivial one. Sites whose conserved quantities cannot live together on one system are
refused.

# Fields

- `sites::Vector{<:AbstractSite}`: the sites
- `pure_indices::Vector{Index}`: the indices for pure representations
- `mixed_indices::Vector{Index}`: the indices for mixed representations

# Examples

    System(10, Qubit())
    System([Qubit(), Spin(1), Qubit(), Boson(5)])

# Indexation

    system[i]                  # site i
    SysIndex{Pure}(system, i)  # pure index of site i
    SysIndex{Mixed}(system, i) # mixed index of site i
"""
struct System
    sites::Vector{<:AbstractSite}
    pure_indices::Vector{Index}
    mixed_indices::Vector{Index}
end

"""
    is_charged(sites)
    is_charged(system)

whether a system of those sites carries quantum numbers, that is whether any of them conserves
something. One such site makes every index of the system charged, the others taking a trivial
charge, since an MPS cannot mix the two kinds.
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
    strong_names(sites)
    strong_names(::System)

the names of the quantities the sites, or the sites of the system, conserve strongly, see
`strong`
"""
strong_names(sites::AbstractVector{<:AbstractSite}) =
    unique(reduce(vcat, map(strong_names, sites); init = String[]))
strong_names(system::System) = strong_names(system.sites)

symmetries(system::System) =
    Conserved(unique(reduce(vcat, [ symmetries(s).names for s in system.sites ];
                            init = Tuple{String, Bool}[])))

"""
    weaken(::System)
    weaken(::System, target)

the same sites, each keeping of what it conserves only what `target` names, and as strongly as
`target` asks. `target` is written as `conserve` is given, or as `symmetries` returns it.
Without a target, the system goes one level down: every strong quantity becomes weak or, when
none is strong, every quantity is dropped.

A quantity can be dropped or made weak, but a target asking for a quantity the system does not
conserve, or for a weak one strongly, is refused. A target equal to what the system conserves
gives it back unchanged. This is the system `weaken(::State)` puts a state on.

# Examples

    s = System(4, Electron(conserve = (strong(Ntot), 2Sz)))
    weaken(s)          # conserves (2Sz, Ntot)
    weaken(s, 2Sz)     # conserves 2Sz only
    weaken(s, ())      # conserves nothing
"""
function weaken(system::System, target::Conserved)
    check_target(symmetries(system), target, "this system")
    if target.names == symmetries(system).names
        return system
    end
    return System([ weaken(s, target) for s in system.sites ])
end

weaken(system::System, spec) = weaken(system, Conserved(spec_names(spec)))

weaken(system::System) = weaken(system, one_step_down(symmetries(system)))

getindex(s::System, i...) = s.sites[i...]

show(io::IO, s::System) = print(io, "System($(s.sites))")

"""
    SysIndex{Pure}(system, i)
    SysIndex{Mixed}(system, i)

the ITensor index of site `i` of `system`, for pure or for mixed representations. Given a
collection of sites, as `1:n`, it gives their indices.

# Examples

    SysIndex{Pure}(system, 1)
    SysIndex{Mixed}(system, 1:length(system))
"""
struct SysIndex{R <: PM}
    SysIndex{Pure}(s::System, i::Int...) = s.pure_indices[i...]
    SysIndex{Mixed}(s::System, i::Int...) = s.mixed_indices[i...]
    SysIndex{R}(s::System, is) where R = map(i->SysIndex{R}(s, i), is)
 end

"""
    length(::System)

the number of sites of the system.
"""
length(system::System) = length(system.sites)

"""
    sim(::System)

a copy of the system with the same sites and new indices
"""
sim(system::System) =
    System(system.sites, sim.(system.pure_indices), sim.(system.mixed_indices))

"""
    ::System ⊗ ::System
    tensor(::System, ::System)

the tensor product of two systems: the sites of the first followed by those of the second,
with their indices. When the two share an index, as in `S ⊗ S`, the second is given new ones.
Both must carry charges or neither: the product of a charged system and an uncharged one is
built from its sites with `System`. Sites whose conserved quantities cannot live together on
one system are refused.
"""
function (sys1::System ⊗ sys2::System)
    if is_charged(sys1) ≠ is_charged(sys2)
        error("cannot take the tensor product of a system carrying charges and one carrying " *
              "none, build it from its sites with System instead")
    end
    check_charges([sys1.sites; sys2.sites])
    if !isdisjoint(sys1.pure_indices, sys2.pure_indices) ||
       !isdisjoint(sys1.mixed_indices, sys2.mixed_indices)
        sys2 = sim(sys2)
    end
    return System(
        [sys1.sites ; sys2.sites],
        [sys1.pure_indices ; sys2.pure_indices],
        [sys1.mixed_indices ; sys2.mixed_indices])
end

tensor(sys1::System, sys2::System) = sys1 ⊗ sys2


"""
    check_index(system, i, a)

refuse the site `i` of the operator `a` when `system` has no such site
"""
function check_index(system::System, i::Int, a)
    n = length(system)
    if i < 1 || i > n
        error("$a acts on site $i, which the system does not have: it has $n sites")
    end
    return nothing
end

"""
    check_indices(system, op)

check that every site the indexed operator `op` acts on is a site of `system`, naming the
factor at fault. Nothing else compares `X(10)` with the size of the system: without this check,
made where an operator enters `expect_norm`, `PreMPO` or `apply`, each would fail on a
`BoundsError` of an internal array.
"""
function check_indices(system::System, a::AtIndex)
    for i in a.index
        check_index(system, i, a)
    end
    return nothing
end

check_indices(system::System, a::Multi_F) =
    (check_index(system, a.start, a); check_index(system, a.stop, a))
check_indices(system::System, a::ComOp) =
    (check_index(system, first(com_sites(a)), a); check_index(system, last(com_sites(a)), a))
check_indices(system::System, a::Union{SumOp, ProdOp}) =
    foreach(x -> check_indices(system, x), a.subs)
check_indices(system::System, a::ScalarOp) = check_indices(system, a.arg)
check_indices(system::System, a::Evolver) = check_indices(system, a.arg)
check_indices(system::System, a::Vector) = foreach(x -> check_indices(system, x), a)
check_indices(::System, _) = nothing

"""
    check_one_site(a, what)

refuse a factor acting on several sites at once, which `what`, an MPO or `expect`, cannot
place, since both place one site factors only. Such a factor is what `simplify` leaves of an
operator of several sites with no expression to be replaced by, one defined by a matrix or a
function and created without its sites.
"""
function check_one_site(a, what)
    if a isa AtIndex && length(a.index) > 1
        error("$a acts on several sites at once, which $what cannot place: apply it as a gate, " *
              "write it as an expression of one site operators, or create it with the sites " *
              "it acts on, Operator{N}(name, m, type, sites...), which splits it into one")
    end
    return nothing
end
