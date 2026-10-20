# Systems, the sites a state lives on, compared by their content, and the ITensor indices of
# their sites, whose identities derive from that content. The later layers also use `is_charged`.

export System, siteind, siteinds

"""
    System(sites::Vector{<:AbstractSite})
    System(n::Int, site::AbstractSite)

a quantum system, made of its sites. `System(n, site)` is made of `n` copies of `site`. Two
systems of the same sites are equal, and have the same indices. Sites whose conserved
quantities cannot live together on one system are refused.

# Examples

    System(10, Qubit())
    System([Qubit(), Spin(1), Qubit(), Boson(5)])
"""
struct System
    sites::Vector{<:AbstractSite}
    # the hash of every site, from which the indices take their identities
    key::UInt
    # read once off the sites, see `is_charged`, `strong_names` and `next_fermion`
    charged::Bool
    strong::Vector{String}
    fermions::Vector{Int}
    # the systems it weakens to, see `weaken`, made once each, under `lock`, so that the states
    # of this system share theirs, and the operators prepared to be measured on its states, see
    # `prepare_op`, with whether they are Hermitian, those asked for once leaving their hash in
    # `seen`, see `kept_for!`
    weakened::Dict{Conserved, System}
    prepared::Dict{IndexedOp{IsOperator}, Any}
    hermitian::Dict{IndexedOp{IsOperator}, Bool}
    seen::Set{UInt}
    lock::ReentrantLock
    function System(sites::Vector{<:AbstractSite})
        foreach(check_site, unique(sites))
        check_charges(sites)
        charged = any(s -> !isempty(conserved(s).quantities), sites)
        strong = unique(String[ x for s in sites for x in strong_names(s) ])
        n = length(sites)
        fermions = fill(n + 1, n + 1)
        for k in n:-1:1
            fermions[k] = any(isodd, fermion_parities(sites[k])) ? k : fermions[k + 1]
        end
        # Julia's hash of a long vector reads only some of its elements
        key = foldl((h, s) -> hash(s, h), sites; init = hash(length(sites)))
        return new(copy(sites), key, charged, strong, fermions, Dict{Conserved, System}(),
                   Dict{IndexedOp{IsOperator}, Any}(), Dict{IndexedOp{IsOperator}, Bool}(),
                   Set{UInt}(), ReentrantLock())
    end
end

System(size::Int, a::AbstractSite) = System(fill(a, size))

"""
    prepared_limit

how many values a dictionary of a system keeps, and how many keys asked for once it remembers,
see `kept_for!`: one more empties it, which bounds its memory. It changes how fast, never what
is measured.
"""
const prepared_limit = 1000

"""
    kept_for!(make, system, cache, key; worth)

the value of `key` in the dictionary `cache` of `system`, made by `make()` when it is not there.
A value `worth` keeping is kept the second time its key is asked for, the first leaving only the
hash of the key in `system.seen`, so that an operator asked for once, as one changing at every
step, is never kept. It is made outside the lock of the system, so that the threads measuring
other operators wait for nothing, two threads making the same value keeping the first.
"""
function kept_for!(make, system::System, cache::AbstractDict, key; worth = _ -> true)
    found = lock(() -> get(cache, key, nothing), system.lock)
    if !isnothing(found)
        return found
    end
    v = make()
    if !worth(v)
        return v
    end
    h = hash(key)
    return lock(system.lock) do
        if !(h in system.seen)
            if length(system.seen) ≥ prepared_limit
                empty!(system.seen)
            end
            push!(system.seen, h)
            return v
        end
        if length(cache) ≥ prepared_limit
            empty!(cache)
        end
        return get!(cache, key, v)
    end
end

(a::System == b::System) = a.sites == b.sites
hash(a::System, h::UInt) = hash(a.sites, hash(System, h))

getindex(s::System, i...) = s.sites[i...]

"""
    length(::System)

the number of sites of the system.
"""
length(system::System) = length(system.sites)

show(io::IO, s::System) = print(io, "System($(s.sites))")

"""
    is_charged(system)

whether the system carries quantum numbers, that is whether any of its sites conserves
something, which makes every index of the system charged
"""
is_charged(system::System) = system.charged

"""
    strong_names(::System)

the names of the quantities the sites of the system conserve strongly
"""
strong_names(system::System) = system.strong

"""
    next_fermion(system, k)

the first fermionic site of `system` from site `k` on, `length(system) + 1` when there is none
"""
next_fermion(system::System, k::Int) = system.fermions[k]

symmetries(system::System) =
    Conserved(unique(Tuple{String, Bool}[ x for s in system.sites for x in symmetries(s).names ]))

"""
    weaken(::System)
    weaken(::System, target)

the same sites, each keeping of what it conserves only what `target` names, as strongly as
`target` asks, written as `conserve` is given or as `symmetries` returns it. Without a target,
one level down: strong quantities made weak, or, when none is strong, all dropped.

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
    return lock(system.lock) do
        get!(() -> System([ weaken(s, target) for s in system.sites ]), system.weakened, target)
    end
end

weaken(system::System, spec) = weaken(system, Conserved(spec_names(spec)))

weaken(system::System) = weaken(system, one_step_down(symmetries(system)))

"""
    ::System ⊗ ::System
    tensor(::System, ::System)

the system of the sites of the first followed by those of the second.
"""
(sys1::System ⊗ sys2::System) = System([sys1.sites ; sys2.sites])

"""
    index_id(system, i, R)

the identity of the index of site `i` of `system` in the representation `R`, a hash of the
content of the system, of the position and of the representation
"""
index_id(system::System, i::Int, ::Type{R}) where R = UInt64(hash((i, nameof(R)), system.key))

"""
    site_index(system, i)

the index of site `i` of `system` in the pure representation, in the basis of the site. When
another site conserves something, a site conserving nothing takes a single sector of charge
zero. The sectors come one per basis state: merged by charge, they would reorder the basis
whenever equal charges are not contiguous, as those of `parity(N)` on a boson, 0, 1, 0, 1.
"""
function site_index(system::System, i::Int)
    site = system[i]
    n = dim(site)
    # `nameof`: `string(typeof(site))` may carry the module, and ITensors cuts a tag at 16
    # characters
    tg = "$(nameof(typeof(site))),Site"
    x = if !is_charged(system)
        Index(n; tags = tg)
    elseif isempty(conserved(site).quantities)
        Index(QN() => n; tags = tg)
    else
        Index([ q => 1 for q in basis_charges(site) ]...; tags = tg)
    end
    return Index(index_id(system, i, Pure), space(x), dir(x), tags(x), plev(x))
end

"""
    mixed_space_index(i, site)

an index of the mixed representation combining the ket index `i` with the bra of the same
site, daggered so that ``|m\\rangle\\langle n|`` carries the difference of their charges, and
starred under a strong symmetry, which keeps the two charges apart
"""
mixed_space_index(i::Index, site::AbstractSite) =
    addtags(combinedind(combiner(i, dag(bra_index(i, site)'); tags = tags(i))), "Mixed")

"""
    siteind(Pure, system, i)
    siteind(Mixed, system, i)
    siteinds(R, system)

the ITensor index of site `i` of `system`, for pure or for mixed representations, or those of
all its sites.

# Examples

    siteind(Pure, system, 1)
    siteinds(Mixed, system)
"""
siteind(::Type{Pure}, system::System, i::Int) = site_index(system, i)

function siteind(::Type{Mixed}, system::System, i::Int)
    x = mixed_space_index(site_index(system, i), system[i])
    return Index(index_id(system, i, Mixed), space(x), dir(x), tags(x), plev(x))
end

siteinds(::Type{R}, system::System) where R = [ siteind(R, system, i) for i in 1:length(system) ]
