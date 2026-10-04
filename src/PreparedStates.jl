# Prepared pure states with an exact MPS of small bond dimension, built from their tensors: the
# GHZ state, the Dicke states, of which the W state, and the states of singlets on pairs.

export ghz_state

"""
    local_vector(site, st)

the vector of the local pure state `st` of `site`, a name or a vector
"""
local_vector(site::AbstractSite, st::String) = state(site, st)
local_vector(::AbstractSite, st::Vector{<:Number}) = st

"""
    mps_state(system, tensors)

the pure state of `system` whose MPS has the tensors `tensors`, each an array `A[l, s, r]` in the
basis of its site, the first of left dimension 1 and the last of right dimension 1, normalized.
On a system carrying charges, the charge of each link is read off the tensors, each of their
elements adding the charge of its basis state to that of its left link, and a state of no
definite charge is refused: a superposition of different charges has no MPS of charged indices.
"""
function mps_state(system::System, tensors::Vector{<:AbstractArray{<:Number, 3}})
    n = length(system)
    s = [ SysIndex{Pure}(system, k) for k in 1:n ]
    links = if is_charged(system)
        # the charge of each basis state of each site, one block of its index each, and none on
        # a site conserving nothing beside sites that do
        charge(k, j) = isempty(conserved(system[k])) ? QN() : basis_charges(system[k])[j]
        partial = QN[ QN() ]
        map(1:n) do k
            A = tensors[k]
            next = Vector{Union{Nothing, QN}}(nothing, size(A, 3))
            for l in axes(A, 1), j in axes(A, 2), r in axes(A, 3)
                if !iszero(A[l, j, r])
                    q = partial[l] + charge(k, j)
                    if isnothing(next[r])
                        next[r] = q
                    elseif next[r] ≠ q
                        error("the state has no definite charge, which a system conserving " *
                              "it cannot hold")
                    end
                end
            end
            # a link state no element reaches carries nothing: any charge will do
            partial = QN[ something(q, QN()) for q in next ]
            # the charges of a link are those of its left part taken out, so that each tensor
            # but the last has no flux and the last the charge of the state
            Index([ -q => 1 for q in partial ]...; tags = "Link,l=$k")
        end
    else
        [ Index(size(tensors[k], 3); tags = "Link,l=$k") for k in 1:n ]
    end
    edge = is_charged(system) ? Index(QN() => 1; tags = "Link,l=0") : Index(1; tags = "Link,l=0")
    its = map(1:n) do k
        l = k == 1 ? edge : dag(links[k-1])
        t = ITensor(tensors[k], l, s[k], links[k])
        # the dimensions of one on the edges, which an MPS does not have
        if k == 1
            t *= onehot(dag(edge) => 1)
        end
        if k == n
            t *= onehot(dag(links[n]) => 1)
        end
        t
    end
    return normalize(State{Pure}(system, MPS(its)))
end

"""
    ghz_state(system, states...)

the GHZ state of `states`, local pure states given by their names or their vectors: the
superposition with equal weights of the product states in which every site is in the same one,
``(|aa\\dots a\\rangle + |bb\\dots b\\rangle + \\dots)/\\sqrt{m}`` for orthonormal states, normalized
in any case. It is exact, of bond dimension the number of states. On sites conserving
something, the product states must have the same charge.

# Examples

    ghz_state(System(10, Qubit()), "Up", "Dn")
    ghz_state(System(6, Qudit(3)), "0", "1", "2")
"""
function ghz_state(system::System, states...)
    if isempty(states)
        error("a GHZ state needs one local state at least")
    end
    n = length(system)
    m = length(states)
    tensors = map(1:n) do k
        vs = [ local_vector(system[k], st) for st in states ]
        A = zeros(promote_type(map(eltype, vs)...), k == 1 ? 1 : m, dim(system[k]),
                  k == n ? 1 : m)
        # added rather than set: on a single site, every state falls on the same element
        for b in 1:m
            A[k == 1 ? 1 : b, :, k == n ? 1 : b] += vs[b]
        end
        A
    end
    return mps_state(system, tensors)
end
