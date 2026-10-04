# Prepared pure states with an exact MPS of small bond dimension, built from their tensors: the
# GHZ state, the Dicke states, of which the W state, and the states of singlets on pairs.

export ghz_state, dicke_state, w_state, dimer_state

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

"""
    dicke_state(system, k, a, b)

the Dicke state of `k` sites in the local pure state `b` and the others in `a`, given by their
names or their vectors: the superposition with equal weights of the product states that place
`b` on `k` sites and `a` on the others, ``\\binom{n}{k}^{-1/2} \\sum |a\\dots b \\dots
a\\rangle`` for orthonormal states, normalized in any case. It is exact, of bond dimension
`k + 1` at most, and holds the charge of its product states, which share it.

# Examples

    dicke_state(System(10, Qubit()), 3, "Up", "Dn")
    dicke_state(System(10, Qubit(conserve = N)), 5, "Up", "Dn")
"""
function dicke_state(system::System, k::Int, a, b)
    n = length(system)
    if !(0 ≤ k ≤ n)
        error("$k sites in $b do not fit on $n sites")
    end
    # the numbers of sites in b on the left of a link that can still reach k
    reach(j) = max(0, k - (n - j)):min(j, k)
    tensors = map(1:n) do j
        va, vb = local_vector(system[j], a), local_vector(system[j], b)
        left, right = reach(j - 1), reach(j)
        A = zeros(promote_type(eltype(va), eltype(vb)), length(left), dim(system[j]),
                  length(right))
        for (p, c) in enumerate(left), (q, d) in enumerate(right)
            if d == c
                A[p, :, q] = va
            elseif d == c + 1
                A[p, :, q] = vb
            end
        end
        A
    end
    return mps_state(system, tensors)
end

"""
    w_state(system, a, b)

the W state, the Dicke state of one site in `b` and the others in `a`, see `dicke_state`.

# Examples

    w_state(System(10, Qubit()), "Up", "Dn")
"""
w_state(system::System, a, b) = dicke_state(system, 1, a, b)

"""
    dimer_state(system, pairs, a, b; others, limits = Limits())

the state of a singlet ``(|ab\\rangle - |ba\\rangle)/\\sqrt{2}`` on each pair `(i, j)` of
`pairs`, `a` on `i` in the first term, of the local pure states `a` and `b`, given by their
names or their vectors and orthonormal: `"Up"` and `"Dn"` for spins one half. The sites in no
pair take the local state `others`, one for all of them or a vector of one for each, which
they then need.

Each singlet comes from a gate of two sites applied to ``|ab\\rangle``, on sites apart as well,
so that pairs of neighbours give the Majumdar-Ghosh state, of bond dimension 2, and nested
pairs `(i, n + 1 - i)` the rainbow state, of bond dimension ``2^{n/2}`` in the middle: `limits`
constrains the truncations made as the gates are applied.

# Examples

    dimer_state(System(10, Qubit()), [ (i, i + 1) for i in 1:2:9 ], "Up", "Dn")
    dimer_state(System(9, Spin(1/2)), [ (i, i + 1) for i in 1:2:7 ], "1/2", "-1/2";
                others = "1/2")
"""
function dimer_state(system::System, pairs::AbstractVector{Tuple{Int, Int}}, a, b;
                     others = nothing, limits::Limits = Limits())
    n = length(system)
    paired = reduce(vcat, [ [i, j] for (i, j) in pairs ]; init = Int[])
    if !allunique(paired)
        error("a site is in two pairs of $pairs")
    end
    if any(k -> !(1 ≤ k ≤ n), paired)
        error("the pairs $pairs reach beyond the $n sites of the system")
    end
    rest = setdiff(1:n, paired)
    states = Vector{Any}(undef, n)
    for (i, j) in pairs
        states[i], states[j] = a, b
    end
    if !isempty(rest)
        if isnothing(others)
            which = length(rest) == 1 ? "site $(only(rest)) is in no pair: give its" :
                                        "sites $(join(rest, ", ")) are in no pair: give their"
            error("$which state with others")
        end
        local_states = others isa AbstractVector && !(eltype(others) <: Number) ? others :
                       fill(others, length(rest))
        if length(local_states) ≠ length(rest)
            error("others gives $(length(local_states)) states for the $(length(rest)) sites " *
                  "in no pair")
        end
        states[rest] = local_states
    end
    st = State{Pure}(system, states)
    gates = map(pairs) do (i, j)
        va, vb = local_vector(system[i], a), local_vector(system[i], b)
        wa, wb = local_vector(system[j], a), local_vector(system[j], b)
        if abs(va' * vb) > 1e-12 || abs(wa' * wb) > 1e-12
            error("a singlet of $a and $b needs them orthonormal")
        end
        # on the basis of the two sites, the gate takes |ab⟩ to the singlet and |ba⟩ to the
        # triplet of no charge, and leaves alone what is orthogonal to both
        x, y = kron(va, wb), kron(vb, wa)
        U = I + ((x - y) / sqrt(2) - x) * x' + ((x + y) / sqrt(2) - y) * y'
        named(U, "Singlet", system[i], system[j])(i, j)
    end
    return isempty(gates) ? st : apply(prod(gates), st; limits)
end
