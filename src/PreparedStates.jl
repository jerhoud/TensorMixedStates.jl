# Prepared states with an exact MPS of small bond dimension, built from their tensors: the GHZ
# state, the Dicke states, of which the W state, the states of singlets on pairs, and the fully
# mixed state of a sector, and the states given by the tensors of their MPS or by their dense
# vector or density matrix.

export mps_state, dense_state, ghz_state, superposition, mixture, dicke_state, w_state,
       dimer_state, fully_mixed

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

# Examples

    # the cluster state of a chain of qubits, each link carrying the state of the site on its
    # left, and a sign when both sites are in "Dn"
    A = zeros(2, 2, 2)
    for l in 1:2, s in 1:2
        A[l, s, s] = l == s == 2 ? -1 : 1
    end
    mps_state(System(10, Qubit()), [ A[1:1, :, :], fill(A, 8)..., sum(A; dims = 3) ])
"""
function mps_state(system::System, tensors::Vector{<:AbstractArray{<:Number, 3}})
    n = length(system)
    if length(tensors) ≠ n
        error("$(length(tensors)) tensors cannot make the MPS of $n sites")
    end
    for k in 1:n
        A = tensors[k]
        if size(A, 2) ≠ dim(system[k])
            error("tensor $k has $(size(A, 2)) states for a site of dimension $(dim(system[k]))")
        end
        left = k == 1 ? 1 : size(tensors[k-1], 3)
        if size(A, 1) ≠ left
            error("tensor $k has a left dimension of $(size(A, 1)) where $left is expected")
        end
    end
    if size(tensors[n], 3) ≠ 1
        error("the last tensor has a right dimension of $(size(tensors[n], 3)) where 1 is expected")
    end
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
    edge = is_charged(system) ? Index(QN() => 1; tags = "Link,l=0") :
                                Index(1; tags = "Link,l=0")
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
    dense_mps(a, idx, sites; combiners = ITensor[])

the MPS on `sites` of the array `a` on the indices `idx`, once the combiners `combiners` have
gathered these indices into the sites, cut at the rounding as `sum_cutoff` cuts a sum
"""
function dense_mps(a::AbstractArray, idx, sites; combiners = ITensor[])
    if !has_definite_flux(a, idx)
        error("the state has no definite charge, which a system conserving it cannot hold")
    end
    t = charged_itensor(a, idx)
    for c in combiners
        t *= c
    end
    return MPS(t, sites; cutoff = sum_cutoff(Limits()))
end

"""
    dense_state(system, ψ; limits = Limits())
    dense_state(system, ρ; limits = Limits())

the pure state of `system` of vector `ψ`, or the mixed state of density matrix `ρ`, on the
basis of the product states of the sites ordered as `kron` orders them, the first site varying
the slowest, normalized. Its MPS comes from successive decompositions of the whole vector, exact
up to rounding and truncated to `limits` when they are given, for a system small enough to be
written down: to compare with an exact diagonalization, or to take a state from another code.
On a system carrying charges, a state of no definite charge is refused.

# Examples

    dense_state(System(2, Qubit()), [1, 0, 0, 1])               # (|↑↑⟩ + |↓↓⟩)/√2
    dense_state(System(2, Qubit()), [1 0 0 0 ; 0 0 0 0 ; 0 0 0 0 ; 0 0 0 1])
"""
function dense_state(system::System, ψ::AbstractVector{<:Number}; limits::Limits = Limits())
    s = [ SysIndex{Pure}(system, k) for k in 1:length(system) ]
    if length(ψ) ≠ prod(dim, s)
        error("a vector of $(length(ψ)) elements cannot be a state of $(prod(dim, s)) basis states")
    end
    a, idx = on_legs(reshape(collect(ψ), :, 1), s, Index[])
    st = normalize(State{Pure}(system, dense_mps(a, idx, s)))
    return limits == Limits() ? st : truncate(st; limits)
end

function dense_state(system::System, ρ::AbstractMatrix{<:Number}; limits::Limits = Limits())
    n = length(system)
    s = [ SysIndex{Pure}(system, k) for k in 1:n ]
    d = prod(dim, s)
    if size(ρ) ≠ (d, d)
        error("a matrix of size $(size(ρ)) cannot be a density matrix of $d basis states")
    end
    mixers = [ mixer(s[k], SysIndex{Mixed}(system, k), system[k]) for k in 1:n ]
    a, idx = on_legs(Matrix(ρ), s, [ dag(b') for (b, _) in mixers ])
    sites = [ SysIndex{Mixed}(system, k) for k in 1:n ]
    st = normalize(State{Mixed}(system, dense_mps(a, idx, sites; combiners = last.(mixers))))
    return limits == Limits() ? st : truncate(st; limits)
end

"""
    product_sum(type, system, terms)

the sum of the product states of `terms`, pairs `c => states` of a coefficient and the local
states of a product state, one for every site or one for all of them, in the representation
`type`, `Pure()` or `Mixed()`. Its MPS has a channel per term of nonzero coefficient, each
carrying the charge of the sites of its product state on the left of the link, so that these
terms must have the same charge on a system conserving something.
"""
function product_sum(type::PM, system::System, terms)
    n = length(system)
    if isempty(terms)
        error("a sum of product states needs one term at least")
    end
    per_site(st) = st isa AbstractVector && !(eltype(st) <: Number) ? st : fill(st, n)
    locals = map(terms) do (_, states)
        sts = per_site(states)
        if length(sts) ≠ n
            error("a product state of $(length(sts)) local states cannot be one of $n sites")
        end
        [ make_one_state(type, system, k, sts[k]) for k in 1:n ]
    end
    kept = findall(!iszero ∘ first, terms)
    if isempty(kept)
        error("a sum of product states needs a nonzero term at least")
    end
    coefs = first.(terms[kept])
    locals = locals[kept]
    charged = hasqns(locals[1][1])
    if charged && !allequal(sum(flux, ts) for ts in locals)
        error("the terms have different charges: the sum has no definite charge, which a " *
              "system conserving it cannot hold")
    end
    # as the links of a product state, the charge of the sites on the left of each, daggered
    links = map(1:n-1) do k
        charged ? dag(Index([ sum(flux, ts[1:k]) => 1 for ts in locals ]...; tags = "Link,l=$k")) :
                  Index(length(kept); tags = "Link,l=$k")
    end
    its = map(1:n) do k
        sum(enumerate(locals)) do (t, ts)
            x = k == 1 ? coefs[t] * ts[k] : ts[k]
            if k > 1
                x *= onehot(dag(links[k-1]) => t)
            end
            if k < n
                x *= onehot(links[k] => t)
            end
            x
        end
    end
    return MPS(its)
end

"""
    superposition(system, terms; limits = Limits())

the pure state ``\\sum_k c_k |s^k_1 s^k_2 \\dots s^k_n\\rangle`` of `terms`, pairs `c => states`
of a coefficient and the local pure states of a product state, given by their names or their
vectors, one for every site or one for all of them, normalized. It is exact, of bond dimension
the number of terms of nonzero coefficient, truncated to `limits` when they are given, and built
at once rather than term by term. On a system conserving something, the product states of these
terms must have the same charge.

# Examples

    superposition(System(4, Qubit()),
                  [1 => ["Up", "Dn", "Up", "Dn"], -1 => ["Dn", "Up", "Dn", "Up"]])
    superposition(System(6, Qubit()), [1 => "Up", 1 => "Dn"])     # the GHZ state
"""
function superposition(system::System, terms::AbstractVector{<:Pair}; limits::Limits = Limits())
    st = normalize(State{Pure}(system, product_sum(Pure(), system, terms)))
    return limits == Limits() ? st : truncate(st; limits)
end

"""
    mixture(system, terms; limits = Limits())

the mixed state ``\\sum_k p_k \\rho^k_1 \\otimes \\dots \\otimes \\rho^k_n`` of `terms`, pairs
`p => states` of a weight, real and not negative, and the local states of a product state, pure
or mixed, given by their names, their vectors or their density matrices, one for every site or
one for all of them, normalized to trace one. It is exact, of bond dimension the number of
terms of nonzero weight, truncated to `limits` when they are given: a classical ensemble of product states.

# Examples

    mixture(System(4, Qubit()), [0.5 => ["Up", "Dn", "Up", "Dn"], 0.5 => ["Dn", "Up", "Dn", "Up"]])
"""
function mixture(system::System, terms::AbstractVector{<:Pair}; limits::Limits = Limits())
    for (p, _) in terms
        if !(p isa Real) || p < 0
            error("the weights of a mixture are real and not negative, and $p is not")
        end
    end
    st = normalize(State{Mixed}(system, product_sum(Mixed(), system, terms)))
    return limits == Limits() ? st : truncate(st; limits)
end

"""
    ghz_state(system, states...)

the GHZ state of `states`, local pure states given by their names or their vectors: the
superposition with equal weights of the product states in which every site is in the same one,
``(|aa\\dots a\\rangle + |bb\\dots b\\rangle + \\dots)/\\sqrt{m}`` for orthonormal states,
normalized in any case. It is exact, of bond dimension the number of states. On sites conserving
something, the product states must have the same charge.

# Examples

    ghz_state(System(10, Qubit()), "Up", "Dn")
    ghz_state(System(6, Qudit(3)), "0", "1", "2")
"""
ghz_state(system::System, states...) =
    if isempty(states)
        error("a GHZ state needs one local state at least")
    else
        superposition(system, [ 1 => st for st in states ])
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
dicke_state(system::System, k::Int, a, b) = dicke_mps(system, k, a, b, ones(length(system)))

"""
    dicke_mps(system, k, a, b, amplitudes)

the state of `dicke_state`, each site in `b` weighted by its amplitude in `amplitudes`
"""
function dicke_mps(system::System, k::Int, a, b, amplitudes::AbstractVector{<:Number})
    n = length(system)
    if !(0 ≤ k ≤ n)
        error("$k sites in $b do not fit on $n sites")
    end
    # the numbers of sites in b on the left of a link that can still reach k
    reach(j) = max(0, k - (n - j)):min(j, k)
    tensors = map(1:n) do j
        va, vb = local_vector(system[j], a), amplitudes[j] * local_vector(system[j], b)
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
    w_state(system, a, b, amplitudes)

the W state, the Dicke state of one site in `b` and the others in `a`, see `dicke_state`. With
`amplitudes`, one for each site, the superposition ``\\sum_j \\phi_j |a \\dots a b_j a \\dots
a\\rangle`` of the states in which site `j` alone is in `b`, normalized: a wave packet of a
single excitation, of bond dimension 2.

# Examples

    w_state(System(10, Qubit()), "Up", "Dn")
    w_state(System(40, Qubit()), "Up", "Dn", [ exp(-(j - 10)^2 / 8 + im * π / 2 * j) for j in 1:40 ])
"""
w_state(system::System, a, b) = dicke_state(system, 1, a, b)

function w_state(system::System, a, b, amplitudes::AbstractVector{<:Number})
    if length(amplitudes) ≠ length(system)
        error("$(length(amplitudes)) amplitudes for $(length(system)) sites")
    end
    if all(iszero, amplitudes)
        error("a wave packet needs a nonzero amplitude at least")
    end
    return dicke_mps(system, 1, a, b, amplitudes)
end

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

"""
    fully_mixed(system, quantities...)

the fully mixed state of the sector where each quantity takes its value, given as pairs
`op => value`, the sum of `op` over the sites taking the value: the projector on the basis
states of the sector divided by their number, the state at infinite temperature of a sector,
from which `Thermalize` gives the canonical thermal state of a hamiltonian conserving these
quantities. An operator is diagonal on the
basis of every site, as a number of particles or a ``S^z``, and the sites need not conserve
it. Without quantities, it is the fully mixed state of the system, which a system conserving
something strongly refuses, its sectors being apart: a sector that fixes what it conserves
strongly is accepted.

It is exact, its bond dimension the number of values the quantities take on the sites on the
left of a link that can still reach the sector, `N + 1` at most for `N => N`.

# Examples

    fully_mixed(System(10, Fermion()), N => 4)
    fully_mixed(System(8, Electron(conserve = (strong(Ntot), 2Sz))), Ntot => 8, Sz => 0)
"""
function fully_mixed(system::System, quantities::Pair...)
    n = length(system)
    if isempty(quantities)
        return State{Mixed}(system, "FullyMixed")
    end
    # the values of the quantities on each basis state of each site, rounded so that sums of
    # the same values are equal, half integers and integers being exact anyway
    key(x) = round(x; digits = 10)
    values = map(1:n) do k
        ms = [ matrix(op, system[k]) for (op, _) in quantities ]
        for (m, (op, _)) in zip(ms, quantities)
            if !isdiag(m)
                error("$op is not diagonal on the basis of $(system[k])")
            end
        end
        [ Tuple(key(real(m[j, j])) for m in ms) for j in 1:dim(system[k]) ]
    end
    target = Tuple(key(Float64(v)) for (_, v) in quantities)
    add(a, b) = key.(a .+ b)
    # the values the sites on the right of each link can still bring to the target
    needed = Vector{Set{Tuple}}(undef, n + 1)
    needed[n+1] = Set{Tuple}([target])
    for k in n:-1:1
        needed[k] = Set{Tuple}(key.(v .- q) for v in needed[k+1] for q in values[k])
    end
    if !(Tuple(0. for _ in quantities) in needed[1])
        error("no basis state of the system has " *
              join(("$op = $v" for (op, v) in quantities), ", "))
    end
    # the local tensor of each basis state |s><s|, on the mixed index, with its charge
    projector(d, j) = [ a == j && b == j ? 1. : 0. for a in 1:d, b in 1:d ]
    locals = [ [ make_one_state(Mixed(), system, k, projector(dim(system[k]), j))
                 for j in 1:dim(system[k]) ] for k in 1:n ]
    charged = hasqns(locals[1][1])
    charge(t) = charged ? flux(t) : nothing
    # the states of each link: the values of the quantities on the sites on its left, with
    # their charge, which differs between paths when the system conserves something else
    states = Vector{Vector{Tuple}}(undef, n + 1)
    states[1] = [ (Tuple(0. for _ in quantities), charged ? QN() : nothing) ]
    terms = Vector{Vector{Tuple{Int, Int, Int}}}(undef, n)
    for k in 1:n
        next = Tuple[]
        ts = Tuple{Int, Int, Int}[]
        for (l, (v, q)) in enumerate(states[k]), j in 1:dim(system[k])
            w = add(v, values[k][j])
            if w in needed[k+1]
                r = (w, charged ? q + charge(locals[k][j]) : nothing)
                i = findfirst(==(r), next)
                if isnothing(i)
                    push!(next, r)
                    i = length(next)
                end
                push!(ts, (l, j, i))
            end
        end
        states[k+1], terms[k] = next, ts
    end
    if length(states[n+1]) > 1
        error("the sector does not fix the charges the system conserves strongly, which a " *
              "mixed state cannot spread over")
    end
    # links as those of a product state, carrying the charge of the sites on their left
    links = [ charged ? dag(Index([ q => 1 for (_, q) in states[k+1] ]...; tags = "Link,l=$k")) :
                        Index(length(states[k+1]); tags = "Link,l=$k") for k in 1:n-1 ]
    its = map(1:n) do k
        sum(terms[k]) do (l, j, r)
            t = locals[k][j]
            if k > 1
                t *= onehot(dag(links[k-1]) => l)
            end
            if k < n
                t *= onehot(links[k] => r)
            end
            t
        end
    end
    return normalize(State{Mixed}(system, MPS(its)))
end
