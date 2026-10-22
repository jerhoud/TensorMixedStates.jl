# Prepared states with an exact MPS: the GHZ, Dicke and W states, the fully mixed state of a
# sector, and the states given by their MPS tensors, their dense vector or their density
# matrix; and the dense vector of a pure state.

export mps_state, dense_state, ghz_state, superposition, mixture, dicke_state, w_state,
       fully_mixed

"""
    local_vector(site, st)

the vector of the local pure state `st` of `site`, refused when it is mixed
"""
function local_vector(site::AbstractSite, st)
    v = local_state(site, st)
    if !(v isa Vector)
        error("$(repr(st)) is a density matrix on $site, where a pure state is needed")
    end
    return v
end

"""
    mps_state(system, tensors)

the pure state of `system`, normalized, whose MPS has the tensors `tensors`, each an array
`A[l, s, r]` in the basis of its site, the first of left dimension 1 and the last of right
dimension 1. On a system carrying charges, a state of no definite charge is refused.

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
    s = siteinds(Pure, system)
    tensors = collect(tensors)
    links = if is_charged(system)
        # none on a site conserving nothing, whose index has a single sector of charge zero
        charge(k, j) = isempty(conserved(system[k]).quantities) ? QN() : basis_charges(system[k])[j]
        partial = Union{Nothing, QN}[ QN() ]
        map(1:n) do k
            # the elements leaving a link state no element reaches multiply a zero amplitude:
            # they are dropped, so that the charge any of them would carry is not checked
            unreached = findall(isnothing, partial)
            if !isempty(unreached)
                tensors[k] = copy(tensors[k])
                tensors[k][unreached, :, :] .= 0
            end
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
            partial = next
            # the charges of the left part taken out: each tensor but the last has no flux; a
            # link state no element reaches carries nothing, any charge will do
            Index([ -something(q, QN()) => 1 for q in partial ]...; tags = "Link,l=$k")
        end
    else
        [ Index(size(tensors[k], 3); tags = "Link,l=$k") for k in 1:n ]
    end
    edge = is_charged(system) ? Index(QN() => 1; tags = "Link,l=0") :
                                Index(1; tags = "Link,l=0")
    its = map(1:n) do k
        l = k == 1 ? edge : dag(links[k-1])
        t = ITensor(tensors[k], l, s[k], links[k])
        # an MPS has no edge links
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
    dense_mps(a, idx, sites, limits; combiners = ITensor[])

the MPS on `sites` of the array `a` on the indices `idx`, once the combiners `combiners` have
gathered these indices into the sites, truncated to `limits`
"""
function dense_mps(a::AbstractArray, idx, sites, limits::Limits; combiners = ITensor[])
    if !has_definite_flux(a, idx)
        error("the state has no definite charge, which a system conserving it cannot hold")
    end
    t = charged_itensor(a, idx)
    for c in combiners
        t *= c
    end
    return MPS(t, sites; truncation(limits)...)
end

"""
    dense_state(system, ψ; limits = Limits())
    dense_state(system, ρ; limits = Limits())

the pure state of `system` of vector `ψ`, or the mixed state of density matrix `ρ`,
normalized, on the basis of product states ordered as `kron` orders them, the first site
varying the slowest. Its MPS is decomposed from the whole array and truncated to `limits`, by
default `cutoff = 1e-26`. On a system carrying charges, a state of no definite charge is refused.

# Examples

    dense_state(System(2, Qubit()), [1, 0, 0, 1])               # (|↑↑⟩ + |↓↓⟩)/√2
    dense_state(System(2, Qubit()), [1 0 0 0 ; 0 0 0 0 ; 0 0 0 0 ; 0 0 0 1])
"""
function dense_state(system::System, ψ::AbstractVector{<:Number}; limits::Limits = Limits())
    s = siteinds(Pure, system)
    if length(ψ) ≠ prod(dim, s)
        error("a vector of $(length(ψ)) elements cannot be a state of $(prod(dim, s)) basis states")
    end
    a, idx = on_legs(reshape(collect(ψ), :, 1), s, Index[])
    return normalize(State{Pure}(system, dense_mps(a, idx, s, limits)))
end

function dense_state(system::System, ρ::AbstractMatrix{<:Number}; limits::Limits = Limits())
    n = length(system)
    s = siteinds(Pure, system)
    d = prod(dim, s)
    if size(ρ) ≠ (d, d)
        error("a matrix of size $(size(ρ)) cannot be a density matrix of $d basis states")
    end
    mixers = [ mixer(s[k], siteind(Mixed, system, k), system[k]) for k in 1:n ]
    a, idx = on_legs(Matrix(ρ), s, [ dag(b') for (b, _) in mixers ])
    sites = siteinds(Mixed, system)
    return normalize(State{Mixed}(system, dense_mps(a, idx, sites, limits;
                                                    combiners = last.(mixers))))
end

"""
    Vector(state::State{Pure})

the vector of the pure `state`, not normalized, on the basis of the product states ordered as
`kron` orders them, the first site varying the slowest, as `dense_state` takes it and
`matrix(op, system)` acts on it. It is refused beyond `dense_max^2` elements, as many as the
largest dense matrix TMS builds. The density matrix of a mixed state is given by
`reduced_density_matrix(ρ, 1:n)`.

# Examples

    Vector(State{Pure}(System(2, Qubit()), "+"))      # [0.5, 0.5, 0.5, 0.5]
"""
function Vector(state::State{Pure})
    system = state.system
    # a Float64, which does not wrap around as an Int does on many sites
    d = prod(float(dim(s)) for s in system.sites)
    if d > dense_max^2
        error("the vector of a state of $(length(system)) sites has " *
              "$(Printf.@sprintf("%.0f", d)) elements, beyond the $(dense_max^2) TMS builds")
    end
    # site by site, the vector of the sites so far for each value of the link on their right: the
    # tensor of all the sites, of n indices, would be compiled anew for each n, in 10 to 20 s
    st = state.state
    n = length(st)
    v = ones(1, 1)
    for k in 1:n
        left = k > 1 ? commonind(st[k-1], st[k]) : nothing
        right = k < n ? commonind(st[k], st[k+1]) : nothing
        idx = [ i for i in (left, siteind(Pure, system, k), right) if !isnothing(i) ]
        a = reshape(Array(dense(st[k]), idx...), size(v, 2), :)
        # the new site varying the fastest
        b = reshape(v * a, size(v, 1), dim(system[k]), :)
        v = reshape(permutedims(b, (2, 1, 3)), :, size(b, 3))
    end
    return vec(v)
end

"""
    product_sum(R, system, terms)

the MPS, in the representation `type`, of the sum of the product states of `terms`, pairs
`c => states`, with a channel per term of nonzero coefficient; on a charged system the terms
must have the same charge
"""
function product_sum(::Type{R}, system::System, terms) where R
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
        [ make_one_state(R, system, k, sts[k]) for k in 1:n ]
    end
    kept = findall(!iszero ∘ first, terms)
    if isempty(kept)
        error("a sum of product states needs a nonzero term at least")
    end
    coefs = first.(terms[kept])
    locals = locals[kept]
    charged = hasqns(locals[1][1])
    # the charge of the sites up to each one, for each term
    partial = charged ? [ accumulate(+, flux.(ts)) for ts in locals ] : nothing
    if charged && !allequal(last.(partial))
        error("the terms have different charges: the sum has no definite charge, which a " *
              "system conserving it cannot hold")
    end
    # each link carries, daggered, the charge of the sites on its left
    links = map(1:n-1) do k
        charged ? dag(Index([ p[k] => 1 for p in partial ]...; tags = "Link,l=$k")) :
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
the number of terms of nonzero coefficient, truncated to `limits` when they are given, by
default `cutoff = 0`. On a system conserving something, the product states of these terms must
have the same charge.

# Examples

    superposition(System(4, Qubit()),
                  [1 => ["Up", "Dn", "Up", "Dn"], -1 => ["Dn", "Up", "Dn", "Up"]])
    superposition(System(6, Qubit()), [1 => "Up", 1 => "Dn"])     # the GHZ state
"""
function superposition(system::System, terms::AbstractVector{<:Pair}; limits::Limits = Limits())
    st = State{Pure}(system, product_sum(Pure, system, terms))
    if limits ≠ Limits()
        st = State(st, truncate(st.state; truncation(limits; cutoff = 0.)...))
    end
    return normalize(st)
end

"""
    mixture(system, terms; limits = Limits())

the mixed state ``\\sum_k p_k \\rho^k_1 \\otimes \\dots \\otimes \\rho^k_n`` of `terms`, pairs
`p => states` of a weight, real and not negative, and the local states of a product state, pure
or mixed, given by their names, their vectors or their density matrices, one for every site or
one for all of them, normalized to trace one. It is exact, of bond dimension the number of
terms of nonzero weight, truncated to `limits` when they are given, by default `cutoff = 0`.

# Examples

    mixture(System(4, Qubit()), [0.5 => ["Up", "Dn", "Up", "Dn"], 0.5 => ["Dn", "Up", "Dn", "Up"]])
"""
function mixture(system::System, terms::AbstractVector{<:Pair}; limits::Limits = Limits())
    for (p, _) in terms
        if !(p isa Real) || p < 0
            error("the weights of a mixture are real and not negative, and $p is not")
        end
    end
    st = State{Mixed}(system, product_sum(Mixed, system, terms))
    if limits ≠ Limits()
        st = State(st, truncate(st.state; truncation(limits; cutoff = 0.)...))
    end
    return normalize(st)
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
    # one element type for every tensor, which mps_state takes
    vs = [ (local_vector(system[j], a), amplitudes[j] * local_vector(system[j], b)) for j in 1:n ]
    T = mapreduce(((va, vb),) -> promote_type(eltype(va), eltype(vb)), promote_type, vs;
                  init = Bool)
    tensors = map(1:n) do j
        va, vb = vs[j]
        left, right = reach(j - 1), reach(j)
        A = zeros(T, length(left), dim(system[j]), length(right))
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
    w_state(System(40, Qubit()), "Up", "Dn",
            [ exp(-(j - 10)^2 / 8 + im * π / 2 * j) for j in 1:40 ])
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
    fully_mixed(system, quantities...)

the fully mixed state of the sector given by `quantities`, pairs `op => value` where the sum
of `op` over the sites takes the value: the projector on the basis states of the sector divided
by their number, from which `Thermalize` gives the canonical thermal state of a Hamiltonian
conserving these quantities. Each `op` must be diagonal on the basis of every site, as a number
of particles or an ``S^z``, but the sites need not conserve it. Without quantities, it is the
fully mixed state of the system, which a system conserving something strongly refuses: give a
sector fixing what it conserves strongly. It is exact, of bond dimension `N + 1` at most for
`N => N`.

# Examples

    fully_mixed(System(10, Fermion()), N => 4)
    fully_mixed(System(8, Electron(conserve = (strong(Ntot), 2Sz))), Ntot => 8, Sz => 0)
"""
function fully_mixed(system::System, quantities::Pair...)
    n = length(system)
    if isempty(quantities)
        return State{Mixed}(system, "FullyMixed")
    end
    # rounded so that sums of the same values compare equal
    key(x) = round(x; digits = 10)
    values = map(1:n) do k
        ms = [ matrix(op, system[k]) for (op, _) in quantities ]
        for (m, (op, _)) in zip(ms, quantities)
            if !isdiag(m)
                error("$(obs_name(op)) is not diagonal on the basis of $(system[k])")
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
              join(("$(obs_name(op)) = $v" for (op, v) in quantities), ", "))
    end
    projector(d, j) = [ a == j && b == j ? 1. : 0. for a in 1:d, b in 1:d ]
    locals = [ [ make_one_state(Mixed, system, k, projector(dim(system[k]), j))
                 for j in 1:dim(system[k]) ] for k in 1:n ]
    charged = hasqns(locals[1][1])
    charge(t) = charged ? flux(t) : nothing
    # a link state is the values of the quantities on its left with their charge, which differs
    # between paths when the system conserves something else
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
    # each link carries, daggered, the charge of the sites on its left
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
