# The measurements of a state that are not expectation values: its purity and its hermiticity,
# the entropies and the reduced density matrix of some of its sites, the negativity between two
# parts, the fidelity of two states and the variance of an operator.

export purity, hermiticity, hermitianize, renyi2, mutual_info_renyi2, entanglement_entropy,
       entanglement_by_sector, reduced_density_matrix, vonneumann_entropy, log_negativity,
       fidelity, hs_fidelity, variance

"""
    purity(::State)

the purity ``\\mathrm{tr}(\\rho^\\dagger\\rho) / (\\mathrm{tr}\\,\\rho)^2``, which is
``\\mathrm{tr}(\\rho^2)`` for a density matrix of trace one: 1 for a pure state, less for a
mixed one, above 1 only for a matrix that is not positive or not Hermitian.
"""
purity(::State{Pure}) = 1.
purity(state::State{Mixed}) = (norm(state.state) / nonzero_trace(state, "purity"))^2

"""
    adjoint_tensor(state, i)

the tensor on site `i` of the adjoint of the density matrix of `state`, which conserves nothing
strongly: the ket and the bra of the site exchanged by renaming, as a combiner in the other
order would ask for a mixed index of the opposite charge
"""
function adjoint_tensor(state::State{Mixed}, i::Int)
    system = state.system
    j = siteind(Pure, system, i)
    c = last(mixer(j, siteind(Mixed, system, i), system[i]))
    return replaceinds(dag(state.state[i]) * c, (dag(j), j'), (dag(j'), j)) * c
end

"""
    dag(::State)

the state whose density matrix is the adjoint of that of the state, computed once and kept with
it. A pure representation is refused.
"""
dag(::State{Pure}) = error("dag is meaningless on pure representations")

function dag(state::State{Mixed})
    a = cached!(state.preobs.adjoint, state) do a, state
        system = state.system
        st = strong_names(system)
        mps = if isempty(st)
            MPS([ adjoint_tensor(state, i) for i in 1:length(state) ])
        else
            # `conj`: the directions stay, only the charges are relabelled, and each element
            # |x⟩⟨y| sent to |y⟩⟨x|
            relab = relabeller(i -> adjoint_index(i, st))
            maps = element_maps(system, system, t -> relabel(t, relab); swap = true)
            MPS([ relabel(conj(state.state[i]), relab) * m for (i, m) in enumerate(maps) ])
        end
        push!(a, State(state, mps))
    end
    return only(a)
end

"""
    hermitianize(state; limits = Limits())

the state whose density matrix is the Hermitian part ``(\\rho + \\rho^\\dagger)/2`` of that of
`state`, the sum truncated by `limits`. A pure representation is returned as it is.
"""
hermitianize(state::State{Pure}; limits::Limits = Limits()) = state
hermitianize(state::State{Mixed}; limits::Limits = Limits()) = 0.5 * +(state, dag(state); limits)

"""
    hermiticity(::State)

how Hermitian the density matrix is, from 0 when it is anti-Hermitian to 1 when it is
Hermitian, as it should be:
``1/2 + \\mathrm{Re}\\,\\mathrm{tr}(\\rho^2)/(2\\,\\mathrm{tr}(\\rho^\\dagger\\rho))``. 1 on a
pure representation.
"""
hermiticity(::State{Pure}) = 1.
hermiticity(state::State{Mixed}) =
    0.5 + 0.5 * real(dot(state.state, dag(state).state)) / nonzero_norm(state, "hermiticity")^2

"""
    cut_svd(state, cut)

the factors `U` and `S` of the singular value decomposition of the state across the cut
between sites `cut` and `cut + 1`, `U` holding the sites up to `cut`; a cut the state does not
have is refused
"""
function cut_svd(state::State, cut::Int)
    n = length(state)
    if !(1 ≤ cut ≤ n)
        error("cannot compute the entanglement entropy at site $cut of a $n site state, the " *
              "cut is on the right of a site so it must be between 1 and $n")
    end
    nonzero_norm(state, "an entanglement entropy")
    s = orthogonalize(state.state, cut)
    U, S = svd(s[cut], (linkinds(s, cut - 1)..., siteinds(s, cut)...))
    return U, S
end

"""
    shannon_entropy(p)

the entropy ``-\\sum_i p_i \\log p_i`` of the probabilities `p`, a zero probability, which a
`mindim` above the Schmidt rank gives, adding nothing rather than NaN
"""
shannon_entropy(p) = -sum(x * log(x) for x in p if x > 0)

"""
    entropy_spectrum(S)

the entanglement entropy and spectrum read off the singular values `S` of a cut, as
`entanglement_entropy` returns them
"""
function entropy_spectrum(S::ITensor)
    # sorted: with conserved quantities the diagonal runs block by block
    sp = sort([ S[i, i]^2 for i in 1:dim(S, 1) ]; rev = true)
    sp /= sum(sp)
    return (shannon_entropy(sp), sp)
end

"""
    entanglement_entropy(state, cut::Int)

the entanglement entropy across the cut between sites `cut` and `cut + 1`, and the spectrum
it is computed from, the eigenvalues of the reduced density matrix of the sites up to `cut`
in decreasing order. `cut` runs from 1 to the number of sites, the last one giving 0. On a
mixed representation, this is the operator space entanglement entropy (OSEE).

# Examples

    ee, spectrum = entanglement_entropy(state, 3)   # cut between sites 3 and 4
"""
entanglement_entropy(state::State, cut::Int) = entropy_spectrum(cut_svd(state, cut)[2])

"""
    entanglement_by_sector(state::State{Pure}, cut::Int)

the entanglement across the cut between sites `cut` and `cut + 1`, resolved by the charge of
the sites up to `cut`: a `Dict` from each charge, a `QN` (`using ITensors: QN`), to
`(weight, entropy, spectrum)`, the probability of that charge and the entropy and decreasing
spectrum of the reduced density matrix restricted to it and normalised. The entanglement
entropy is `Σ weight * entropy - Σ weight * log(weight)`, the second sum being the number
entropy. A state conserving nothing has the single sector `QN()`; a mixed representation is
refused.

# Examples

    using ITensors: QN
    sectors = entanglement_by_sector(state, 3)   # cut between sites 3 and 4
    sectors[QN("N", 2)].weight
"""
function entanglement_by_sector(state::State{Pure}, cut::Int)
    U, S = cut_svd(state, cut)
    u = commonind(U, S)
    if !hasqns(u)
        ee, sp = entropy_spectrum(S)
        return Dict(QN() => (weight = 1.0, entropy = ee, spectrum = sp))
    end
    # U carrying no flux, what flows through `u` is what the sites up to `cut` hold
    if flux(U) ≠ QN()
        error("bug: the left factor of a cut carries $(flux(U)), so its sectors cannot be read")
    end
    squares = Dict{QN, Vector{Float64}}()
    k = 0
    for (q, d) in space(u)
        append!(get!(squares, dir(u) == ITensors.In ? q : -q, Float64[]),
                [ S[k + i, k + i]^2 for i in 1:d ])
        k += d
    end
    total = sum(sum, values(squares))
    sectors = Dict{QN, NamedTuple{(:weight, :entropy, :spectrum),
                                  Tuple{Float64, Float64, Vector{Float64}}}}()
    for (q, sq) in squares
        w = sum(sq)
        # a sector kept with nothing in it
        if w == 0
            continue
        end
        sp = sort(sq / w; rev = true)
        sectors[q] = (weight = w / total, entropy = shannon_entropy(sp), spectrum = sp)
    end
    return sectors
end

entanglement_by_sector(::State{Mixed}, ::Int) =
    error("entanglement_by_sector takes a pure state: on a mixed one the sectors of a link are " *
          "differences between the charges of the ket and of the bra")

"""
    renyi2(::State)
    renyi2(::State, positions::AbstractVector{Int})

the Rényi entropy of order 2, ``-\\log \\mathrm{tr}(\\rho^2)``, of the state, 0 on a pure
representation; given positions, that of the state reduced to those sites, 0 when there is
none. On a pure state, positions forming a block at an end, `1:k` or `k:n`, read the Schmidt
spectrum of the cut, as cheap as `entanglement_entropy`; any others mix the state first, which
is much more expensive.

# Examples

    renyi2(state, 1:3)
"""
renyi2(::State{Pure}) = 0.
renyi2(state::State{Mixed}) = -log(purity(state))

function renyi2(state::State, positions::AbstractVector{Int})
    check_positions(state, positions, "renyi2")
    kept = sort!(unique(collect(positions)))
    if isempty(kept)
        return 0.0
    elseif length(kept) == length(state)
        return renyi2(state)
    end
    return part_renyi2(state, kept)
end

# positions in a vector of another element type, `[]` or `Any[1, 3]`
renyi2(state::State, positions::AbstractVector) = renyi2(state, Vector{Int}(positions))

"""
    part_renyi2(state, kept)

the Rényi entropy of order 2 of the state reduced to the increasing positions `kept`, some of
the sites and not all: a block at an end of a pure state from the Schmidt spectrum of its cut,
read on the state carrying the signs of the fermionic sites traced out, whose plain reduction is
the fermionic one, see `trace_signs`
"""
part_renyi2(state::State{Mixed}, kept::Vector{Int}) =
    renyi2(partial_trace(weak_form(state), kept; keep = true))

function part_renyi2(state::State{Pure}, kept::Vector{Int})
    n = length(state)
    cut = kept == 1:last(kept) ? last(kept) : kept == first(kept):n ? first(kept) - 1 : nothing
    if isnothing(cut)
        return part_renyi2(mix(state), kept)
    end
    return -log(sum(abs2, last(entanglement_entropy(trace_signs(state, kept), cut))))
end

"""
    mutual_info_renyi2(state, positions::AbstractVector{Int})

the Rényi-2 analogue of the mutual information between two parts of the state,
``S_2(A) + S_2(B) - S_2(A \\cup B)``, which can be negative on a mixed state: part A at
`positions`, part B the rest. An empty part gives 0. On a pure state it is twice the `renyi2`
of part A, at the same cost, but on fermions, whose state may have no definite parity, where
it is the sum of the `renyi2` of the two parts.

# Examples

    mutual_info_renyi2(state, 1:3)
"""
function mutual_info_renyi2(state::State, positions::AbstractVector{Int})
    check_positions(state, positions, "mutual_info_renyi2")
    kept = sort!(unique(collect(positions)))
    if isempty(kept) || length(kept) == length(state)
        return 0.0
    end
    return parts_renyi2(state, kept)
end

# positions in a vector of another element type
mutual_info_renyi2(state::State, positions::AbstractVector) =
    mutual_info_renyi2(state, Vector{Int}(positions))

"""
    parts_renyi2(state, kept)

``S_2(A) + S_2(B) - S_2(A \\cup B)`` for part A at the increasing positions `kept` and B the
rest, neither of them empty; ``2 S_2(A)`` on a pure state, but on fermions, where the signs of
the sites traced out make ``S_2(B)`` another value on a state of no definite parity
"""
function parts_renyi2(state::State{Mixed}, kept::Vector{Int})
    w = weak_form(state)
    return renyi2(partial_trace(w, kept; keep = true)) + renyi2(partial_trace(w, kept)) - renyi2(w)
end

function parts_renyi2(state::State{Pure}, kept::Vector{Int})
    a = part_renyi2(state, kept)
    n = length(state)
    return next_fermion(state.system, 1) > n ? 2a : a + part_renyi2(state, setdiff(1:n, kept))
end

"""
    unit_trace(m)

the matrix `m` divided by its trace, refused when that trace is zero up to rounding, where a
density matrix of trace one cannot be formed
"""
function unit_trace(m::AbstractMatrix)
    t = tr(m)
    if abs(t) ≤ rounding_tol * norm(m)
        error("reduced_density_matrix needs a state of nonzero trace, its result having trace one")
    end
    return m / t
end

"""
    reduced_density_matrix(state, positions)

the density matrix of the sites at `positions`, the others traced out, as a matrix of trace
one on the basis of the product states of those sites, taken in their order, ordered as `kron`
orders them, the first site varying the slowest, as for `dense_state`. On fermionic sites it is
written in the Jordan-Wigner basis of the sites kept alone, so that `tr(ρ * matrix(op, ...))`
is the expectation value of `op` placed on them. It is meant for a few sites, the matrix
having the square of their dimension as its number of elements. A state of trace zero is
refused.

# Examples

    ρ = reduced_density_matrix(state, [2, 5])
"""
function reduced_density_matrix(state::State{Mixed}, positions::AbstractVector{Int})
    check_positions(state, positions, "reduced_density_matrix")
    if isempty(positions)
        return ones(1, 1)
    end
    r = partial_trace(weak_form(state), positions; keep = true)
    system = r.system
    t = prod(r.state)
    js, bs = Index[], Index[]
    # each mixed index split back into its ket and its bra
    for i in 1:length(system)
        j = siteind(Pure, system, i)
        b, c = mixer(j, siteind(Mixed, system, i), system[i])
        t *= dag(c)
        push!(js, j)
        push!(bs, b)
    end
    d = prod(dim, js)
    return unit_trace(reshape(Array(dense(t), reverse(js)..., reverse(dag.(prime.(bs)))...), d, d))
end

# the ket contracted with the bra, the sites kept left open on both sides
function reduced_density_matrix(state::State{Pure}, positions::AbstractVector{Int})
    check_positions(state, positions, "reduced_density_matrix")
    if isempty(positions)
        return ones(1, 1)
    end
    kept = sort(unique(positions))
    state = trace_signs(state, kept)
    st = state.state
    t = get_left(state, first(kept))
    for k in first(kept):last(kept)
        if k > first(kept)
            t *= st[k]
        end
        if !insorted(k, kept)
            t *= identity_at(state, k)
        end
        t *= k == last(kept) ? get_right(state, k) : dag(st[k]')
    end
    js = [ siteind(Pure, state.system, k) for k in kept ]
    d = prod(dim, js)
    return unit_trace(reshape(Array(dense(t), reverse(js)..., reverse(dag.(prime.(js)))...), d, d))
end

# positions in a vector of another element type, `[]` or `Any[1, 3]`
reduced_density_matrix(state::State, positions::AbstractVector) =
    reduced_density_matrix(state, Vector{Int}(positions))

"""
    vonneumann_entropy(state, positions)

the von Neumann entropy ``-\\mathrm{tr}\\,\\rho_A \\log \\rho_A`` of the sites at `positions`,
``\\rho_A`` being their `reduced_density_matrix`, 0 when there is none; the eigenvalues of
zero or below, which an approximate state may have, are left out. On a pure state and the
sites `1:cut`, `entanglement_entropy(state, cut)` gives the same, at a far lower cost and for
any number of sites.

# Examples

    vonneumann_entropy(state, [2, 5])
"""
function vonneumann_entropy(state::State, positions::AbstractVector)
    check_positions(state, positions, "vonneumann_entropy")
    λ = eigvals(Hermitian(reduced_density_matrix(state, positions)))
    return sum((-x * log(x) for x in λ if x > 0); init = 0.0)
end

"""
    log_negativity(state, a, b)

the logarithmic negativity ``\\log \\|\\rho^{T_B}\\|_1`` between the sites at `a` and those at
`b`, ``\\rho`` being their `reduced_density_matrix` and ``T_B`` the partial transpose on `b`,
zero for a separable state and 0 when a part is empty. The partial transpose in the
Jordan-Wigner basis is not the fermionic one: it is taken on a part holding no fermionic
site, and refused when both hold one.

# Examples

    log_negativity(state, [1, 2], [5])
"""
function log_negativity(state::State, a::AbstractVector, b::AbstractVector)
    check_positions(state, a, "log_negativity")
    check_positions(state, b, "log_negativity")
    a, b = sort(unique(Vector{Int}(a))), sort(unique(Vector{Int}(b)))
    if !isempty(intersect(a, b))
        error("log_negativity was given sites in both parts: $(intersect(a, b))")
    end
    if isempty(a) || isempty(b)
        return 0.0
    end
    system = state.system
    fermionic(part) = any(i -> next_fermion(system, i) == i, part)
    if fermionic(b)
        if fermionic(a)
            error("log_negativity needs a part without fermionic sites, whose partial " *
                  "transpose is the usual one")
        end
        a, b = b, a
    end
    u = sort([a; b])
    k = length(u)
    dims = [ dim(system[i]) for i in u ]
    # kets then bras, the last site first, as kron orders the basis
    t = reshape(reduced_density_matrix(state, u), reverse(dims)..., reverse(dims)...)
    perm = collect(1:2k)
    for (p, i) in enumerate(u)
        if i in b
            perm[k + 1 - p], perm[2k + 1 - p] = 2k + 1 - p, k + 1 - p
        end
    end
    d = prod(dims)
    return log(sum(abs, eigvals(Hermitian(reshape(permutedims(t, perm), d, d)))))
end

"""
    fidelity(a, b)

the fidelity of two states of the same system, between 0 and 1, whatever their norms and
traces: ``|\\langle a | b \\rangle|^2`` on two pure representations,
``\\langle \\psi | \\rho | \\psi \\rangle`` on a pure and a mixed one, in either order. Two
mixed representations are refused, their Uhlmann fidelity being out of reach for a matrix
product state: use `hs_fidelity`.

# Examples

    fidelity(state, ground_state)
"""
fidelity(a::State{Pure}, b::State{Pure}) =
    abs2(inner(a, b)) / (nonzero_trace(a, "fidelity") * nonzero_trace(b, "fidelity"))

fidelity(p::State{Pure}, r::State{Mixed}) =
    real(inner(mix(p), r)) / (nonzero_trace(p, "fidelity") * nonzero_trace(r, "fidelity"))

fidelity(r::State{Mixed}, p::State{Pure}) = fidelity(p, r)

fidelity(::State{Mixed}, ::State{Mixed}) =
    error("no fidelity between two mixed representations, it needs the spectrum of a " *
          "density operator. Use hs_fidelity")

"""
    hs_fidelity(a::State{Mixed}, b::State{Mixed})

the normalised Hilbert-Schmidt overlap of two mixed states of the same system,
``\\mathrm{tr}(ab)/\\sqrt{\\mathrm{tr}(a^2)\\mathrm{tr}(b^2)}``, 1 exactly when they are
proportional. This is **not** the Uhlmann fidelity, see `fidelity`.
"""
hs_fidelity(a::State{Mixed}, b::State{Mixed}) =
    real(inner(a, b)) / (nonzero_norm(a, "hs_fidelity") * nonzero_norm(b, "hs_fidelity"))

"""
    variance(state::State{Pure}, op)
    variance(state::State{Pure}, mpo::MPO)

the variance ``\\langle H^2\\rangle - \\langle H\\rangle^2`` of the operator `op`, or of the
operator of `mpo`, on the pure `state`, which need not be normalised, an operator changing the
charges by no definite amount being measured on the state conserving nothing, as by `expect`. A
mixed representation is refused.

# Examples

    variance(ground_state, hamiltonian)
"""
function variance(state::State{Pure}, mpo::MPO)
    st = state.state
    n2 = nonzero_trace(state, "variance")
    # the bra primed: unprimed, `inner` takes a deprecated fallback
    e = real(inner(prime(st), mpo, st)) / n2
    e2 = real(inner(mpo, st, mpo, st)) / n2
    return e2 - e^2
end

function variance(state::State{Pure}, op)
    mpo = try
        make_mpo(state, op)
    catch e
        e isa FluxError || rethrow()
        # an operator changing the charges by no definite amount is measured on the state
        # conserving nothing, as expect does
        return variance(plain_form(state), op)
    end
    return variance(state, mpo)
end

variance(::State{Mixed}, _) =
    error("variance needs a pure representation. On a mixed one steady_state returns " *
          "the equivalent for a Lindbladian")
