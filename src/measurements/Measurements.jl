# The measurements of a state that are not expectation values: its purity and its hermiticity,
# the entropies and the reduced density matrix of some of its sites, the negativity between two
# parts, the fidelity of two states and the variance of an operator.

export purity, hermiticity, hermitianize, renyi2, mutual_info_renyi2, entanglement_entropy,
       entanglement_spectrum, entanglement_by_sector, reduced_density_matrix, von_neumann_entropy,
       log_negativity,
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
refuse a cut a state of `n` sites does not have, `cut` being the site on its left
"""
function check_cut(n::Int, cut::Int)
    if !(1 ≤ cut ≤ n)
        error("cannot compute the entanglement entropy at site $cut of a $n site state, the " *
              "cut is on the right of a site so it must be between 1 and $n")
    end
end

"""
    cut_svd(state, cut)

the factors `U` and `S` of the singular value decomposition of the state across the cut
between sites `cut` and `cut + 1`, `U` holding the sites up to `cut`; a cut the state does not
have is refused
"""
function cut_svd(state::State, cut::Int)
    check_cut(length(state), cut)
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
shannon_entropy(p) = sum((-x * log(x) for x in p if x > 0); init = 0.0)

"""
    entropy_spectrum(s)

the entanglement entropy and spectrum read off the singular values `s` of a cut, a vector or the
diagonal `ITensor` of a decomposition, as `entanglement_entropy` returns them
"""
function entropy_spectrum(s::AbstractVector)
    # sorted: with conserved quantities the diagonal runs block by block
    sp = sort(abs2.(s); rev = true)
    sp /= sum(sp)
    return (shannon_entropy(sp), sp)
end

entropy_spectrum(S::ITensor) = entropy_spectrum([ S[i, i] for i in 1:dim(S, 1) ])

"""
    entanglement_entropy(state, cut::Int)

the entanglement entropy across the cut between sites `cut` and `cut + 1`, that of the
spectrum `entanglement_spectrum` gives. `cut` runs from 1 to the number of sites, the last one
giving 0. On a mixed representation, this is the operator space entanglement entropy (OSEE).

# Examples

    entanglement_entropy(state, 3)   # cut between sites 3 and 4
"""
entanglement_entropy(state::State, cut::Int) = first(entropy_spectrum(cut_svd(state, cut)[2]))

"""
    entanglement_spectrum(state, cut::Int)

the entanglement spectrum across the cut between sites `cut` and `cut + 1`, the eigenvalues of
the reduced density matrix of the sites up to `cut` in decreasing order, whose entropy
`entanglement_entropy` gives; on a mixed representation, the spectrum of the operator space
entanglement entropy

# Examples

    entanglement_spectrum(state, 3)   # cut between sites 3 and 4
"""
entanglement_spectrum(state::State, cut::Int) = last(entropy_spectrum(cut_svd(state, cut)[2]))

# the charges of a sector, as the pairs `name => charge` of the quantities conserved in the order
# of their names, and what is known of it
const SectorCharges = Tuple{Vararg{Pair{String, Int}}}
const Sector = NamedTuple{(:weight, :entropy, :spectrum), Tuple{Float64, Float64, Vector{Float64}}}

"""
    Sectors

the sectors of the entanglement across a cut, see `entanglement_by_sector`: a dictionary from
their charges to what is known of each, in the order of the charges, a sector being read by
its charges as `fully_mixed` takes them
"""
struct Sectors <: AbstractDict{SectorCharges, Sector}
    # the name and the modulus of each quantity conserved, in the order of their names
    quantities::Vector{Tuple{String, Int}}
    entries::Vector{Pair{SectorCharges, Sector}}
end

Base.length(s::Sectors) = length(s.entries)
Base.iterate(s::Sectors, k::Int = 1) = k > length(s.entries) ? nothing : (s.entries[k], k + 1)

"""
    sector_key(s, charges)

the charges `charges`, pairs `op => charge` in any order, `op` an operator or its name, as a key
of `s`, refused unless they give an integer charge to each quantity conserved
"""
function sector_key(s::Sectors, charges::Tuple)
    names = first.(s.quantities)
    given = Dict{String, Any}()
    for p in charges
        if !(p isa Pair)
            error("a sector is read by its charges, sectors[N => 2], not by $(repr(p))")
        end
        given[first(p) isa AbstractString ? String(first(p)) : obs_name(first(p))] = last(p)
    end
    if length(given) ≠ length(charges) || sort(collect(keys(given))) ≠ names
        error("a sector gives the charge of each quantity conserved, " *
              (isempty(names) ? "and nothing is" : join(names, ", ")) *
              ", got $(join(sort(collect(keys(given))), ", "))")
    end
    return Tuple(map(s.quantities) do (n, m)
        v = given[n]
        if !(v isa Real && isinteger(v))
            error("the charge of $n is an integer, not $v")
        end
        n => (m > 1 ? mod(Int(v), m) : Int(v))
    end)
end

function Base.get(s::Sectors, charges::Tuple, default)
    key = sector_key(s, charges)
    k = findfirst(e -> first(e) == key, s.entries)
    return isnothing(k) ? default : last(s.entries[k])
end

Base.get(s::Sectors, charges::Pair, default) = get(s, (charges,), default)

Base.get(::Sectors, ::QN, default) =
    error("a sector is read by its charges, as fully_mixed takes them, sectors[N => 2], where " *
          "1.x took a QN")

"""
    sector_charges(q, quantities)

the charges of the sector of charge `q`, the key of `Sectors` for the `quantities` conserved,
a quantity `q` does not hold having the charge zero
"""
function sector_charges(q::QN, quantities)
    vals = Dict(n => v for (n, v, _) in qn_components(q))
    return Tuple(n => (m > 1 ? mod(get(vals, n, 0), m) : get(vals, n, 0)) for (n, m) in quantities)
end

"""
the name and the modulus of each quantity the sites of `system` conserve, in the order of their
names, as `Sectors` holds them
"""
conserved_quantities(system::System) =
    sort(unique([ (q.name, q.modulus) for site in system.sites for q in conserved(site).quantities ]))

"""
the `Sectors` of the `quantities` conserved from the squares of the singular values of each
sector, `squares`, a sector of nothing in it left out
"""
function entanglement_sectors(quantities, squares::AbstractDict)
    total = sum(sum, values(squares))
    entries = Pair{SectorCharges, Sector}[]
    for (c, sq) in squares
        w = sum(sq)
        if w == 0
            continue
        end
        sp = sort(sq / w; rev = true)
        push!(entries, c => (weight = w / total, entropy = shannon_entropy(sp), spectrum = sp))
    end
    return Sectors(quantities, sort!(entries; by = first))
end

"""
    entanglement_by_sector(state::State{Pure}, cut::Int)

the entanglement across the cut between sites `cut` and `cut + 1`, resolved by the charge of
the sites up to `cut`: a dictionary from the charges of each sector, the pairs `name => charge`
of the quantities conserved in the order of their names, to `(weight, entropy, spectrum)`, the
probability of that charge and the entropy and decreasing spectrum of the reduced density
matrix restricted to it and normalised. A sector is read by its charges as `fully_mixed` takes
them, in any order, `sectors[N => 2]` or `sectors[Ntot => 4, 2Sz => 0]`. The entanglement
entropy is `Σ weight * entropy - Σ weight * log(weight)`, the second sum being the number
entropy. A state conserving nothing has the single sector `()`; a mixed representation is
refused.

# Examples

    sectors = entanglement_by_sector(state, 3)   # cut between sites 3 and 4
    sectors[N => 2].weight
"""
function entanglement_by_sector(state::State{Pure}, cut::Int)
    U, S = cut_svd(state, cut)
    u = commonind(U, S)
    if !hasqns(u)
        return entanglement_sectors(Tuple{String, Int}[],
                                    Dict(() => [ S[i, i]^2 for i in 1:dim(S, 1) ]))
    end
    # U carrying no flux, what flows through `u` is what the sites up to `cut` hold
    if flux(U) ≠ QN()
        error("bug: the left factor of a cut carries $(flux(U)), so its sectors cannot be read")
    end
    quantities = conserved_quantities(state.system)
    squares = Dict{SectorCharges, Vector{Float64}}()
    k = 0
    for (q, d) in space(u)
        c = sector_charges(dir(u) == ITensors.In ? q : -q, quantities)
        append!(get!(squares, c, Float64[]), [ S[k + i, k + i]^2 for i in 1:d ])
        k += d
    end
    return entanglement_sectors(quantities, squares)
end

entanglement_by_sector(::State{Mixed}, ::Int) =
    error("entanglement_by_sector takes a pure state: on a mixed one the sectors of a link are " *
          "differences between the charges of the ket and of the bra")

"""
    renyi2(state)
    renyi2(state, positions::AbstractVector{Int})

the Rényi entropy of order 2, ``-\\log \\mathrm{tr}(\\rho^2)``, of the state, 0 on a pure
representation; given positions, that of the state reduced to those sites, 0 when there is
none. On a pure state, positions forming a block at an end, `1:k` or `k:n`, read the Schmidt
spectrum of the cut, as cheap as `entanglement_entropy`; any others mix the state first, which
is much more expensive.

# Examples

    renyi2(state, 1:3)
"""
renyi2(::State{Pure}) = 0.
renyi2(state::AbstractState) = -log(purity(state))

function renyi2(state::AbstractState, positions::AbstractVector{Int})
    check_positions(state, positions, "renyi2")
    kept = sort!(collect(positions))
    if isempty(kept)
        return 0.0
    elseif length(kept) == length(state)
        return renyi2(state)
    end
    return part_renyi2(state, kept)
end

# positions in a vector of another element type, `[]` or `Any[1, 3]`
renyi2(state::AbstractState, positions::AbstractVector) = renyi2(state, Vector{Int}(positions))

"""
    part_renyi2(state, kept)

the Rényi entropy of order 2 of the state reduced to the increasing positions `kept`, some of
the sites and not all: a block at an end of a pure state from the Schmidt spectrum of its cut,
read on the state carrying the signs of the fermionic sites traced out, whose plain reduction is
the fermionic one, see `trace_signs`; on another representation, from its
`reduced_density_matrix`, of trace one
"""
part_renyi2(state::AbstractState, kept::Vector{Int}) =
    -log(sum(abs2, reduced_density_matrix(state, kept)))

part_renyi2(state::State{Mixed}, kept::Vector{Int}) =
    renyi2(partial_trace(weak_form(state); keep = kept))

function part_renyi2(state::State{Pure}, kept::Vector{Int})
    n = length(state)
    cut = kept == 1:last(kept) ? last(kept) : kept == first(kept):n ? first(kept) - 1 : nothing
    if isnothing(cut)
        return part_renyi2(mix(state), kept)
    end
    return -log(sum(abs2, entanglement_spectrum(trace_signs(state, kept), cut)))
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
function mutual_info_renyi2(state::AbstractState, positions::AbstractVector{Int})
    check_positions(state, positions, "mutual_info_renyi2")
    kept = sort!(collect(positions))
    if isempty(kept) || length(kept) == length(state)
        return 0.0
    end
    return parts_renyi2(state, kept)
end

# positions in a vector of another element type
mutual_info_renyi2(state::AbstractState, positions::AbstractVector) =
    mutual_info_renyi2(state, Vector{Int}(positions))

"""
    parts_renyi2(state, kept)

``S_2(A) + S_2(B) - S_2(A \\cup B)`` for part A at the increasing positions `kept` and B the
rest, neither of them empty; ``2 S_2(A)`` on a pure state, but on fermions, where the signs of
the sites traced out make ``S_2(B)`` another value on a state of no definite parity
"""
parts_renyi2(state::AbstractState, kept::Vector{Int}) =
    part_renyi2(state, kept) + part_renyi2(state, setdiff(1:length(state), kept)) - renyi2(state)

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
    check_reduced(state, kept)

refuse the density matrix of the sites at the positions `kept` of `state` when its dimension,
the product of theirs, is beyond `dense_max`, rather than exhausting memory
"""
function check_reduced(state::AbstractState, kept)
    # a Float64, which does not wrap around as an Int does on many sites
    d = prod(float(dim(state.system[k])) for k in kept)
    if d > dense_max
        error("reduced_density_matrix of $(length(kept)) sites needs a dense matrix of dimension " *
              "$(Printf.@sprintf("%.0f", d)), beyond the $dense_max TMS builds")
    end
end

"""
    reduced_density_matrix(state, positions)

the density matrix of the sites at `positions`, the others traced out, as a matrix of trace
one on the basis of the product states of those sites, given in increasing order, ordered as
`kron` orders them, the first site varying the slowest, as for `dense_state`. On fermionic sites
it is written in the Jordan-Wigner basis of the sites kept alone, so that
`tr(ρ * matrix(op, ...))` is the expectation value of `op` placed on them. It is meant for a few
sites, the matrix having the square of their dimension as its number of elements, and refused
beyond the dimension `dense_max`, 4096. A state of trace zero is refused, and so are positions
out of increasing order, which `permute_sites` gives.

# Examples

    ρ = reduced_density_matrix(state, [2, 5])
"""
function reduced_density_matrix(state::State{Mixed}, positions::AbstractVector{Int})
    check_positions(state, positions, "reduced_density_matrix")
    check_increasing(positions, "reduced_density_matrix")
    if isempty(positions)
        return ones(1, 1)
    end
    check_reduced(state, positions)
    return unit_trace(full_matrix(partial_trace(weak_form(state); keep = positions)))
end

"""
the environment of the sites up to `k` of the pure `state`, ket with bra, as a matrix on the link
after site `k`, of dimension 1 for `k = 0`
"""
function left_link_matrix(state::State{Pure}, k::Int)
    if k == 0
        return ones(1, 1)
    end
    st = state.state
    link = commonind(st[k], st[k+1])
    return Array(dense(get_left(state, k) * identity_at(state, k) * dag(st[k]')), link, link')
end

"""
the environment of the sites after `k` of the pure `state`, ket with bra, as a matrix on the
link after site `k`, of dimension 1 for the last site
"""
function right_link_matrix(state::State{Pure}, k::Int)
    st = state.state
    if k == length(st)
        return ones(1, 1)
    end
    link = commonind(st[k], st[k+1])
    return Array(dense(get_right(state, k + 1) * identity_at(state, k + 1) * st[k+1]), link, link')
end

"""
the elements of the density matrix of the sites `kept` of the pure `state`, in the order of
`mixed_order`, ket and bra contracted site by site, those of the sites kept left open
"""
function kept_pairs(state::State{Pure}, kept::Vector{Int})
    st = state.state
    e = left_link_matrix(state, first(kept) - 1)
    # the ket link, the bra link, the pairs of the sites kept so far
    v = reshape(e, size(e)..., 1)
    for k in first(kept):last(kept)
        a = mps_array(st, k, siteind(Pure, state.system, k))
        χl, d, χr = size(a)
        np = size(v, 3)
        # the ket of the site taken: its value, the ket link, the bra link, the pairs
        x = reshape(transpose(reshape(a, χl, :)) * reshape(v, χl, :), d, χr, χl, np)
        if insorted(k, kept)
            y = transpose(conj(reshape(a, χl, :))) * reshape(permutedims(x, (3, 1, 2, 4)), χl, :)
            # the pair of the site the fastest of the pairs, its ket faster
            v = reshape(permutedims(reshape(y, d, χr, d, χr, np), (4, 2, 3, 1, 5)), χr, χr, :)
        else
            y = transpose(conj(reshape(a, χl * d, χr))) *
                reshape(permutedims(x, (3, 1, 2, 4)), χl * d, :)
            v = permutedims(reshape(y, χr, χr, np), (2, 1, 3))
        end
    end
    e = right_link_matrix(state, last(kept))
    return transpose(reshape(v, :, size(v, 3))) * vec(e)
end

"""
the density matrix, not normalised, of the consecutive sites `ks` of the pure `state`: ``M B``,
`M` the block of the kets, its links open, and `B` that of the bras with the environments
"""
function block_matrix(state::State{Pure}, ks::UnitRange{Int})
    el = left_link_matrix(state, first(ks) - 1)
    er = right_link_matrix(state, last(ks))
    χl, χr = size(el, 1), size(er, 1)
    m = contracted(state.state, ks, [ siteind(Pure, state.system, k) for k in ks ];
                   start = Matrix{Float64}(I, χl, χl))
    d = length(m) ÷ (χl * χr)
    m = reshape(m, d, χl, χr)
    # the bras with the environment on the left, then with that on the right
    b = el * reshape(permutedims(conj(m), (2, 1, 3)), χl, :)
    b = reshape(reshape(b, χl * d, χr) * transpose(er), χl, d, χr)
    return reshape(m, d, :) * transpose(reshape(permutedims(b, (2, 1, 3)), d, :))
end

function reduced_density_matrix(state::State{Pure}, positions::AbstractVector{Int})
    check_positions(state, positions, "reduced_density_matrix")
    check_increasing(positions, "reduced_density_matrix")
    if isempty(positions)
        return ones(1, 1)
    end
    kept = collect(positions)
    check_reduced(state, kept)
    state = trace_signs(state, kept)
    if kept == first(kept):last(kept)
        return unit_trace(block_matrix(state, first(kept):last(kept)))
    end
    return unit_trace(pairs_matrix(kept_pairs(state, kept), [ dim(state.system[k]) for k in kept ]))
end

# positions in a vector of another element type, `[]` or `Any[1, 3]`; integer positions reaching
# this method have no method of their own, which converting them again would not find
function reduced_density_matrix(state::AbstractState, positions::AbstractVector)
    if positions isa AbstractVector{Int}
        throw(MethodError(reduced_density_matrix, (state, positions)))
    end
    return reduced_density_matrix(state, Vector{Int}(positions))
end

"""
    von_neumann_entropy(state, positions)

the von Neumann entropy ``-\\mathrm{tr}\\,\\rho_A \\log \\rho_A`` of the sites at `positions`,
``\\rho_A`` being their `reduced_density_matrix`, 0 when there is none; a negative eigenvalue,
which an approximate state may have, is taken as zero, and refused beyond `1e-8` of the largest. On a pure state and the
sites `1:cut`, `entanglement_entropy(state, cut)` gives the same, at a far lower cost and for
any number of sites.

# Examples

    von_neumann_entropy(state, [2, 5])
"""
function von_neumann_entropy(state::AbstractState, positions::AbstractVector)
    check_positions(state, positions, "von_neumann_entropy")
    if isempty(positions)
        return 0.0
    end
    λ = eigvals(Hermitian(reduced_density_matrix(state, sort(Vector{Int}(positions)))))
    if first(λ) < -approximation_tol * maximum(abs, λ)
        error("the reduced density matrix has the negative eigenvalue " *
              "$(round(first(λ); sigdigits = 2)): the state is not a positive density matrix, as " *
              "an approximation can leave it (truncation, time step, unconverged search)")
    end
    return shannon_entropy(λ)
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
function log_negativity(state::AbstractState, a::AbstractVector, b::AbstractVector)
    check_positions(state, a, "log_negativity")
    check_positions(state, b, "log_negativity")
    a, b = sort(Vector{Int}(a)), sort(Vector{Int}(b))
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
    pure_on_mixed(p, r)

``\\langle p|r|p\\rangle`` for the pure state `p` and the density matrix `r` of one system,
contracted site by site with the tensors of `mix(p)` conjugated, their links left apart: the
contraction holds ``\\chi_p^2\\chi_r`` numbers where the MPS of `mix(p)` has bonds of ``\\chi_p^2``
"""
function pure_on_mixed(p::State{Pure}, r::State{Mixed})
    sys = r.system
    starred = relabeller(i -> star(i, strong_names(sys)))
    st = hasqns(p.state) ? p.state : dense(p.state)
    t = ITensor(1.)
    for i in 1:length(sys)
        j = siteind(Pure, sys, i)
        # the combiner of `mix`, the bra starred whole
        c = combinerto(siteind(Mixed, sys, i), j, dag(starred(j')))
        t = t * r.state[i] * dag(c) * dag(st[i]) * relabel(st[i]', starred)
    end
    return scalar(t)
end

"""
    fidelity(a, b)

the fidelity of two states of the same sites, taken on what both conserve, between 0 and 1,
whatever their norms and traces: ``|\\langle a | b \\rangle|^2`` on two pure representations,
``\\langle \\psi | \\rho | \\psi \\rangle`` on a pure and a mixed one, in either order. Two
mixed representations are refused, their Uhlmann fidelity being out of reach for a matrix
product state: use `hs_fidelity`.

# Examples

    fidelity(state, ground_state)
"""
fidelity(a::State{Pure}, b::State{Pure}) =
    abs2(inner(a, b)) / (nonzero_trace(a, "fidelity") * nonzero_trace(b, "fidelity"))

function fidelity(p::State{Pure}, r::State{Mixed})
    p, r = on_common_system(p, r)
    return real(pure_on_mixed(p, r)) /
           (nonzero_trace(p, "fidelity") * nonzero_trace(r, "fidelity"))
end

fidelity(r::State{Mixed}, p::State{Pure}) = fidelity(p, r)

fidelity(::State{Mixed}, ::State{Mixed}) =
    error("no fidelity between two mixed representations, it needs the spectrum of a " *
          "density operator. Use hs_fidelity")

"""
    hs_fidelity(a::State{Mixed}, b::State{Mixed})

the normalised Hilbert-Schmidt overlap of two mixed states of the same sites, taken on what
both conserve, ``\\mathrm{tr}(ab)/\\sqrt{\\mathrm{tr}(a^2)\\mathrm{tr}(b^2)}``, 1 exactly when
they are proportional. This is **not** the Uhlmann fidelity, see `fidelity`.
"""
hs_fidelity(a::State{Mixed}, b::State{Mixed}) =
    real(inner(a, b)) / (nonzero_norm(a, "hs_fidelity") * nonzero_norm(b, "hs_fidelity"))

"""
    variance(state, op)
    variance(state::State{Pure}, mpo::MPO)

the variance ``\\langle H^\\dagger H\\rangle - |\\langle H\\rangle|^2`` of the operator `op`,
placed on sites, ``\\langle H^2\\rangle - \\langle H\\rangle^2`` for a Hermitian one, or of
the operator of `mpo`, on `state`, which need not be normalised, an operator changing the
charges by no definite amount being measured on the state conserving nothing, as by `expect`.
On a density matrix ``\\rho`` of trace one, it is
``\\mathrm{tr}(\\rho H^\\dagger H) - |\\mathrm{tr}(\\rho H)|^2``.

# Examples

    variance(ground_state, hamiltonian)
    β^2 * variance(thermal_state, hamiltonian)      # the specific heat
"""
function variance(state::State{Pure}, mpo::MPO)
    st = state.state
    n2 = nonzero_trace(state, "variance")
    # the bra primed: unprimed, `inner` takes a deprecated fallback
    e = inner(prime(st), mpo, st) / n2
    e2 = real(inner(mpo, st, mpo, st)) / n2
    return e2 - abs2(e)
end

variance(state::State{Pure}, op::IndexedOp{IsOperator}) =
    or_on_plain(st -> variance(st, make_mpo(st, op)), state)

# `op` of any type in the signature, so that a representation defining `variance(::MyState, op)`
# stays unambiguous, and checked here
function variance(state::AbstractState, op)
    if !(op isa IndexedOp{IsOperator})
        error("variance takes an operator placed on sites, as X(1) * X(2), or an MPO on a pure " *
              "state")
    end
    return expect(state, op' * op) - abs2(expect(state, op; complex = true))
end
