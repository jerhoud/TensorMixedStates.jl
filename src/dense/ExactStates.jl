# The exact representation of a state of a small system: its whole vector, or its whole density
# matrix, on the basis of the product states ordered as `kron` orders them, the first site
# varying the slowest, the operators acting through their terms, see `dense_acted`.

export ExactState

"""
    ExactState{R}(system, states; max_dense = nothing)
    ExactState(state::State; max_dense = nothing)
    ExactState(system, x)

the state of a small `system` held whole, without truncation: its vector for `R = Pure`, its
density matrix for `R = Mixed`. It is made from a description as `State{R}` takes it, from a
`State`, or from the whole vector or density matrix `x`, normalised, of definite charges on a
system conserving something, what rounding leaves in the sectors they forbid dropped. The
truncation limits and the methods that functions take for a `State` are accepted and ignored;
evolutions are exact up to `1e-10`, and a search stays in the sector of the charges declared of
the state it starts from, the symmetries not declared broken as DMRG does.
`random_state(ExactState{R}, …)`, with the arguments of `random_state(R, …)`, draws a pure one
uniformly in its sector, a mixed one from a uniform purification as for a `State`, whatever the
link dimension, and its `maxlinkdim` is that of the whole space it holds. Its memory grows as
``d^n``, ``d^{2n}`` for a density matrix, and a search or an evolution keeps tens of copies of
it, see the manual.

- `max_dense`: the number of elements of the vector or of the density matrix beyond which it is
  refused, none by default, as for `mix`; a whole array given is taken whatever its size

# Examples

    ExactState{Pure}(System(10, Qubit()), "Up")
    ExactState{Mixed}(System(4, Fermion()), FullyMixed(N => 2))
    ExactState(System(2, Qubit()), [1, 0, 0, 1])
"""
struct ExactState{R <: PM} <: AbstractState
    system::System
    state::VecOrMat{ComplexF64}
    # copied: the state never shares the array of the caller
    function ExactState{R}(system::System, x::AbstractVecOrMat, ::Val{:raw}) where {R <: PM}
        # a Float64, which does not wrap around as an Int does on many sites
        d = prod(float(dim(s)) for s in system.sites)
        if float.(size(x)) ≠ (R === Mixed ? (d, d) : (d,))
            n = Printf.@sprintf("%.0f", d)
            what = R === Mixed ? "a density matrix of $n × $n" : "a vector of $n elements"
            error("an ExactState{$R} of the system holds $what, not an array of size $(size(x))")
        end
        return new{R}(system, Array{ComplexF64}(x))
    end
end

ExactState{R}(args...; max_dense::Union{Nothing, Int} = nothing) where {R <: PM} =
    ExactState(State{R}(args...); max_dense)

# the limits, of the conversion between types of states, without use
ExactState(state::State{Pure}; max_dense::Union{Nothing, Int} = nothing,
           limits::Limits = Limits()) =
    ExactState{Pure}(state.system, Vector(state; max_dense), Val(:raw))

ExactState(state::State{Mixed}; max_dense::Union{Nothing, Int} = nothing,
           limits::Limits = Limits()) =
    ExactState{Mixed}(state.system, full_matrix(state; max_dense), Val(:raw))

function ExactState(system::System, x::AbstractVecOrMat{<:Number})
    R = x isa AbstractVector ? Pure : Mixed
    e = ExactState{R}(system, x, Val(:raw))
    if is_charged(system)
        a, idx, _ = dense_legs(system, x)
        check_definite_flux(a, idx)
        # the sectors the charges declared forbid are identically zero: what rounding leaves
        # there is dropped, as dense_state does
        e = ExactState{R}(system, reshape(in_charge(a, idx), size(x)), Val(:raw))
    end
    return normalize(e)
end

"""
    State(exact::ExactState; limits = Limits())

the MPS of the exact state `exact`, of its norm or its trace, truncated to `limits`.

# Examples

    State(ExactState{Pure}(System(4, Qubit()), "+"); limits = Limits(maxdim = 2))
"""
function State(state::ExactState{R}; limits::Limits = Limits()) where {R <: PM}
    a, idx, combiners = dense_legs(state.system, state.state)
    return State{R}(state.system, dense_mps(a, idx, siteinds(R, state.system), limits; combiners))
end

show(io::IO, s::ExactState{R}) where {R <: PM} = print(io, "ExactState{$R}($(s.system))")

representation(::ExactState{R}) where {R <: PM} = R

"""
    maxlinkdim(::ExactState)

the largest link dimension of an MPS of the system of the state, without truncation: that of
the whole space it holds

# Examples

    maxlinkdim(ExactState{Pure}(System(4, Qubit()), "Up"))      # 4
"""
function maxlinkdim(state::ExactState)
    ds = [ representation(state) === Mixed ? dim(s)^2 : dim(s) for s in state.system.sites ]
    return maximum((min(prod(ds[1:k]), prod(ds[k+1:end])) for k in 1:length(ds)-1); init = 1)
end

# the approximations of an MPS, without use for a state held whole
ignored_options(::ExactState) = (:algo, :limits, :mpo_limits, :randomize)

# as partial_trace of a State, the density matrix of the sites kept
function partial_trace(state::ExactState{Mixed};
                       keep::Union{Nothing, AbstractVector{Int}} = nothing,
                       trace_out::Union{Nothing, AbstractVector{Int}} = nothing)
    kept = kept_positions(state, keep, trace_out)
    return ExactState{Mixed}(System(state.system.sites[kept]), traced_matrix(state, kept),
                             Val(:raw))
end

partial_trace(::ExactState{Pure}; kwargs...) =
    error("partial_trace needs a mixed representation: use mix(state) first")

# for CreateState and ConvertState, see as_representation, the limits without use
as_representation(::Type{R}, state::ExactState{R}; limits::Limits = Limits()) where {R <: PM} =
    state
as_representation(::Type{Mixed}, state::ExactState{Pure}; limits::Limits = Limits()) = mix(state)
as_representation(::Type{Pure}, ::ExactState{Mixed}; limits::Limits = Limits()) =
    error("a mixed state cannot be made pure: ask for Mixed rather than Pure")

trace(state::ExactState{Pure}) = dot(state.state, state.state)
trace(state::ExactState{Mixed}) = tr(state.state)

norm(state::ExactState) = norm(state.state)

dense_bound(state::ExactState) = max(dense_elements, length(state.state))

nonzero_trace(state::ExactState{Mixed}, what) =
    checked_trace(real(trace(state)), norm(state), what)

function normalize(state::ExactState{R}) where {R <: PM}
    t = nonzero_trace(state, "normalize")
    x = R === Mixed ? state.state / t : state.state / sqrt(t)
    return ExactState{R}(state.system, x, Val(:raw))
end

function expect_value(state::ExactState{R}, a::IndexedOp{IsOperator}) where {R <: PM}
    # measured on what conserves nothing, as on a State, an operator of no definite flux too
    system = weaken(state.system, ())
    check_indices(system, a)
    x = state.state
    g = gathered(system, a)
    ts = developed(system, g)
    # on a density matrix, only the elements the terms reach, rather than the whole product
    v = R === Pure ? dot(x, first(gathered_acted(system, g, x))) :
        isnothing(ts) ? tr(first(gathered_acted(system, g, x))) : dense_trace(system, ts, x)
    return v / nonzero_trace(state, "expect")
end

mix(state::ExactState{Mixed}; limits::Limits = Limits(),
    max_dense::Union{Nothing, Int} = nothing) = state
function mix(state::ExactState{Pure}; limits::Limits = Limits(),
             max_dense::Union{Nothing, Int} = nothing)
    check_elements(float(length(state.state))^2,
                   "the density matrix of a state of $(length(state.system)) sites", max_dense)
    return ExactState{Mixed}(state.system, state.state * state.state', Val(:raw))
end

# as for a State, exact: the limits without use
(a::Number * b::ExactState{R}) where {R <: PM} = ExactState{R}(b.system, a * b.state, Val(:raw))
(a::ExactState * b::Number) = b * a
(a::ExactState / b::Number) = inv(b) * a
-(a::ExactState) = -1 * a

function +(a::ExactState{R}, b::ExactState{R}; limits::Limits = Limits()) where {R <: PM}
    if a.system ≠ b.system
        error("the two states are on different systems, $(system_difference(a.system, b.system))")
    end
    return ExactState{R}(a.system, a.state + b.state, Val(:raw))
end

-(a::ExactState{R}, b::ExactState{R}; limits::Limits = Limits()) where {R <: PM} =
    +(a, -b; limits)

truncate(state::ExactState; limits::Limits) = state

"""
    Vector(state::ExactState{Pure}; max_dense = nothing)

a copy of the vector of the exact pure `state`, as `Vector` of a `State` gives it; the density
matrix of a mixed state is given by `reduced_density_matrix(ρ, 1:n)`

- `max_dense`: the number of elements beyond which the vector is refused, none by default

# Examples

    Vector(ExactState{Pure}(System(2, Qubit()), "+"))           # [0.5, 0.5, 0.5, 0.5]
"""
function Vector(state::ExactState{Pure}; max_dense::Union{Nothing, Int} = nothing)
    check_elements(float(length(state.state)), "the vector of a state of " *
                   "$(length(state.system)) sites", max_dense)
    return copy(state.state)
end

dag(::ExactState{Pure}) = error("dag is meaningless on pure representations")
dag(state::ExactState{Mixed}) = ExactState{Mixed}(state.system, state.state', Val(:raw))

purity(::ExactState{Pure}) = 1.
purity(state::ExactState{Mixed}) = (norm(state) / nonzero_trace(state, "purity"))^2

hermiticity(::ExactState{Pure}) = 1.
hermiticity(state::ExactState{Mixed}) =
    0.5 + 0.5 * real(dot(state.state', state.state)) / nonzero_norm(state, "hermiticity")^2

hermitianize(state::ExactState{Pure}; limits::Limits = Limits()) = state
hermitianize(state::ExactState{Mixed}; limits::Limits = Limits()) =
    ExactState{Mixed}(state.system, (state.state + state.state') / 2, Val(:raw))

weaken(state::ExactState{R}, target::Conserved) where {R <: PM} =
    ExactState{R}(weaken(state.system, target), state.state, Val(:raw))
weaken(state::ExactState, spec) = weaken(state, Conserved(spec_names(spec)))
weaken(state::ExactState) = weaken(state, one_step_down(symmetries(state.system)))

"""
refuse two exact states of different sites, whatever they conserve
"""
function same_sites(a::ExactState, b::ExactState)
    if weaken(a.system, ()) ≠ weaken(b.system, ())
        error("the two states are on different systems, $(system_difference(a.system, b.system))")
    end
end

function inner(a::ExactState{R}, b::ExactState{R}) where {R <: PM}
    same_sites(a, b)
    return dot(a.state, b.state)
end

inner(::ExactState{Pure}, ::ExactState{Mixed}) =
    error("no inner product between a pure and a mixed representation: mix the pure one")
inner(::ExactState{Mixed}, ::ExactState{Pure}) =
    error("no inner product between a pure and a mixed representation: mix the pure one")

dot(a::ExactState, b::ExactState) = inner(a, b)

fidelity(a::ExactState{Pure}, b::ExactState{Pure}) =
    abs2(inner(a, b)) / (nonzero_trace(a, "fidelity") * nonzero_trace(b, "fidelity"))

function fidelity(p::ExactState{Pure}, r::ExactState{Mixed})
    same_sites(p, r)
    return real(dot(p.state, r.state * p.state)) /
           (nonzero_trace(p, "fidelity") * nonzero_trace(r, "fidelity"))
end

fidelity(r::ExactState{Mixed}, p::ExactState{Pure}) = fidelity(p, r)

"""
the square root of the Hermitian part of the density matrix of `state` of trace one, as `what`
needs, with a warning when it is not Hermitian beyond `approximation_tol`, refused when it has a
negative eigenvalue beyond it
"""
function density_sqrt(state::ExactState{Mixed}, what)
    ρ = state.state / nonzero_trace(state, what)
    h = (ρ + ρ') / 2
    δ = norm(ρ - h) / norm(h)
    if δ > approximation_tol
        @warn "$what takes the Hermitian part of a density matrix that is not Hermitian by " *
              "$(round(δ; sigdigits = 2))"
    end
    λ, u = eigen(Hermitian(h))
    if first(λ) < -approximation_tol * maximum(abs, λ)
        error("$what needs a positive density matrix, and the state has the eigenvalue " *
              "$(round(first(λ); sigdigits = 2)), as an approximation can leave it " *
              "(truncation, time step, unconverged search)")
    end
    # a defect τ of the state already makes the fidelity uncertain by about √τ, which zeroing a
    # negative eigenvalue above -τ does not exceed; the square root also makes the rounding of a
    # zero eigenvalue √ε: those within the noise of the eigensolver, about √D ε of the largest
    # (measured), are zero, which makes it exact on a pure state. The others are brought back to
    # a trace of one, which keeps the fidelity within [0, 1]
    noise = 3 * sqrt(length(λ)) * eps() * maximum(abs, λ)
    kept = [ x > noise ? x : 0.0 for x in λ ]
    return u * Diagonal(sqrt.(kept / sum(kept))) * u'
end

# the fidelity of Uhlmann and Jozsa, squared as |⟨a|b⟩|² is: (tr √(√a b √a))², the sum of the
# singular values of √a √b
function fidelity(a::ExactState{Mixed}, b::ExactState{Mixed})
    same_sites(a, b)
    return sum(svdvals(density_sqrt(a, "fidelity") * density_sqrt(b, "fidelity")))^2
end

hs_fidelity(a::ExactState{Mixed}, b::ExactState{Mixed}) =
    real(inner(a, b)) / (nonzero_norm(a, "hs_fidelity") * nonzero_norm(b, "hs_fidelity"))

# with a State, converted exactly, as the measurements of an exact state against an MPS ask
inner(a::ExactState, b::State) = inner(a, ExactState(b))
inner(a::State, b::ExactState) = inner(ExactState(a), b)
dot(a::ExactState, b::State) = inner(a, b)
dot(a::State, b::ExactState) = inner(a, b)
fidelity(a::ExactState, b::State) = fidelity(a, ExactState(b))
fidelity(a::State, b::ExactState) = fidelity(ExactState(a), b)
hs_fidelity(a::ExactState, b::State) = hs_fidelity(a, ExactState(b))
hs_fidelity(a::State, b::ExactState) = hs_fidelity(ExactState(a), b)

# as for a State, on what conserves nothing, not normalised
function inner(a::ExactState{Pure}, op::IndexedOp{IsOperator}, b::ExactState{Pure})
    same_sites(a, b)
    system = weaken(b.system, ())
    check_indices(system, op)
    return dot(a.state, operator_acted(system, op, b.state))
end

function inner(a::ExactState{Mixed}, op::IndexedOp{<:Union{IsChannel, IsGenerator}},
               b::ExactState{Mixed})
    same_sites(a, b)
    system = weaken(b.system, ())
    check_indices(system, op)
    return dot(a.state, dense_acted(Mixed, system, lower(system, op), b.state))
end

function inner(::ExactState{Mixed}, op::IndexedOp{IsOperator}, ::ExactState{Mixed})
    name = obs_name(op)
    error("between density matrices the operator $name is ambiguous: write Gate($name) for the " *
          "map ρ ↦ KρK† or Evolver($name) for ρ ↦ Kρ + ρK†")
end

dot(a::ExactState, op::IndexedOp, b::ExactState) = inner(a, op, b)
inner(a::ExactState, op::IndexedOp, b::State) = inner(a, op, ExactState(b))
inner(a::State, op::IndexedOp, b::ExactState) = inner(ExactState(a), op, b)
dot(a::ExactState, op::IndexedOp, b::State) = inner(a, op, b)
dot(a::State, op::IndexedOp, b::ExactState) = inner(a, op, b)

# as for a State, exact: its sites following those of `a`, the operators of `a` on the left on
# fermions, which the order of kron gives
⊗(a::ExactState{R}, b::ExactState{R}) where {R <: PM} =
    ExactState{R}(a.system ⊗ b.system, kron(a.state, b.state), Val(:raw))
⊗(a::ExactState{Pure}, b::ExactState{Mixed}) = mix(a) ⊗ b
⊗(a::ExactState{Mixed}, b::ExactState{Pure}) = a ⊗ mix(b)

# as for a State, exact, the limits without use: site `k` moved to position `p[k]`, each pair of
# sites `p` puts in the other order giving the sign of their exchange on fermions, on the kets
# and on the bras of a density matrix
function permute_sites(state::ExactState{R}, p::AbstractVector{<:Integer};
                       limits::Limits = Limits()) where {R <: PM}
    n = length(state)
    if sort(p) ≠ 1:n
        error("permute_sites needs a permutation of the positions 1 to $n, not $p")
    end
    sites = state.system.sites
    dims = site_dims(state.system)
    # for each basis state of the result, the one of `state` it comes from, and its sign
    q = invperm(p)
    from = in_kron(permutedims(on_sites(1:prod(dims), dims), q))
    signs = pair_signs(sites, [ (k, l) for k in 1:n for l in k+1:n if p[k] > p[l] ])
    sign = isnothing(signs) ? ones(prod(dims)) : signs[from]
    x = state.state
    y = R === Pure ? sign .* x[from] : (sign * sign') .* x[from, from]
    return ExactState{R}(System(sites[q]), y, Val(:raw))
end

"""
the matrix of `state` across the cut between sites `cut` and `cut + 1`: the sites up to `cut` on
the rows, the others on the columns, the ket and the bra of each together on a mixed state; a
cut the state does not have is refused
"""
function cut_matrix(state::ExactState{R}, cut::Int) where {R <: PM}
    system = state.system
    check_cut(length(system), cut)
    nonzero_norm(state, "an entanglement entropy")
    dims = site_dims(system)
    l, r = prod(dims[1:cut]), prod(dims[cut+1:end]; init = 1)
    if R === Pure
        return transpose(reshape(state.state, r, l))
    end
    # the axes of the kets then of the bras, the last site first, as kron orders the basis
    t = reshape(state.state, r, l, r, l)
    return reshape(permutedims(t, (2, 4, 1, 3)), l^2, r^2)
end

entanglement_entropy(state::ExactState, cut::Int) =
    first(entropy_spectrum(svdvals(cut_matrix(state, cut))))

entanglement_spectrum(state::ExactState, cut::Int) =
    last(entropy_spectrum(svdvals(cut_matrix(state, cut))))

function entanglement_by_sector(state::ExactState{Pure}, cut::Int)
    m = cut_matrix(state, cut)
    system = state.system
    if !is_charged(system)
        return entanglement_sectors(Tuple{String, Int}[], Dict(() => abs2.(svdvals(m))))
    end
    quantities = conserved_quantities(system)
    # the charge of each product state of the sites up to `cut`, in the order of `kron`
    qs = [QN()]
    for i in 1:cut
        site = reduce(vcat, [ fill(q, d) for (q, d) in space(siteind(Pure, system, i)) ])
        qs = [ q + s for q in qs for s in site ]
    end
    rows = Dict{SectorCharges, Vector{Int}}()
    for (k, q) in enumerate(qs)
        push!(get!(rows, sector_charges(q, quantities), Int[]), k)
    end
    return entanglement_sectors(quantities,
                                Dict(c => abs2.(svdvals(m[ks, :])) for (c, ks) in rows))
end

entanglement_by_sector(::ExactState{Mixed}, ::Int) =
    error("entanglement_by_sector takes a pure state: on a mixed one the sectors of a link are " *
          "differences between the charges of the ket and of the bra")

"""
the signs ``(-1)^{\\sum n_l n_k}`` of the product states of `system`, in the order of `kron`, `l`
running over the fermionic sites traced out and `k` over those of `kept` on their right, `n`
the parity of a site, which make the plain partial trace the fermionic one, see `trace_signs`;
`nothing` when there is none
"""
trace_sign_vector(system::System, kept::Vector{Int}) =
    pair_signs(system.sites, [ (l, k) for k in kept for l in 1:k-1 if !insorted(l, kept) ])

"""
the array of `state` with the signs of the fermionic sites traced out laid, see
`trace_sign_vector`, its axes in the order of `kron`, the last site first, the permutation of
these axes putting the sites at the increasing positions `kept` before the others, and the
dimensions of the two parts
"""
function kept_first(state::ExactState{R}, kept::Vector{Int}) where {R <: PM}
    system = state.system
    n = length(system)
    dims = site_dims(system)
    s = trace_sign_vector(system, kept)
    x = isnothing(s) ? state.state :
        R === Mixed ? s .* state.state .* transpose(s) : s .* state.state
    rest = setdiff(1:n, kept)
    perm = [ n + 1 - i for i in [reverse(kept); reverse(rest)] ]
    return x, reverse(dims), perm, prod(dims[kept]), prod(dims[rest]; init = 1)
end

"""
the vector of the pure `state` as a matrix ``M``, the sites at the increasing positions `kept`
on its rows, the others on its columns, the fermionic signs laid: their reduced density matrix
is ``M M^\\dagger``, and the squares of its singular values its spectrum
"""
function partition_matrix(state::ExactState{Pure}, kept::Vector{Int})
    x, rdims, perm, a, b = kept_first(state, kept)
    return reshape(permutedims(reshape(x, rdims...), perm), a, b)
end

function reduced_density_matrix(state::ExactState{R}, positions::AbstractVector{Int};
                                max_dense::Union{Nothing, Int} = dense_bound(state)) where {R <: PM}
    check_positions(state, positions, "reduced_density_matrix")
    check_increasing(positions, "reduced_density_matrix")
    if isempty(positions)
        return ones(1, 1)
    end
    kept = collect(positions)
    check_reduced(state, kept, max_dense)
    return unit_trace(traced_matrix(state, kept))
end

"""
the density matrix of the sites at the increasing positions `kept` of `state`, of the trace of
`state`
"""
function traced_matrix(state::ExactState{R}, kept::Vector{Int}) where {R <: PM}
    if R === Pure
        m = partition_matrix(state, kept)
        return m * m'
    end
    x, rdims, perm, a, b = kept_first(state, kept)
    n = length(rdims)
    t = reshape(permutedims(reshape(x, rdims..., rdims...), [perm; perm .+ n]), a, b, a, b)
    return sum(t[:, k, :, k] for k in 1:b)
end

"""
the spectrum of the reduced density matrix of trace one of the sites at the increasing
positions `kept` of the pure `state`, in decreasing order, from the singular values of
`partition_matrix` without forming that matrix; `what` names the caller in the refusals
"""
function kept_spectrum(state::ExactState{Pure}, kept::Vector{Int}, what)
    nonzero_trace(state, what)
    return last(entropy_spectrum(svdvals(partition_matrix(state, kept))))
end

# its reduced density matrix is never formed, which `max_dense` would bound
function von_neumann_entropy(state::ExactState{Pure}, positions::AbstractVector{Int};
                             max_dense::Union{Nothing, Int} = dense_bound(state))
    check_positions(state, positions, "von_neumann_entropy")
    if isempty(positions)
        return 0.0
    end
    kept = sort(Vector{Int}(positions))
    return shannon_entropy(kept_spectrum(state, kept, "von_neumann_entropy"))
end

part_renyi2(state::ExactState{Pure}, kept::Vector{Int}) =
    -log(sum(abs2, kept_spectrum(state, kept, "renyi2")))

"""
the operator `p`, placed on sites, applied to the kets of the density matrix `ρ`, then to its
bras as its adjoint: ``p \\rho p^\\dagger``
"""
function on_both_sides(system::System, p::IndexedOp{IsOperator}, ρ::AbstractMatrix)
    check_indices(system, p)
    g = gathered(system, p)
    return first(gathered_acted(system, g, first(gathered_acted(system, g, ρ))'))'
end

# applied as `apply` applies them to a `State`, the limits of the MPS without use
function apply(a::IndexedOp{IsOperator}, state::ExactState{Pure}; limits::Limits = Limits())
    check_indices(state.system, a)
    x = operator_acted(state.system, a, state.state)
    return ExactState{Pure}(state.system, x, Val(:raw))
end

function apply(a::IndexedOp{IsChannel}, state::ExactState{Pure}; limits::Limits = Limits())
    if !is_pure(a)
        error("$(obs_name(a)) has more than one Kraus operator, which a pure state cannot " *
              "take: apply it to mix(state)")
    end
    return apply(kraus(a), state)
end

apply(a::IndexedOp{IsGenerator}, ::ExactState{Pure}; limits::Limits = Limits()) =
    error("apply takes an operator or a channel, and $(obs_name(a)) is a generator")

apply(a::IndexedOp{IsOperator}, state::ExactState{Mixed}; limits::Limits = Limits()) =
    ExactState{Mixed}(state.system, on_both_sides(state.system, a, state.state), Val(:raw))

function apply(a::IndexedOp{<:Union{IsChannel, IsGenerator}}, state::ExactState{Mixed};
               limits::Limits = Limits())
    check_indices(state.system, a)
    x = dense_acted(Mixed, state.system, lower(state.system, a), state.state)
    return ExactState{Mixed}(state.system, x, Val(:raw))
end

apply(::MPO, ::ExactState; kwargs...) =
    error("an ExactState is given the operator rather than its MPO: apply(op, state)")

# the limits of the MPS have no use on an exact state, see the docstring of ExactState
function project_on(state::ExactState{R}, p::IndexedOp; limits::Limits) where {R <: PM}
    system = weaken(state.system, ())
    x = R === Mixed ? on_both_sides(system, p, state.state) :
                      operator_acted(system, p, state.state)
    return normalize(ExactState{R}(state.system, x, Val(:raw)))
end

"""
the vector `x`, of `extra` columns after the product states of `system`, projected on the total
`T` of the sum that `counting` describes: each site of the sum applies the projectors of its
values to each partial total, those that can still reach `T` kept, see `count_bonds`
"""
function count_projected(system::System, sites, locals, T::Int, x::AbstractVector, extra::Int)
    dims = site_dims(system)
    bonds = count_bonds(length(system), sites, locals, T)
    parts = Dict(0 => x)
    for (k, (values, ps)) in zip(sites, locals)
        a, before = around(dims, k)
        c = before * extra
        next = Dict{Int, Vector{ComplexF64}}()
        for (s, y) in parts, (t, p) in zip(values, ps)
            if insorted(s + t, bonds[k + 1])
                v = axis_acted(p, y, a, dims[k], c)
                next[s + t] = haskey(next, s + t) ? next[s + t] + v : v
            end
        end
        parts = next
    end
    return get(parts, T, zero(x))
end

function project_count(state::ExactState{R}, sites, locals, T; limits::Limits) where {R <: PM}
    system = state.system
    x = state.state
    if R === Pure
        y = count_projected(system, sites, locals, T, x, 1)
    else
        d = size(x, 1)
        kets = reshape(count_projected(system, sites, locals, T, vec(x), d), d, d)
        y = reshape(count_projected(system, sites, locals, T, vec(Matrix(kets')), d), d, d)'
    end
    return normalize(ExactState{R}(system, y, Val(:raw)))
end

"""
the probabilities of the product states of `state`, in the order of `kron`, which `draw` refuses
when one is negative beyond rounding
"""
basis_probabilities(state::ExactState{R}) where {R <: PM} =
    (R === Mixed ? real(diag(state.state)) : abs2.(state.state)) / nonzero_trace(state, "sample")

function sample(state::ExactState; rng = Random.default_rng())
    p = basis_probabilities(state)
    dims = site_dims(state.system)
    c = CartesianIndices(Tuple(reverse(dims)))[draw(x -> p[x + 1], length(p), rand(rng)) + 1]
    return [ c[length(dims) + 1 - i] - 1 for i in eachindex(dims) ]
end

function sample(state::ExactState, pos::Int; rng = Random.default_rng())
    check_positions(state, [pos], "sample")
    dims = site_dims(state.system)
    p = basis_probabilities(state)
    a, c = around(dims, pos)
    q = vec(sum(reshape(p, a, dims[pos], c); dims = (1, 3)))
    return draw(x -> q[x + 1], dims[pos], rand(rng))
end

function write_state(g, state::ExactState{R}) where {R <: PM}
    attributes(g)["representation"] = string(nameof(R))
    g["state"] = state.state
    return nothing
end

function read_state(::Type{<:ExactState}, g, sites)
    name = read(attributes(g)["representation"])
    R = name == "Pure" ? Pure : name == "Mixed" ? Mixed : nothing
    if isnothing(R)
        error("an ExactState of unknown representation \"$name\"")
    end
    return ExactState{R}(System(sites), read(g["state"]), Val(:raw))
end
