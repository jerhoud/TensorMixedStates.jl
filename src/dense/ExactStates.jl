# The exact representation of a state of a small system: its whole vector, or its whole density
# matrix, on the basis of the product states ordered as `kron` orders them, the first site
# varying the slowest, the operators acting through their terms, see `dense_acted`.

"""
    ExactState{R}(system, states; max_dense = nothing)
    ExactState(state::State; max_dense = nothing)
    ExactState(system, x)

the state of a small `system` held whole, without truncation: its vector for `R = Pure`, its
density matrix for `R = Mixed`. It is made from a description as `State{R}` takes it, from a
`State`, or from the whole vector or density matrix `x`, normalised, of definite charges on a
system conserving something. The truncation limits and the methods that functions take for a
`State` are accepted and ignored.

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

ExactState(state::State{Pure}; max_dense::Union{Nothing, Int} = nothing) =
    ExactState{Pure}(state.system, Vector(state; max_dense), Val(:raw))

ExactState(state::State{Mixed}; max_dense::Union{Nothing, Int} = nothing) =
    ExactState{Mixed}(state.system, full_matrix(state; max_dense), Val(:raw))

function ExactState(system::System, x::AbstractVecOrMat{<:Number})
    e = ExactState{x isa AbstractVector ? Pure : Mixed}(system, x, Val(:raw))
    if is_charged(system)
        a, idx, _ = dense_legs(system, x)
        check_definite_flux(a, idx)
    end
    return normalize(e)
end

"""
    State(exact::ExactState; limits = Limits())

the MPS of the exact state `exact`, normalised, truncated to `limits`.
"""
State(state::ExactState; limits::Limits = Limits()) =
    dense_state(state.system, state.state; limits)

show(io::IO, s::ExactState{R}) where R = print(io, "ExactState{$R}($(s.system))")

trace(state::ExactState{Pure}) = dot(state.state, state.state)
trace(state::ExactState{Mixed}) = tr(state.state)

norm(state::ExactState) = norm(state.state)

dense_bound(state::ExactState) = max(dense_elements, length(state.state))

nonzero_trace(state::ExactState{Mixed}, what) =
    checked_trace(real(trace(state)), norm(state), what)

function normalize(state::ExactState{R}) where R
    t = nonzero_trace(state, "normalize")
    x = R === Mixed ? state.state / t : state.state / sqrt(t)
    return ExactState{R}(state.system, x, Val(:raw))
end

function expect_value(state::ExactState{R}, a::IndexedOp{IsOperator}) where R
    # measured on what conserves nothing, as on a State, an operator of no definite flux too
    system = weaken(state.system, ())
    x = state.state
    ax = dense_acted(Pure, system, lower(system, a), x)
    return (R === Mixed ? tr(ax) : dot(x, ax)) / nonzero_trace(state, "expect")
end

mix(state::ExactState{Mixed}; limits::Limits = Limits(),
    max_dense::Union{Nothing, Int} = nothing) = state
function mix(state::ExactState{Pure}; limits::Limits = Limits(),
             max_dense::Union{Nothing, Int} = nothing)
    check_elements(float(length(state.state))^2,
                   "the density matrix of a state of $(length(state.system)) sites", max_dense)
    return ExactState{Mixed}(state.system, state.state * state.state', Val(:raw))
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

function inner(a::ExactState{R}, b::ExactState{R}) where R
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

"""
the matrix of `state` across the cut between sites `cut` and `cut + 1`: the sites up to `cut` on
the rows, the others on the columns, the ket and the bra of each together on a mixed state; a
cut the state does not have is refused
"""
function cut_matrix(state::ExactState{R}, cut::Int) where R
    system = state.system
    check_cut(length(system), cut)
    nonzero_norm(state, "an entanglement entropy")
    dims = [ dim(s) for s in system.sites ]
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
function trace_sign_vector(system::System, kept::Vector{Int})
    n = length(system)
    dims = [ dim(s) for s in system.sites ]
    fermionic = [ next_fermion(system, i) == i for i in 1:n ]
    last_kept = findlast(k -> fermionic[k], kept)
    if isnothing(last_kept) || !any(l -> fermionic[l] && !insorted(l, kept), 1:kept[last_kept]-1)
        return nothing
    end
    odd = [ isodd.(collect(fermion_parities(s))) for s in system.sites ]
    s = ones(prod(dims))
    for (j, c) in enumerate(CartesianIndices(Tuple(reverse(dims))))
        left, flip = false, false
        for i in 1:n
            if fermionic[i]
                p = odd[i][c[n + 1 - i]]
                if insorted(i, kept)
                    flip ⊻= left & p
                else
                    left ⊻= p
                end
            end
        end
        if flip
            s[j] = -1
        end
    end
    return s
end

"""
the array of `state` with the signs of the fermionic sites traced out laid, see
`trace_sign_vector`, its axes in the order of `kron`, the last site first, the permutation of
these axes putting the sites at the increasing positions `kept` before the others, and the
dimensions of the two parts
"""
function kept_first(state::ExactState{R}, kept::Vector{Int}) where R
    system = state.system
    n = length(system)
    dims = [ dim(s) for s in system.sites ]
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
                                max_dense::Union{Nothing, Int} = dense_bound(state)) where R
    check_positions(state, positions, "reduced_density_matrix")
    check_increasing(positions, "reduced_density_matrix")
    if isempty(positions)
        return ones(1, 1)
    end
    kept = collect(positions)
    check_reduced(state, kept, max_dense)
    if R === Pure
        m = partition_matrix(state, kept)
        return unit_trace(m * m')
    end
    x, rdims, perm, a, b = kept_first(state, kept)
    n = length(rdims)
    t = reshape(permutedims(reshape(x, rdims..., rdims...), [perm; perm .+ n]), a, b, a, b)
    return unit_trace(sum(t[:, k, :, k] for k in 1:b))
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
function von_neumann_entropy(state::ExactState{Pure}, positions::AbstractVector;
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
function on_both_sides(system::System, p::IndexedOp, ρ::AbstractMatrix)
    ts = lower(system, p)
    return dense_acted(Pure, system, ts, dense_acted(Pure, system, ts, ρ)')'
end

# applied as `apply` applies them to a `State`, the limits of the MPS without use
function apply(a::IndexedOp{IsOperator}, state::ExactState{Pure}; limits::Limits = Limits())
    check_indices(state.system, a)
    x = dense_acted(Pure, state.system, lower(state.system, a), state.state)
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
    apply(Gate(a), state)

function apply(a::IndexedOp{<:Union{IsChannel, IsGenerator}}, state::ExactState{Mixed};
               limits::Limits = Limits())
    check_indices(state.system, a)
    x = dense_acted(Mixed, state.system, lower(state.system, a), state.state)
    return ExactState{Mixed}(state.system, x, Val(:raw))
end

apply(::MPO, ::ExactState; kwargs...) =
    error("an ExactState is given the operator rather than its MPO: apply(op, state)")

"""
    ExactEvolver

the evolver of an exact state: the terms of each of its parts, see `lower`, their time
functions, `nothing` when it does not depend on time, see `split_terms`, and the duration of the
whole evolution, which bounds its error, see `magnus_evolved`
"""
struct ExactEvolver
    terms::Vector{Vector{Term}}
    functions::Union{Nothing, Vector}
    duration::Union{Nothing, Number}
end

prepared_evolver(state::ExactState, a::Union{IndexedOp, Pair}, mpo_limits::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    prepared_evolver(state, [a], mpo_limits; duration)

# a pure state follows the operator, a density matrix the evolver of an operator
function prepared_evolver(state::ExactState{R}, terms::AbstractVector, ::Limits = Limits();
                          duration::Union{Nothing, Number} = nothing) where {R <: PM}
    system = state.system
    ops, fs = split_terms(terms)
    parts = map(ops) do a
        check_indices(system, a)
        return lower(system, R === Pure ? followed_operator(a) : superoperator(a))
    end
    return ExactEvolver(parts, fs, duration)
end

prepared_evolver(::ExactState, ::Union{PreMPO, MPO}, ::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    error("an ExactState is given the evolver rather than its PreMPO or MPO: " *
          "evolve(evolver, duration, state)")

prepared_evolver(::ExactState, x, ::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) = error("$(repr(x)) is not an evolver")

"""
the exponential ``e^{tA}x`` of the linear map `A`, without its matrix, by Krylov, to
`rounding_tol` of the norm of `x`, refused when it does not converge
"""
function exact_exp(A, t::Number, x::AbstractArray)
    n = norm(x)
    if iszero(n)
        return x
    end
    # KrylovKit takes its precision as absolute: the vector of norm one, its norm put back after
    y, info = KrylovKit.exponentiate(A, t, x / n; tol = rounding_tol, ishermitian = false,
                                     eager = true)
    if info.converged == 0
        error("the exponential of the evolver of an ExactState did not converge in Krylov to " *
              "$rounding_tol: take smaller steps")
    end
    return n * y
end

"""
the array `x` evolved from `time` for `dt` under the evolver `cs -> x -> Σ_j c_j A_j x` that
`act` gives for the values `cs` of its time functions `fs`, by the commutator-free Magnus
integrator of order 4, see `step_coefs`, each exponential by Krylov, in substeps whose error,
measured against two half substeps, keeps that of an evolution of `duration` below
`1000 rounding_tol`
"""
function magnus_evolved(act, fs, time::Number, dt::Number, duration::Number, x::AbstractArray)
    if iszero(norm(x))
        return x
    end
    magnus(t, h, y) = foldl((y, cs) -> exact_exp(act(cs), h, y), step_coefs(fs, t, h, 4);
                            init = y)
    # the substeps as fractions of `dt`, from `s` and of length `f`
    s, f = 0.0, 1.0
    while 1 - s > 4 * eps()
        f = min(f, 1 - s)
        t, h = time + s * dt, f * dt
        whole = magnus(t, h, x)
        halves = magnus(t + h / 2, h / 2, magnus(t, h / 2, x))
        # the local error of order 5 in the length of the substep: the two halves kept err by
        # about a fifteenth of their difference to the whole step
        err = norm(halves - whole) / norm(halves) / 15
        if isnan(err)
            error("the evolver of an ExactState gave NaN within a step")
        end
        # its share of the error of the whole evolution, kept ten times above the noise of the
        # exponentials, `rounding_tol` each
        allowed = max(1000 * rounding_tol * abs(h / duration), 10 * rounding_tol)
        if err ≤ allowed
            x, s = halves, s + f
        end
        f *= clamp(0.9 * (allowed / max(err, eps()))^(1 / 5), 0.25, 2.0)
        if f < 1e-9
            error("an evolver depending on time varies too fast for the substeps of an " *
                  "ExactState: take smaller steps")
        end
    end
    return x
end

# the method of an exact state is its own, those of the MPS without use
function evolve_step(::Union{Tdvp, ApproxW}, pre::ExactEvolver, state::ExactState{R},
                     time::Number, dt::Number, step::Int, limits::Limits) where {R <: PM}
    system = state.system
    act(cs) = x -> sum(c * dense_acted(R, system, ts, x) for (c, ts) in zip(cs, pre.terms))
    # a step taken alone, outside evolve, is the whole evolution
    duration = something(pre.duration, dt)
    x = isnothing(pre.functions) ? exact_exp(act(ones(length(pre.terms))), dt, state.state) :
                                   magnus_evolved(act, pre.functions, time, dt, duration,
                                                  state.state)
    return ExactState{R}(system, x, Val(:raw))
end

"""
the linear map of the operator, the channel or the generator of terms `ts` on the arrays of the
exact states of representation `R`, see `dense_acted`
"""
acting(::Type{R}, system::System, ts::Vector{Term}) where R = x -> dense_acted(R, system, ts, x)

"""
a bound of the operator norm of the operator of terms `ts`: the sum over its terms of `|coef|`
times the product of the norms of their blocks, a Jordan-Wigner string having norm one
"""
terms_bound(ts::Vector{Term}) =
    sum(abs(t.coef) * prod((opnorm(b.m) for b in t.blocks); init = 1.0) for t in ts; init = 0.0)

"""
the eigenvalue of the linear map `A` of the smallest real part, `which = :SR`, or of the largest,
`:LR`, and its eigenvector, searched by Krylov from `x` to `rounding_tol`, refused when it does
not converge; `what` names the caller
"""
function exact_eigen(A, x::AbstractArray, which::Symbol, hermitian::Bool, what::String)
    vals, vecs, info = KrylovKit.eigsolve(A, x / norm(x), 1, which; ishermitian = hermitian,
                                          tol = rounding_tol)
    if info.converged == 0
        error("$what of an ExactState did not converge in Krylov to $rounding_tol")
    end
    return first(vals), first(vecs)
end

"""
the operator `a` of a search of an exact state, refused when it is the PreMPO or the MPO of one,
`what` naming the search
"""
function searched_operator(system::System, a, what::String)
    if !(a isa IndexedOp)
        error("$what of an ExactState is given the operator rather than its PreMPO or MPO")
    end
    check_indices(system, a)
    return a
end

# the method and the limits of the MPS without use, `algo` not needed; the states of
# `orthogonal_to` projected out exactly, a weight without use
function ground_state(hamiltonian, state::ExactState{Pure};
                      algo::Union{Nothing, SearchAlgorithm} = nothing, limits::Limits = Limits(),
                      mpo_limits::Limits = Limits(), orthogonal_to::AbstractVector = [],
                      weight::Union{Nothing, Real} = nothing)
    system = state.system
    ts = lower(system, searched_operator(system, hamiltonian, "ground_state"))
    H = acting(Pure, system, ts)
    ms = map(orthogonal_to) do m
        e = m isa ExactState ? m : ExactState(m)
        same_sites(state, e)
        return e.state
    end
    M = isempty(ms) ? zeros(ComplexF64, length(state.state), 0) : reduce(hcat, ms)
    project(v) = isempty(ms) ? v : v - M * (M \ v)
    # the states given lifted above the spectrum, so that the lowest is orthogonal to them
    c = terms_bound(ts) + 1
    A = isempty(ms) ? H : v -> (w = project(v); project(H(w)) + c * (v - w))
    x = project(state.state)
    if iszero(norm(x))
        error("ground_state is started from a state the states of orthogonal_to span")
    end
    _, v = exact_eigen(A, x, :SR, true, "ground_state")
    v = project(v)
    return (energy = real(dot(v, H(v))) / real(dot(v, v)),
            state = normalize(ExactState{Pure}(system, v, Val(:raw))))
end

ground_state(::Any, ::ExactState{Mixed}; kwargs...) =
    error("ground_state searches a pure state: search it pure and mix it")

# the steady state as the eigenvector of the Lindbladian of the largest real part, zero, and
# `residual` the eigenvalue of L†L that the search of a State reaches
function steady_state(lindbladian, state::ExactState{Mixed};
                      algo::Union{Nothing, SearchAlgorithm} = nothing, limits::Limits = Limits(),
                      mpo_limits::Limits = Limits())
    system = state.system
    a = searched_operator(system, lindbladian, "steady_state")
    L = acting(Mixed, system, lower(system, superoperator(a)))
    _, x = exact_eigen(L, state.state, :LR, false, "steady_state")
    t = tr(x)
    if abs(t) ≤ rounding_tol * norm(x)
        error("the eigenvector found by steady_state has no trace: the steady state is not " *
              "unique, as in each sector of a quantity the Lindbladian conserves")
    end
    ρ = ExactState{Mixed}(system, x / t, Val(:raw))
    δ = 1 - hermiticity(ρ)
    if δ > 1e-6
        @warn "the steady state found has a HermiticityDeviation of $(round(δ; sigdigits = 2)): " *
              "the steady state is not unique, see steady_state"
    end
    return (residual = (norm(L(ρ.state)) / norm(ρ))^2, state = ρ)
end

steady_state(::Any, ::ExactState{Pure}; kwargs...) =
    error("steady_state searches a density matrix: mix the state")

# as thermalize of a State, each step an exponential by Krylov, the method and the limits of the
# MPS without use
function thermalize(hamiltonian, beta::Real, state::ExactState{Mixed};
                    algo::EvolutionAlgorithm = Tdvp(), nsteps::Union{Nothing, Int} = nothing,
                    beta_step::Union{Nothing, Real} = nothing, limits::Limits = Limits(),
                    mpo_limits::Limits = Limits(), step::Int = 1)
    nsteps = evolution_steps("thermalize", beta, nsteps, beta_step, "beta_step")
    if nsteps == 0
        return (log_trace = 0., state = step == 1 ? normalize(state) : state)
    end
    system = state.system
    a = searched_operator(system, hamiltonian, "thermalize")
    # ρ ↦ Hρ + ρH, whose exponential at -β/2 is e^{-βH/2} ρ e^{-βH/2}
    K = acting(Mixed, system, lower(system, superoperator(a)))
    x = (step == 1 ? normalize(state) : state).state
    log_trace = 0.
    for k in 0:nsteps - 1
        x = exact_exp(K, -beta / nsteps / 2, x)
        t = real(tr(x))
        if !(t > 0)
            error("the trace of the state reached $t at step $(step + k) of thermalize")
        end
        x /= t
        log_trace += log(t)
    end
    return (; log_trace, state = ExactState{Mixed}(system, x, Val(:raw)))
end

thermalize(::Any, ::Real, ::ExactState{Pure}; kwargs...) =
    error("a thermal state needs a mixed representation: mix the state first")

# the limits of the MPS have no use on an exact state, see the docstring of ExactState
function project_on(state::ExactState{R}, p::IndexedOp; limits::Limits) where R
    system = weaken(state.system, ())
    x = R === Mixed ? on_both_sides(system, p, state.state) :
                      dense_acted(Pure, system, lower(system, p), state.state)
    return normalize(ExactState{R}(state.system, x, Val(:raw)))
end

"""
the vector `x`, of `extra` columns after the product states of `system`, projected on the total
`T` of the sum that `counting` describes: each site of the sum applies the projectors of its
values to each partial total, those that can still reach `T` kept, see `count_bonds`
"""
function count_projected(system::System, sites, locals, T::Int, x::AbstractVector, extra::Int)
    dims = [ dim(s) for s in system.sites ]
    bonds = count_bonds(length(system), sites, locals, T)
    parts = Dict(0 => x)
    for (k, (values, ps)) in zip(sites, locals)
        a, c = prod(dims[k+1:end]; init = 1), prod(dims[1:k-1]; init = 1) * extra
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

function project_count(state::ExactState{R}, sites, locals, T; limits::Limits) where R
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
basis_probabilities(state::ExactState{R}) where R =
    (R === Mixed ? real(diag(state.state)) : abs2.(state.state)) / nonzero_trace(state, "sample")

function sample(state::ExactState; rng = Random.default_rng())
    p = basis_probabilities(state)
    dims = [ dim(s) for s in state.system.sites ]
    c = CartesianIndices(Tuple(reverse(dims)))[draw(x -> p[x + 1], length(p), rand(rng)) + 1]
    return [ c[length(dims) + 1 - i] - 1 for i in eachindex(dims) ]
end

function sample(state::ExactState, pos::Int; rng = Random.default_rng())
    check_positions(state, [pos], "sample")
    dims = [ dim(s) for s in state.system.sites ]
    p = basis_probabilities(state)
    a, c = prod(dims[pos+1:end]; init = 1), prod(dims[1:pos-1]; init = 1)
    q = vec(sum(reshape(p, a, dims[pos], c); dims = (1, 3)))
    return draw(x -> q[x + 1], dims[pos], rand(rng))
end

function write_state(g, state::ExactState{R}) where R
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
