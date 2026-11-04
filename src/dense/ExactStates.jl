# The exact representation of a state of a small system: its whole vector, or its whole density
# matrix, on the basis of the product states ordered as `kron` orders them, the first site
# varying the slowest, the operators acting through their terms, see `dense_acted`.

"""
    ExactState{R}(system, states)
    ExactState(state::State)
    ExactState(system, x)

the state of a small `system` held whole, without truncation: its vector for `R = Pure`, its
density matrix for `R = Mixed`. It is made from a description as `State{R}` takes it, from a
`State`, or from the whole vector or density matrix `x`, normalised, of definite charges on a
system conserving something.

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
        limit = R === Mixed ? dense_max : dense_max^2
        if d > limit
            error("an ExactState{$R} of $(length(system)) sites has the dimension " *
                  "$(Printf.@sprintf("%.0f", d)), beyond the $limit TMS builds")
        end
        shape = R === Mixed ? (Int(d), Int(d)) : (Int(d),)
        if size(x) ≠ shape
            what = R === Mixed ? "a density matrix of $(Int(d)) × $(Int(d))" :
                                 "a vector of $(Int(d)) elements"
            error("an ExactState{$R} of the system holds $what, not an array of size $(size(x))")
        end
        return new{R}(system, Array{ComplexF64}(x))
    end
end

ExactState{R}(args...) where {R <: PM} = ExactState(State{R}(args...))

ExactState(state::State{Pure}) = ExactState{Pure}(state.system, Vector(state), Val(:raw))

ExactState(state::State{Mixed}) = ExactState{Mixed}(state.system, full_matrix(state), Val(:raw))

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

mix(state::ExactState{Mixed}) = state
mix(state::ExactState{Pure}) = ExactState{Mixed}(state.system, state.state * state.state', Val(:raw))

dag(::ExactState{Pure}) = error("dag is meaningless on pure representations")
dag(state::ExactState{Mixed}) = ExactState{Mixed}(state.system, state.state', Val(:raw))

purity(::ExactState{Pure}) = 1.
purity(state::ExactState{Mixed}) = (norm(state) / nonzero_trace(state, "purity"))^2

hermiticity(::ExactState{Pure}) = 1.
hermiticity(state::ExactState{Mixed}) =
    0.5 + 0.5 * real(dot(state.state', state.state)) / nonzero_norm(state, "hermiticity")^2

hermitianize(state::ExactState{Pure}) = state
hermitianize(state::ExactState{Mixed}) =
    ExactState{Mixed}(state.system, (state.state + state.state') / 2, Val(:raw))

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

fidelity(::ExactState{Mixed}, ::ExactState{Mixed}) =
    error("no fidelity between two mixed representations: use hs_fidelity")

hs_fidelity(a::ExactState{Mixed}, b::ExactState{Mixed}) =
    real(inner(a, b)) / (nonzero_norm(a, "hs_fidelity") * nonzero_norm(b, "hs_fidelity"))

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
    odd = map(1:n) do i
        if !fermionic[i]
            return Bool[]
        end
        f = matrix(F, system[i])
        if !isdiag(f)
            error("the parity F of the site $i is not diagonal on its basis, which " *
                  "reduced_density_matrix needs")
        end
        return real(diag(f)) .< 0
    end
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

function reduced_density_matrix(state::ExactState{R}, positions::AbstractVector{Int}) where R
    check_positions(state, positions, "reduced_density_matrix")
    check_increasing(positions, "reduced_density_matrix")
    if isempty(positions)
        return ones(1, 1)
    end
    kept = collect(positions)
    check_reduced(state, kept)
    system = state.system
    n = length(system)
    dims = [ dim(s) for s in system.sites ]
    s = trace_sign_vector(system, kept)
    x = isnothing(s) ? state.state :
        R === Mixed ? s .* state.state .* transpose(s) : s .* state.state
    # the axes of the sites, the last first as kron orders them, those kept before the others
    rest = setdiff(1:n, kept)
    perm = [ n + 1 - i for i in [reverse(kept); reverse(rest)] ]
    a, b = prod(dims[kept]), prod(dims[rest]; init = 1)
    if R === Pure
        m = reshape(permutedims(reshape(x, reverse(dims)...), perm), a, b)
        return unit_trace(m * m')
    end
    t = reshape(permutedims(reshape(x, reverse(dims)..., reverse(dims)...), [perm; perm .+ n]),
                a, b, a, b)
    return unit_trace(sum(t[:, k, :, k] for k in 1:b))
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
