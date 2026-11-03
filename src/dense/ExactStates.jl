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
