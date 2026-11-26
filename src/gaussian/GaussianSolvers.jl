# Evolutions, ground states and gates of the Gaussian backend: a Gaussian unitary turns the covariance
# by its rotation R, M → R M Rᵀ, read exactly from its quadratic generator.

import .GaussianKernels: propagator, quadratic_ground_state, antisymmetric

# whether the quadratic form k couples the ancilla to the system: a linear term, which changes the
# physical parity
couples_ancilla(k::AbstractMatrix) = any(!iszero, k[1:2, 3:end])

# k real and c imaginary of L = c + (1/4) Σ k_ab γ_a γ_b, as for L = -iH with H Hermitian, the
# imaginary part of k and the real part of c below rounding_tol times the size of L taken as the
# rounding of the expansion of the operator (5.6e-17 for an exactly Hermitian matrix with complex
# entries), refused above with its size
function unitary_form(c::Number, k::AbstractMatrix, what::String)
    defect = max(maximum(abs, imag(k); init = 0.0), abs(real(c)))
    if defect > rounding_tol * (norm(k) + abs(c))
        throw(ArgumentError("$what: this operator departs from it by a part of size $(round(defect; sigdigits = 2))"))
    end
    im * imag(c), real(k)
end

# the state turned by the Gaussian unitary of rotation R, which keeps the parity when `definite`: the
# covariance of one component, the rotation accumulated by a sum, its norm unchanged
function turned(state::GaussianState{Pure}, R::AbstractMatrix, definite::Bool)
    if isnothing(state.sum)
        return single(state, antisymmetric(R * state.cov * transpose(R)); definite = state.definite && definite)
    end
    s = as_sum(state)
    evolve!(s, R)
    superposed(state, s; trace = state.trace)
end

"""
the evolver of a Gaussian state: L = c + (1/4) Σ k_ab γ_a γ_b, k real, the rotations exp(k dt) already
computed, and whether L keeps the physical parity
"""
struct GaussianEvolver
    k::Matrix{Float64}
    rotations::Dict{Float64, Matrix{Float64}}
    even::Bool
end

prepared_evolver(state::GaussianState{Pure}, a::Union{IndexedOp, Pair}, mpo_limits::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    prepared_evolver(state, [a], mpo_limits; duration)

function prepared_evolver(state::GaussianState{Pure}, terms::AbstractVector, ::Limits = Limits();
                          duration::Union{Nothing, Number} = nothing)
    ops, fs = split_terms(terms)
    if !isnothing(fs)
        throw(ArgumentError("the Gaussian backend evolves under an evolver that does not depend on time, for now"))
    end
    if !all(op -> op isa IndexedOp{IsOperator}, ops)
        throw(ArgumentError("the Gaussian backend evolves a pure state under an operator, -iH"))
    end
    c, k, _ = gaussian_generator(state.system, reduce(+, ops); axes = state.axes)
    _, h = unitary_form(c, k, "the Gaussian backend evolves a pure state under -iH with H Hermitian, for now")
    GaussianEvolver(h, Dict{Float64, Matrix{Float64}}(), !couples_ancilla(k))
end

function evolve_step(::Union{Tdvp, ApproxW}, pre::GaussianEvolver, state::GaussianState{Pure},
                     time::Number, dt::Number, step::Int, limits::Limits)
    if iszero(dt)
        return state
    end
    if !(dt isa Real)
        throw(ArgumentError("the Gaussian backend evolves by real time steps, not $dt"))
    end
    R = get!(() -> propagator(pre.k, dt), pre.rotations, Float64(dt))
    turned(state, R, pre.even)
end

"""
    ground_state(hamiltonian, state::GaussianState; kwargs...)

the ground state of a quadratic Hamiltonian, from its normal modes (`quadratic_ground_state`), with a
line of log when it is degenerate within 1e-10 of the width of its spectrum; `algo`, `limits` and
`mpo_limits` are without use.
"""
function ground_state(hamiltonian, state::GaussianState{Pure};
                      algo::Union{Nothing, SearchAlgorithm} = nothing, limits::Limits = Limits(),
                      mpo_limits::Limits = Limits(), orthogonal_to::AbstractVector = [],
                      weight::Union{Nothing, Real} = nothing)
    system = state.system
    if !isempty(orthogonal_to)
        throw(ArgumentError("the Gaussian backend searches the ground state only, not states orthogonal to others"))
    end
    if is_charged(system)
        throw(ArgumentError("the Gaussian backend searches the ground state of a system without conserved " *
                            "quantities, for now"))
    end
    H = followed_operator(searched_operator(system, hamiltonian, "ground_state"))
    c, k, _ = gaussian_generator(system, H; axes = state.axes)
    # H = c + (1/4) Σ k γγ Hermitian is -iH = -ic + (1/4) Σ -ik γγ anti-Hermitian
    _, h = unitary_form(-im * c, -im * k, "ground_state needs a Hermitian Hamiltonian")
    g = quadratic_ground_state(h, 1)
    ε = g.energies
    # the gap within the even states: two modes more, or the lowest occupied one moved up
    gap = g.occupied ? ε[2] - ε[1] : ε[1] + ε[2]
    width = sum(ε)
    if gap < 1e-10 * width
        @info "ground_state of a GaussianState: the ground level is degenerate within " *
              "$(round(gap; sigdigits = 2)) for a spectrum of width $(round(width; sigdigits = 2)), " *
              "one of its states is returned"
    end
    found = GaussianState{Pure}(system, state.axes, state.layout, g.cov, nothing, !couples_ancilla(k), 1.0,
                                state.max_terms)
    (energy = real(c) + g.energy, state = found)
end

# the generator A of a gate written exp(A), placed on its sites, nothing for another gate
function gate_generator(g::AtIndex)
    op = g.op
    if op isa Operator && op.expr isa ExpOp
        return AtIndex(op.expr.arg, g.index)
    elseif op isa ExpOp
        return AtIndex(op.arg, g.index)
    end
    nothing
end
gate_generator(::IndexedOp) = nothing

"""
    apply(op, state::GaussianState{Pure})

the gates of `op` applied to the state, a product applied factor by factor from the right, a factor
changing the norm only: a gate written exp(A), A quadratic and anti-Hermitian (the rotations of
Circuits, `exp` written by the user), turns every component; any other operator, a polynomial in
Majoranas, has each of its monomials applied exactly to each component, which multiplies their
number, refused beyond `max_terms`.
"""
function apply(op::IndexedOp{IsOperator}, state::GaussianState{Pure}; limits::Limits = Limits())
    c, s = scalarcoef(op), scalararg(op)
    if c != 1
        # a factor in front of the gates changes the norm only, as for an ExactState
        return scaled(apply(s, state; limits), c)
    end
    if s isa ProdOp
        return foldr((g, st) -> apply(g, st; limits), s.subs; init = state)
    end
    A = gate_generator(s)
    if !isnothing(A)
        c0, k, _ = gaussian_generator(state.system, A; axes = state.axes)
        _, h = unitary_form(c0, k, "the Gaussian backend applies exp(A) with A anti-Hermitian, unitary gates")
        return turned(state, propagator(h, 1), !couples_ancilla(k))
    end
    p = polynomial(state.system, state.layout, lower(weaken(state.system, ()), s))
    if length(p) == 1 && isnothing(state.sum)
        α, c = only(p.terms)
        D = ones(size(state.cov, 1))
        D[α] .= -1
        # c c(α) multiplies the norm by |c|, c(α) being unitary
        return single(state, D .* state.cov .* transpose(D); trace = abs2(c) * state.trace)
    end
    applied_poly(state, p)
end

# the state times `c`
function scaled(state::GaussianState{Pure}, c::Number)
    if isnothing(state.sum)
        return single(state, state.cov; trace = abs2(c) * state.trace)
    end
    s = as_sum(state)
    s.coefs .*= Ref(LogComplex(c))
    superposed(state, s; trace = abs2(c) * state.trace)
end

# the polynomial p applied to the state: each monomial to each component, exactly, with its phase
function applied_poly(state::GaussianState{Pure}, p::MajoranaPoly)
    K = ncomponents(state) * length(p)
    if K > state.max_terms
        throw(ArgumentError("the operator would give the state $K Gaussian components, above max_terms = " *
                            "$(state.max_terms): raise it with GaussianState(…; max_terms)"))
    end
    s = as_sum(state)
    fix!(s)
    comps, coefs = PureGaussian[], LogComplex[]
    for (g, a) in zip(s.comps, s.coefs), (α, c) in p.terms
        push!(comps, apply_monomial(α, g))
        push!(coefs, a * c)
    end
    s.comps, s.coefs = comps, coefs
    relink!(s)
    superposed(state, s)
end
