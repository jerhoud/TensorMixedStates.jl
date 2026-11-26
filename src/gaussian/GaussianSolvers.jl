# Evolutions, ground states and gates of the Gaussian backend: a Gaussian unitary turns the covariance
# by its rotation R, M → R M Rᵀ, read exactly from its quadratic generator.

import .GaussianKernels: propagator, quadratic_ground_state, antisymmetric

# whether the quadratic form k couples the ancilla to the system: a linear term, which changes the
# physical parity
couples_ancilla(k::AbstractMatrix) = any(!iszero, k[1:2, 3:end])

# k real and c imaginary of L = c + (1/4) Σ k_ab γ_a γ_b, as for L = -iH with H Hermitian, the
# imaginary part of k and the real part of c below rounding_tol times the size of L taken as the
# rounding of the expansion of the operator (5.6e-17 for an exactly Hermitian matrix with complex
# entries), refused above with its size; with `lossy`, the real part of c kept, as for the
# -iH - Γ/2 of a trajectory whose Σ L†L is the constant Γ
function unitary_form(c::Number, k::AbstractMatrix, what::String; lossy::Bool = false)
    defect = max(maximum(abs, imag(k); init = 0.0), lossy ? 0.0 : abs(real(c)))
    if defect > rounding_tol * (norm(k) + abs(c))
        throw(ArgumentError("$what: this operator departs from it by a part of size $(round(defect; sigdigits = 2))"))
    end
    (lossy ? complex(c) : im * imag(c)), real(k)
end

# the state turned by the Gaussian unitary of rotation R, which keeps the parity when `definite`: the
# covariance of one component, the rotation accumulated by a sum, its norm unchanged
function turned(state::GaussianState{<:Ket}, R::AbstractMatrix, definite::Bool)
    if isnothing(state.sum)
        return single(state, antisymmetric(R * state.cov * transpose(R)); definite = state.definite && definite)
    end
    s = as_sum(state)
    evolve!(s, R)
    superposed(state, s; trace = state.trace)
end

# the physical parity of a state of one component and definite parity, nothing for another: that of the
# ancilla on the even states, ⟨-iγ₁γ₂⟩ = M₁₂ = ±1
state_parity(state::GaussianState{<:Ket}) =
    isnothing(state.sum) && state.definite ? Int(sign(state.cov[1, 2])) : nothing

"""
the evolver of a Gaussian state under L(t) = Σ_j f_j(t) L_j, nothing for the functions of a constant
evolver: its operators compiled for each physical parity of the states it steps (nothing when they do
not need it), and the rotations of the constant one by step; `lossy` for the ket of a trajectory, whose
L may have a constant real part, which changes its norm
"""
struct GaussianEvolver
    system::System
    axes::Dict{Int, Symbol}
    ops::Vector
    fs::Union{Nothing, Vector}
    duration::Union{Nothing, Number}
    lossy::Bool
    # whether the evolver is quadratic only in the sector of the parity of the state (a ring of qubits)
    sectors::Bool
    forms::Dict{Any, Any}
    rotations::Dict{Any, Matrix{Float64}}
end

evolver_what(pre::GaussianEvolver) =
    pre.lossy ? "the Gaussian backend evolves a trajectory under -iH - Γ/2 with H Hermitian and Γ the sum of " *
                "the L†L of its jumps a constant (jumps by Paulis or unitaries), for now" :
                "the Gaussian backend evolves a pure state under -iH with H Hermitian, for now"

# the evolver compiled for the parity σ of the states (nothing: none needed): for a constant one, the
# real part of c, h real and whether it keeps the parity; for one depending on time, the (c_j, k_j) and
# whether it keeps the parity
function compiled(pre::GaussianEvolver, σ)
    get!(pre.forms, σ) do
        if isnothing(pre.fs)
            c, k, _ = gaussian_generator(pre.system, reduce(+, pre.ops); axes = pre.axes, parity = σ)
            c0, h = unitary_form(c, k, evolver_what(pre); pre.lossy)
            return real(c0), h, !couples_ancilla(k)
        end
        forms = [gaussian_generator(pre.system, op; axes = pre.axes, parity = σ) for op in pre.ops]
        [(c, k) for (c, k, _) in forms], !any(couples_ancilla(k) for (_, k, _) in forms)
    end
end

prepared_evolver(state::GaussianState, a::Union{IndexedOp, Pair}, mpo_limits::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    prepared_evolver(state, [a], mpo_limits; duration)

prepared_evolver(state::GaussianState{Pure}, terms::AbstractVector, ::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    gaussian_evolver(state, terms, duration; lossy = false)

# a trajectory draws its jumps, see `trajectory_evolver`, its ket following the operators of
# `ket_evolver`
prepared_evolver(state::GaussianState{Trajectory}, terms::AbstractVector, mpo_limits::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    trajectory_evolver(state, terms, mpo_limits; duration)

ket_evolver(state::GaussianState{Trajectory}, terms::AbstractVector, ::Limits;
            duration::Union{Nothing, Number} = nothing) =
    gaussian_evolver(state, terms, duration; lossy = true)

# the evolver of the state under the terms, see `GaussianEvolver`
function gaussian_evolver(state::GaussianState, terms::AbstractVector, duration::Union{Nothing, Number}; lossy::Bool)
    ops, fs = split_terms(terms)
    if !all(op -> op isa IndexedOp{IsOperator}, ops)
        throw(ArgumentError("the Gaussian backend evolves a pure state under an operator, -iH"))
    end
    system = state.system
    pre(sectors) = GaussianEvolver(system, state.axes, ops, fs, duration, lossy, sectors, Dict{Any, Any}(),
                                   Dict{Any, Matrix{Float64}}())
    quadratic = all(op -> degree(polynomial(system, state.layout, lower(system, op))) <= 2, ops)
    # the sectors of the parity when the whole evolver conserves it and needs them
    whole = polynomial(system, state.layout, lower(system, reduce(+, ops)))
    p = pre(!quadratic && conserves_parity(whole))
    # compiled now for the state, so that an error comes at once
    compiled(p, step_parity(p, state))
    p
end

# the parity the evolver is compiled for to step the state: that of the state when it needs the sectors,
# the error of the operator that is not quadratic when the state has none
function step_parity(pre::GaussianEvolver, state::GaussianState{<:Ket})
    if !pre.sectors
        return nothing
    end
    σ = state_parity(state)
    if isnothing(σ)
        compiled(pre, nothing)
    end
    σ
end

function evolve_step(::Union{Tdvp, ApproxW}, pre::GaussianEvolver, state::GaussianState{<:Ket},
                     time::Number, dt::Number, step::Int, limits::Limits)
    if iszero(dt)
        return state
    end
    if !(dt isa Real)
        throw(ArgumentError("the Gaussian backend evolves by real time steps, not $dt"))
    end
    σ = step_parity(pre, state)
    if isnothing(pre.fs)
        loss, h, even = compiled(pre, σ)
        R = get!(() -> propagator(h, dt), pre.rotations, (σ, Float64(dt)))
        return lost(turned(state, R, even), loss * dt)
    end
    forms, even = compiled(pre, σ)
    R, loss = magnus_rotation(pre, forms, time, dt)
    lost(turned(state, R, even), loss)
end

# the state with its norm multiplied by e^l
lost(state::GaussianState, l::Real) = iszero(l) ? state : scaled(state, exp(l))

# the rotation R of the step from `time` to `time + dt`, dR/dt = (Σ_j f_j(t) k_j) R, and the integral of
# the real part of Σ_j f_j(t) c_j, the logarithm of the factor of the norm, by the commutator-free
# Magnus integrator of order 4 on the substeps of `substepped`, as for an ExactState, that logarithm
# carried as the last diagonal element so that the substeps control it too; the sum anti-Hermitian at
# each substep, its terms need not be
function magnus_rotation(pre::GaussianEvolver, forms, time::Number, dt::Number)
    m = size(last(first(forms)), 1)
    generator(cs) = unitary_form(sum(c * c0 for (c, (c0, _)) in zip(cs, forms)),
                                 sum(c * k for (c, (_, k)) in zip(cs, forms)), evolver_what(pre); pre.lossy)
    function stage(X, cs, h)
        c0, g = generator(cs)
        Y = copy(X)
        Y[1:m, 1:m] = propagator(g, h) * X[1:m, 1:m]
        Y[m+1, m+1] += h * real(c0)
        Y
    end
    magnus(t, h, X) = foldl((X, cs) -> stage(X, cs, h), step_coefs(pre.fs, t, h, 4); init = X)
    X = zeros(m + 1, m + 1)
    X[1:m, 1:m] = Matrix{Float64}(I, m, m)
    X = substepped(magnus, time, dt, something(pre.duration, dt), X, "a GaussianState")
    X[1:m, 1:m], X[m+1, m+1]
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
    generator = try
        gaussian_generator(system, H; axes = state.axes)
    catch e
        if !(e isa ArgumentError)
            rethrow()
        end
        nothing
    end
    if isnothing(generator)
        # quadratic in each sector of the physical parity, which it conserves (a ring), or the error
        found = [sector_ground_state(system, H, state, σ) for σ in (1, -1)]
        if any(isnothing, found)
            gaussian_generator(system, H; axes = state.axes)
        end
        a, b = found
        if abs(a.energy - b.energy) < 1e-10 * max(a.width, b.width)
            @info "ground_state of a GaussianState: the ground levels of the two parities differ by " *
                  "$(round(abs(a.energy - b.energy); sigdigits = 2)), one of them is returned"
        end
        best = a.energy <= b.energy ? a : b
        degenerate_info(best.levels)
        return (energy = best.energy, state = best.state)
    end
    c, k, _ = generator
    # H = c + (1/4) Σ k γγ Hermitian is -iH = -ic + (1/4) Σ -ik γγ anti-Hermitian
    _, h = unitary_form(-im * c, -im * k, "ground_state needs a Hermitian Hamiltonian")
    g = quadratic_ground_state(h, 1)
    degenerate_info(g)
    found = GaussianState{Pure}(system, state.axes, state.layout, g.cov, nothing, !couples_ancilla(k), 1.0,
                                state.max_terms)
    (energy = real(c) + g.energy, state = found)
end

# the line of log of a ground level degenerate within 1e-10 of the width of the spectrum, the gap taken
# within the states of the same parity: two modes more, or the lowest occupied one moved up
function degenerate_info(g)
    ε = g.energies
    if length(ε) < 2
        return nothing
    end
    gap = g.occupied ? ε[2] - ε[1] : ε[1] + ε[2]
    width = sum(ε)
    if gap < 1e-10 * width
        @info "ground_state of a GaussianState: the ground level is degenerate within " *
              "$(round(gap; sigdigits = 2)) for a spectrum of width $(round(width; sigdigits = 2)), " *
              "one of its states is returned"
    end
    nothing
end

# the ground state of H among the states of physical parity σ, nothing when H is not quadratic there:
# the system's modes in their ground state of parity σ, the ancilla in that of parity σ, all even
function sector_ground_state(system::System, H, state::GaussianState{Pure}, σ::Int)
    generator = try
        gaussian_generator(system, H; axes = state.axes, parity = σ)
    catch e
        if !(e isa ArgumentError)
            rethrow()
        end
        return nothing
    end
    c, k, M = generator
    if couples_ancilla(k)
        return nothing
    end
    _, h = unitary_form(-im * c, -im * k, "ground_state needs a Hermitian Hamiltonian")
    g = quadratic_ground_state(h[3:end, 3:end], σ)
    cov = zeros(M, M)
    cov[1, 2], cov[2, 1] = σ, -σ
    cov[3:end, 3:end] = g.cov
    found = GaussianState{Pure}(system, state.axes, state.layout, cov, nothing, true, 1.0, state.max_terms)
    (energy = real(c) + g.energy, width = sum(g.energies), levels = g, state = found)
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
    apply(op, state::GaussianState)

the gates of `op` applied to the state, a product applied factor by factor from the right, a factor
changing the norm only: a gate written exp(A), A quadratic and anti-Hermitian (the rotations of
Circuits, `exp` written by the user), turns every component; any other operator, a polynomial in
Majoranas, has each of its monomials applied exactly to each component, which multiplies their
number, refused beyond `max_terms`. A trajectory draws the Kraus operators of a channel, and takes
the weight an operator keeps of it, as for a State.
"""
apply(op::IndexedOp{IsOperator}, state::GaussianState{Pure}; limits::Limits = Limits()) = gaussian_applied(op, state)

apply(a::IndexedOp{IsOperator}, state::GaussianState{Trajectory}; limits::Limits = Limits(),
      rng::AbstractRNG = Random.default_rng()) = trajectory_applied(a, state, limits)

apply(a::IndexedOp{IsChannel}, state::GaussianState{Trajectory}; limits::Limits = Limits(),
      rng::AbstractRNG = Random.default_rng()) = trajectory_applied(a, state, limits, rng)

ket_applied(a::IndexedOp{IsOperator}, state::GaussianState{Trajectory}; limits::Limits = Limits()) =
    gaussian_applied(a, state)

# the operator applied to the ket of the state, see `apply`
function gaussian_applied(op::IndexedOp{IsOperator}, state::GaussianState)
    c, s = scalarcoef(op), scalararg(op)
    if c != 1
        # a factor in front of the gates changes the norm only, as for an ExactState
        return scaled(gaussian_applied(s, state), c)
    end
    if s isa ProdOp
        return foldr((g, st) -> gaussian_applied(g, st), s.subs; init = state)
    end
    A = gate_generator(s)
    generator = nothing
    if !isnothing(A)
        generator = try
            gaussian_generator(state.system, A; axes = state.axes)
        catch e
            if !(e isa ArgumentError)
                rethrow()
            end
            # A not quadratic: exp(A) as a polynomial, below
            nothing
        end
    end
    if !isnothing(generator)
        c0, k, _ = generator
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
function scaled(state::GaussianState{<:Ket}, c::Number)
    if isnothing(state.sum)
        return single(state, state.cov; trace = abs2(c) * state.trace)
    end
    s = as_sum(state)
    s.coefs .*= Ref(LogComplex(c))
    superposed(state, s; trace = abs2(c) * state.trace)
end

(c::Number * state::GaussianState) = scaled(state, c)
(state::GaussianState * c::Number) = scaled(state, c)
(state::GaussianState / c::Number) = scaled(state, inv(c))

# the pure state the trajectory holds, normalised
as_pure(state::GaussianState{Trajectory}) = with_representation(Pure, scaled(state, 1 / sqrt(nonzero_trace(state, "as_pure"))))

# the polynomial p applied to the state: each monomial to each component, exactly, with its phase;
# `cause` names what multiplies the components in the refusal beyond max_terms
function applied_poly(state::GaussianState{<:Ket}, p::MajoranaPoly; cause::String = "the operator")
    if iszero(p)
        # the null state
        return scaled(state, 0)
    end
    K = ncomponents(state) * length(p)
    if K > state.max_terms
        throw(ArgumentError("$cause would give the state $K Gaussian components, above max_terms = " *
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

import .GaussianKernels: occupation_probability, project_covariance!, project!, MajoranaPoly

# the modes and occupations of a projector that is a product of projectors on the occupations of modes,
# Π_m (1/2 + σ_m (i/2) γ_{2m-1}γ_{2m}), σ_m = 1 for occupied (-iγγ = 1 - 2n), as m => occupied; nothing
# for another operator. Decided by the exact equality of the polynomials: a projector whose expansion
# carries rounding takes the general path, exact too.
function occupation_projector(p::MajoranaPoly)
    singles = [α for α in keys(p.terms) if length(α) == 2 && isodd(α[1]) && α[2] == α[1] + 1]
    r = length(singles)
    if r == 0 || length(p) != 2^r
        return nothing
    end
    occupied = [(α[1] + 1) ÷ 2 => p.terms[α] == im / 2^r for α in singles]
    if !all(p.terms[α] == im / 2^r || p.terms[α] == -im / 2^r for α in singles)
        return nothing
    end
    product = prod(MajoranaPoly(Int[] => 0.5, [2m - 1, 2m] => (o ? im : -im) / 2) for (m, o) in occupied)
    product == p ? occupied : nothing
end

# the state projected on these occupations of modes, exactly, not normalised; an error when that has
# probability zero
function projected(state::GaussianState{<:Ket}, occupations)
    if isnothing(state.sum)
        M = copy(state.cov)
        t = state.trace
        for (m, o) in occupations
            s = o ? -1.0 : 1.0
            q = occupation_probability(M, m, s)
            if q == 0
                throw(ArgumentError("the state cannot be projected on a result of probability zero"))
            end
            project_covariance!(M, m, s, q)
            t *= q
        end
        # the occupations commute with the parity
        return single(state, M; trace = t)
    end
    s = as_sum(state)
    for (m, o) in occupations
        project!(s, m, o)
    end
    superposed(state, s)
end

function normalized(state::GaussianState{<:Ket})
    nonzero_trace(state, "project_on")
    isnothing(state.sum) ? single(state, state.cov; trace = 1.0) : scaled(state, 1 / sqrt(state.trace))
end

"""
    project_on(state::GaussianState, p; limits)

the state projected by `p` and normalised: exactly and without new components for a product of
projectors on the occupations of modes (N of a fermion, (1 ± σ)/2 along the axis of a qubit), as a
polynomial otherwise, see `apply`.
"""
function project_on(state::GaussianState{<:Ket}, p::IndexedOp; limits::Limits = Limits())
    poly = polynomial(state.system, state.layout, lower(weaken(state.system, ()), p))
    occupations = occupation_projector(poly)
    if !isnothing(occupations)
        return normalized(projected(state, occupations))
    end
    cause = "this measurement, outside the occupations of the modes,"
    normalized(applied_poly(state, poly; cause))
end

project_count(::GaussianState, sites, locals, T; limits::Limits = Limits()) =
    throw(ArgumentError("the Gaussian backend does not project on a total of operators of one site for now; " *
                        "measure the sites one by one"))

# drawn from the probabilities, without projecting the state
sample(state::GaussianState{<:Ket}, pos::Int; rng = Random.default_rng()) =
    round(Int, sample(state, Level(pos); rng))

# the sites drawn one after the other, each collapsing the state
function sample(state::GaussianState{<:Ket}; rng = Random.default_rng())
    out = Int[]
    for pos in eachindex(state.system.sites)
        r = collapse(state, pos; rng)
        push!(out, r.result)
        state = r.state
    end
    out
end
