# Evolutions, ground states and gates of the Gaussian backend: a Gaussian unitary turns the covariance
# by its rotation R, M → R M Rᵀ, read exactly from its quadratic generator.

import .GaussianKernels: propagator, quadratic_ground_state, antisymmetric, evolve_branches!, project,
                         project!, occupation_probability, nonunitary_step, linear_applied, linear_annihilates

const nonunitary_sum_what = "the Gaussian backend evolves a superposition of Gaussian states under unitary operators " *
                            "only: under the jumps of a trajectory that are not unitary (projectors, linear operators), " *
                            "start from a state of one component"

# whether the quadratic form k couples the ancilla to the system: a linear term, which changes the
# physical parity
couples_ancilla(k::AbstractMatrix) = any(!iszero, k[1:2, 3:end])

# k real and c imaginary of L = c + (1/4) Σ k_ab γ_a γ_b, as for L = -iH with H Hermitian, the
# imaginary part of k and the real part of c below rounding_tol times the size of L taken as the
# rounding of the expansion of the operator (5.6e-17 for an exactly Hermitian matrix with complex
# entries), refused above with its size; with `lossy`, for the ket of a trajectory, c and k as they are,
# k real when its imaginary part is that rounding, so that a unitary step keeps the norm exactly
function unitary_form(c::Number, k::AbstractMatrix, what::String; lossy::Bool = false)
    scale = rounding_tol * (norm(k) + abs(c))
    if lossy
        return complex(c), maximum(abs, imag(k); init = 0.0) > scale ? Matrix{ComplexF64}(k) : real(k)
    end
    defect = max(maximum(abs, imag(k); init = 0.0), abs(real(c)))
    if defect > scale
        throw(ArgumentError("$what: this operator departs from it by a part of size $(round(defect; sigdigits = 2))"))
    end
    im * imag(c), real(k)
end

# exp(k t), the rotation of the Gaussian operator of generator k over t, complex orthogonal for k complex
step_operator(k::AbstractMatrix, t::Real) = eltype(k) <: Real ? propagator(k, t) : exp(k * t)

# the state turned by the Gaussian operator of rotation R, which keeps the parities of the components when
# `keeps` (it does not couple the ancilla): the covariance of one component, the rotation accumulated by
# a sum, its norm unchanged for R real; for R complex, a non-unitary step of one component, its trace
# multiplied by the factor on its norm (`nonunitary_step`)
function turned(state::GaussianState{<:Ket}, R::AbstractMatrix, keeps::Bool)
    parities = keeps ? state.parities : nothing
    if !(eltype(R) <: Real)
        if !isnothing(state.sum)
            throw(ArgumentError(nonunitary_sum_what))
        end
        M, lognorm = nonunitary_step(state.cov, R)
        return single(state, M; parities, trace = state.trace * exp(lognorm))
    end
    if isnothing(state.sum)
        return single(state, antisymmetric(R * state.cov * transpose(R)); parities)
    end
    s = as_sum(state)
    evolve!(s, R)
    superposed(state, s; trace = state.trace, parities)
end

# the physical parity of the state when known and the same for all its components, nothing otherwise: that
# of the ancilla on the even states, ⟨-iγ₁γ₂⟩ = M₁₂ = ±1
function state_parity(state::GaussianState{<:Ket})
    p = state.parities
    isnothing(p) || !allequal(p) ? nothing : first(p)
end

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
    rotations::Dict{Any, Any}
end

evolver_what(pre::GaussianEvolver) =
    pre.lossy ? "the Gaussian backend evolves a trajectory under a quadratic operator" :
                "the Gaussian backend evolves a pure state under -iH with H Hermitian, for now"

# the evolver compiled for the parity σ of the states (nothing: none needed): for a constant one, c (its
# real part the loss of a trajectory, its imaginary part the phase between sectors), h real and whether it
# keeps the parity; for one depending on time, the (c_j, k_j) and whether it keeps the parity
function compiled(pre::GaussianEvolver, σ)
    get!(pre.forms, σ) do
        if isnothing(pre.fs)
            c, k, _ = gaussian_generator(pre.system, reduce(+, pre.ops); axes = pre.axes, parity = σ)
            c0, h = unitary_form(c, k, evolver_what(pre); pre.lossy)
            return c0, h, !couples_ancilla(k)
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
    checked_jumps(state, trajectory_evolver(state, terms, mpo_limits; duration))

# the evolver `e` of the trajectory `state`, refused when one of its jumps does not keep the state Gaussian
function checked_jumps(state::GaussianState, e)
    if e isa JumpEvolver
        for d in e.decays, j in d.jumps
            if !gaussian_jump(state, j)
                throw(ArgumentError("the Gaussian backend draws the jumps that keep a state Gaussian: unitary gates of " *
                                    "quadratic generators, products of projectors on occupations, linear operators " *
                                    "(c, c†, σ⁻, …); $(obs_name(j)) is none of them"))
            end
        end
    end
    e
end

# whether the jump j keeps a Gaussian state Gaussian, applied exactly by `gaussian_applied`
function gaussian_jump(state::GaussianState, j::IndexedOp{IsOperator})
    s = scalararg(j)
    if !isnothing(quadratic_gate(state, s))
        return true
    end
    p = polynomial(state.system, state.layout, lower(weaken(state.system, ()), s))
    m = maximum((isempty(α) ? 0 : last(α) for α in keys(p.terms)); init = 0)
    length(p) == 1 || !isnothing(occupation_projector(p)) || !isnothing(monomial_linear(p, m))
end

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
                                   Dict{Any, Any}())
    quadratic = all(op -> degree(polynomial(system, state.layout, lower(system, op))) <= 2, ops)
    # the sectors of the parity when the whole evolver conserves it and needs them
    whole = polynomial(system, state.layout, lower(system, reduce(+, ops)))
    p = pre(!quadratic && conserves_parity(whole))
    # compiled now for the parities the state needs, so that an error comes at once
    σs = step_parities(p, state)
    if length(σs) == 2 && !isnothing(fs)
        throw(ArgumentError(branches_what))
    end
    for σ in σs
        compiled(p, σ)
    end
    p
end

# the parities the evolver is compiled for to step the state: none needed, that of the state, or both
# sectors when the components differ or are not known
function step_parities(pre::GaussianEvolver, state::GaussianState{<:Ket})
    if !pre.sectors
        return [nothing]
    end
    σ = state_parity(state)
    isnothing(σ) ? [1, -1] : [σ]
end

const branches_what = "the Gaussian backend evolves a state lying in both sectors of the parity (a ring) under a " *
                      "constant operator only, for now; for values even in the parity, project it on each sector, " *
                      "(Id + P)/2 and (Id - P)/2 with P the parity, and evolve them apart"

function evolve_step(::Union{Tdvp, ApproxW}, pre::GaussianEvolver, state::GaussianState{<:Ket},
                     time::Number, dt::Number, step::Int, limits::Limits)
    if iszero(dt)
        return state
    end
    if !(dt isa Real)
        throw(ArgumentError("the Gaussian backend evolves by real time steps, not $dt"))
    end
    σs = step_parities(pre, state)
    if length(σs) == 2
        if !isnothing(pre.fs)
            throw(ArgumentError(branches_what))
        end
        return branches_step(pre, state, dt)
    end
    σ = only(σs)
    if isnothing(pre.fs)
        c, h, even = compiled(pre, σ)
        return lost(turned(state, rotation(pre, σ, h, dt), even), real(c) * dt)
    end
    forms, even = compiled(pre, σ)
    R, loss = magnus_rotation(pre, forms, time, dt)
    lost(turned(state, R, even), loss)
end

# the rotation of a step of the constant evolver in the sector σ, kept for the last step only: the steps
# of a trajectory between its jumps are of any length
function rotation(pre::GaussianEvolver, σ, h::AbstractMatrix, dt::Real)
    key = (σ, Float64(dt))
    if !haskey(pre.rotations, key)
        filter!(kv -> first(first(kv)) !== σ, pre.rotations)
        pre.rotations[key] = step_operator(h, dt)
    end
    pre.rotations[key]
end

# the state as a sum of components of known parities: its components cut, where not known, into their
# projections on the two sectors (the occupations of the ancilla, exactly), those of weight √p |a| below
# 1e-13 times the largest of both sectors left out, as `project!` does within one
function in_sectors(state::GaussianState{<:Ket})
    if !isnothing(state.parities)
        return as_sum(state), state.parities
    end
    s = as_sum(state)
    fix!(s)
    parts = [(g, a, σ, occupation_probability(g.cov, 1, Float64(σ))) for (g, a) in zip(s.comps, s.coefs) for σ in (1, -1)]
    logw = [p > 0 ? log(p) / 2 + a.logabs : -Inf for (_, a, _, p) in parts]
    largest = maximum(logw)
    comps, coefs, parities = PureGaussian[], LogComplex[], Int[]
    for ((g, a, σ, _), w) in zip(parts, logw)
        if w > -Inf && w >= largest + log(1e-13)
            p, h = project(g, 1, σ == -1)
            push!(comps, h)
            push!(coefs, a * LogComplex(log(p) / 2, 1))
            push!(parities, σ)
        end
    end
    K = length(comps)
    if K > state.max_terms
        throw(ArgumentError("the sectors of the parity would give the state $K Gaussian components, above " *
                            "max_terms = $(state.max_terms): raise it with GaussianState(…; max_terms)"))
    end
    GaussianSum(comps, coefs, Random.Xoshiro(reference_seed)), parities
end

# a step of a constant evolver on a state whose components lie in both sectors, each under the Gaussian
# unitary of its sector with their relative phase (`evolve_branches!`)
function branches_step(pre::GaussianEvolver, state::GaussianState{<:Ket}, dt::Real)
    if iszero(state.trace)
        # the null state
        return state
    end
    s, parities = in_sectors(state)
    steps = Dict{Int, Tuple{Matrix{Float64}, Float64, ComplexF64}}()
    for σ in (1, -1)
        c, h, even = compiled(pre, σ)
        if !even
            throw(ErrorException("an evolver compiled in a sector couples the ancilla"))
        end
        if !(eltype(h) <: Real)
            throw(ArgumentError(nonunitary_sum_what))
        end
        steps[σ] = (h, Float64(dt), c)
    end
    evolve_branches!(s, parities, steps)
    # the sectors are orthogonal: a loss common to both multiplies the trace, another needs the norm
    losses = [real(steps[σ][3]) for σ in (1, -1)]
    trace = allequal(losses) ? state.trace * exp(2 * first(losses) * dt) : nothing
    superposed(state, s; trace, parities)
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
        Y[1:m, 1:m] = step_operator(g, h) * X[1:m, 1:m]
        Y[m+1, m+1] += h * real(c0)
        Y
    end
    magnus(t, h, X) = foldl((X, cs) -> stage(X, cs, h), step_coefs(pre.fs, t, h, 4); init = X)
    X = zeros(pre.lossy ? ComplexF64 : Float64, m + 1, m + 1)
    X[1:m, 1:m] = I(m)
    X = substepped(magnus, time, dt, something(pre.duration, dt), X, "a GaussianState")
    R = X[1:m, 1:m]
    (isreal(R) ? real(R) : R), real(X[m+1, m+1])
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
    if isnothing(generator) || !couples_ancilla(generator[2])
        # quadratic in each sector of the physical parity, which it conserves (a ring, or a term in the parity,
        # which lands on the ancilla), or the error
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
    found = GaussianState{Pure}(system, state.axes, state.layout, g.cov, nothing, nothing, 1.0, state.max_terms)
    (energy = real(c) + g.energy, state = found)
end

ground_state(::Any, ::GaussianState{Trajectory}; kwargs...) = not_searched()

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
    found = GaussianState{Pure}(system, state.axes, state.layout, cov, nothing, [σ], 1.0, state.max_terms)
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

function apply(a::IndexedOp{IsChannel}, state::GaussianState{Pure}; limits::Limits = Limits())
    if !is_pure(a)
        throw(ArgumentError("$(obs_name(a)) has more than one Kraus operator, which a pure state cannot take: " *
                            "draw it on a trajectory, CreateState{Trajectory}"))
    end
    gaussian_applied(kraus(a), state)
end

apply(a::IndexedOp{IsGenerator}, ::GaussianState; limits::Limits = Limits(), rng::AbstractRNG = Random.default_rng()) =
    throw(ArgumentError("apply takes an operator or a channel, and $(obs_name(a)) is a generator"))

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
    generator = quadratic_gate(state, s)
    if !isnothing(generator)
        c0, k, _ = generator
        _, h = unitary_form(c0, k, "the Gaussian backend applies exp(A) with A anti-Hermitian, unitary gates")
        return turned(state, propagator(h, 1), !couples_ancilla(k))
    end
    p = polynomial(state.system, state.layout, lower(weaken(state.system, ()), s))
    if length(p) == 1 && isnothing(state.sum)
        α, c = only(p.terms)
        # c c(α) multiplies the norm by |c|, c(α) being unitary
        return single(state, flipped_covariance(state.cov, α); parities = flipped_parities(state.parities, [α]),
                      trace = abs2(c) * state.trace)
    end
    occupations = occupation_projector(p)
    if !isnothing(occupations) && isnothing(state.sum)
        return projected(state, occupations; null = true)
    end
    linear = isnothing(state.sum) ? monomial_linear(p, size(state.cov, 1)) : nothing
    if !isnothing(linear)
        β, u = linear
        if linear_annihilates(state.cov, u)
            return scaled(state, 0)
        end
        M, w = linear_applied(state.cov, u)
        return single(state, flipped_covariance(M, β); parities = flipped_parities(state.parities, collect(keys(p.terms))),
                      trace = w * state.trace)
    end
    applied_poly(state, p)
end

# the generator of the gate s, exp(A) with A quadratic, as `gaussian_generator` gives it; nothing for another
# operator or a gate whose A is not quadratic
function quadratic_gate(state::GaussianState, s)
    A = gate_generator(s)
    if isnothing(A)
        return nothing
    end
    try
        gaussian_generator(state.system, A; axes = state.axes)
    catch e
        if !(e isa ArgumentError)
            rethrow()
        end
        nothing
    end
end

# the covariance of c(α)|ψ⟩: the rows and columns α change sign
function flipped_covariance(M::AbstractMatrix, α::AbstractVector{Int})
    D = ones(size(M, 1))
    D[α] .= -1
    D .* M .* transpose(D)
end

# the parity of one component after one of the monomials αs, nothing when they disagree or it is not known
function flipped_parities(parities, αs)
    if isnothing(parities)
        return nothing
    end
    after = unique(flipped(only(parities), α) for α in αs)
    length(after) == 1 ? after : nothing
end

# (β, u) when p = c(β) Σ_a u_a γ_a, its monomials one Majorana longer than their common part β: a linear
# operator behind the strings of Jordan-Wigner and the ancilla (c_i, σ⁻), whose c(β ∪ {a}) is ±c(β)γ_a;
# nothing for another polynomial
function monomial_linear(p::MajoranaPoly, m::Int)
    αs = collect(keys(p.terms))
    if isempty(αs)
        return nothing
    end
    β = reduce(intersect, αs)
    if any(α -> length(α) != length(β) + 1, αs)
        return nothing
    end
    u = zeros(ComplexF64, m)
    for α in αs
        a = only(setdiff(α, β))
        u[a] = p.terms[α] * (iseven(count(>(a), β)) ? 1 : -1)
    end
    β, u
end

# the state times `c`
function scaled(state::GaussianState{<:Ket}, c::Number)
    if isnothing(state.sum)
        return single(state, state.cov; trace = abs2(c) * state.trace)
    end
    s = as_sum(state)
    s.coefs .*= Ref(LogComplex(c))
    superposed(state, s; trace = abs2(c) * state.trace, state.parities)
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
    parities = isnothing(state.parities) ? nothing : [flipped(q, α) for q in state.parities for (α, _) in p.terms]
    superposed(state, s; parities)
end

import .GaussianKernels: occupation_probability, project_covariance!, MajoranaPoly

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

# the state projected on these occupations of modes, exactly, not normalised; when that has probability
# zero, an error, or with `null` on one component the null state
function projected(state::GaussianState{<:Ket}, occupations; null::Bool = false)
    if isnothing(state.sum)
        M = copy(state.cov)
        t = state.trace
        for (m, o) in occupations
            s = o ? -1.0 : 1.0
            q = occupation_probability(M, m, s)
            if q == 0 && null
                return scaled(state, 0)
            elseif q == 0
                throw(ArgumentError("the state cannot be projected on a result of probability zero"))
            end
            project_covariance!(M, m, s, q)
            t *= q
        end
        return single(state, M; trace = t, parities = projected_parities(state.parities, occupations))
    end
    s = as_sum(state)
    parities = state.parities
    for (m, o) in occupations
        kept = project!(s, m, o)
        parities = isnothing(parities) ? nothing : parities[kept]
    end
    superposed(state, s; parities = projected_parities(parities, occupations, length(s.comps)))
end

# the parities of the components after a projection on these occupations of modes, which commute with the
# parity: that of the ancilla (mode 1) when projected, empty or occupied
function projected_parities(parities, occupations, K::Int = 1)
    for (m, o) in occupations
        if m == 1
            return fill(o ? -1 : 1, K)
        end
    end
    parities
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
