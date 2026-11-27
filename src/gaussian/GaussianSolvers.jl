# Evolutions, ground states and gates of the Gaussian backend: a Gaussian unitary turns the covariance
# by its rotation R, M → R M Rᵀ, read exactly from its quadratic generator.

import .GaussianKernels: propagator, quadratic_ground_state, antisymmetric, evolve_branches!, project,
                         project!, occupation_probability, nonunitary_step, nonunitary_pieces, linear_applied,
                         lyapunov_propagator, lyapunov_pieces, purified_covariance

const nonunitary_sum_what = "the Gaussian backend evolves a superposition of Gaussian states under unitary operators " *
                            "only: under the jumps of a trajectory that are not unitary (projectors, linear operators), " *
                            "start from a state of one component"

const nonunitary_parity_what = "the Gaussian backend evolves a state lying in both sectors of the parity (a ring) under " *
                               "unitary operators only: under jumps that are not unitary, project it on each sector, " *
                               "(Id + P)/2 and (Id - P)/2 with P the parity, and evolve them apart"

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

# the state turned by the Gaussian unitary of rotation R, which keeps the parities of the components when
# `keeps` (it does not couple the ancilla): the covariance of one component, the rotation accumulated by
# a sum, its norm unchanged
function turned(state::GaussianState{<:Ket}, R::AbstractMatrix{<:Real}, keeps::Bool)
    parities = keeps ? state.parities : nothing
    if isnothing(state.sum)
        return single(state, antisymmetric(R * state.cov * transpose(R)); parities)
    end
    s = as_sum(state)
    evolve!(s, R)
    superposed(state, s; trace = state.trace, parities)
end

# a gate turning a mixed state: refused when it couples the ancilla, a term odd in the Majoranas that the
# parity-commuting states of the backend cannot follow
function turned(state::GaussianState{Mixed}, R::AbstractMatrix{<:Real}, keeps::Bool)
    if !keeps
        throw(ArgumentError(mixed_odd_what))
    end
    single(state, antisymmetric(R * state.cov * transpose(R)))
end

const mixed_odd_what = "the Gaussian backend holds mixed states that commute with the parity, and this " *
                       "operator has a term odd in the Majoranas (C + C† of a fermion, a field across the axis " *
                       "of the first qubit): apply it to a pure state or to trajectories"

# the state of one component under E^m, E = exp(k) complex orthogonal growing vectors by a bounded factor
# (`nonunitary_pieces`), a QR after each (`nonunitary_step`): its trace multiplied by the factor on the
# norm and by e^{2 loss}, its parities kept when `keeps`
function dissipated(state::GaussianState{<:Ket}, E::AbstractMatrix, m::Int, keeps::Bool, loss::Real)
    if !isnothing(state.sum)
        throw(ArgumentError(nonunitary_sum_what))
    end
    M, lognorm = nonunitary_step(state.cov, E, m)
    single(state, M; parities = keeps ? state.parities : nothing, trace = state.trace * exp(lognorm + 2loss))
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
L may not be anti-Hermitian, which changes its norm; `even` false for a mixed state, compiled as `polynomial`
says
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
    even::Bool
    forms::Dict{Any, Any}
    rotations::Dict{Any, Any}
end

evolver_what(::GaussianEvolver) = "the Gaussian backend evolves a pure or a mixed state under -iH with H Hermitian"

# the evolver compiled for the parity σ of the states (nothing: none needed): for a constant one, c (its
# real part the loss of a trajectory, its imaginary part the phase between sectors), h (complex for the ket
# of a trajectory that is not unitary) and whether it keeps the parity; for one depending on time, the
# (c_j, k_j) and whether it keeps the parity
function compiled(pre::GaussianEvolver, σ)
    get!(pre.forms, σ) do
        if isnothing(pre.fs)
            c, k, _ = gaussian_generator(pre.system, reduce(+, pre.ops); axes = pre.axes, parity = σ, even = pre.even)
            c0, h = unitary_form(c, k, evolver_what(pre); pre.lossy)
            return c0, h, !couples_ancilla(k)
        end
        forms = [gaussian_generator(pre.system, op; axes = pre.axes, parity = σ, even = pre.even) for op in pre.ops]
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

"""
the evolver of a mixed Gaussian state under a constant Lindbladian, -iH with H quadratic and dissipators of
linear operators: dM/dt = X M + M Xᵀ + Y (see `lindblad_form`), and for each length of step its propagator
(`lyapunov_propagator`) over a piece and the number of pieces
"""
struct LyapunovEvolver
    X::Matrix{Float64}
    Y::Matrix{Float64}
    steps::Dict{Float64, Tuple{Matrix{Float64}, Matrix{Float64}, Int}}
end

# under a constant Lindbladian by the equation of Lyapunov; under -iH(t), without dissipators, by the rotations
# of the pure states (Magnus)
function prepared_evolver(state::GaussianState{Mixed}, terms::AbstractVector, ::Limits = Limits();
                          duration::Union{Nothing, Number} = nothing)
    ops, fs = split_terms(terms)
    if !isnothing(fs)
        if !all(op -> op isa IndexedOp{IsOperator}, ops)
            throw(ArgumentError("the Gaussian backend evolves a mixed state under dissipators that are constant, " *
                                "without terms depending on time: draw trajectories, as " *
                                "State{Trajectory}(system, \"FullyMixed\")"))
        end
        pre = gaussian_evolver(state, terms, duration; lossy = false)
        if !last(compiled(pre, nothing))
            throw(ArgumentError(mixed_odd_what))
        end
        return pre
    end
    LyapunovEvolver(lindblad_form(state, ops)..., Dict{Float64, Tuple{Matrix{Float64}, Matrix{Float64}, Int}}())
end

# X and Y of dM/dt = X M + M Xᵀ + Y for the Lindbladian of `ops` on a mixed state: k, the generator of -iH,
# turns M (dM/dt = k M + M kᵀ), and the dissipators of rates r and operators L = Σ_a u_a γ_a make
# B = Σ r ū uᵀ, of dM/dt = 4 Im B - 2 (M Re B + Re B M) by ⟨γ_a γ_b⟩ = δ_ab + i M_ab; so X = k - 2 Re B and
# Y = 4 Im B, the equation of Bravyi and König (2012, eq. (32)) for the Lindbladian of TMS, theirs having a
# factor 2 on the dissipators
function lindblad_form(state::GaussianState{Mixed}, ops)
    system = state.system
    n = 2 + state.layout.count
    hs = []
    B = zeros(ComplexF64, n, n)
    for a in ops
        pure, ds = a isa IndexedOp{IsGenerator} ? unravelled(a) : ([a], [])
        append!(hs, pure)
        for (r, x) in ds
            for j in (x isa IndexedOp{IsChannel} ? kraus_operators(system, x) : [x])
                # a Kraus operator proportional to the identity has a null dissipator
                if !(scalararg(j) isa IdentityOp)
                    u = mixed_jump(state, j, x)
                    B .+= r .* conj(u) .* transpose(u)
                end
            end
        end
    end
    k = zeros(n, n)
    if !isempty(hs)
        c, k0, _ = gaussian_generator(system, sum(hs); axes = state.axes, even = false)
        _, k = unitary_form(c, k0, "the Gaussian backend evolves a mixed state under -iH, H Hermitian, and dissipators")
        if couples_ancilla(k)
            throw(ArgumentError(mixed_odd_what))
        end
    end
    k - 2 * real(B), 4 * imag(B)
end

# the coefficients u of the operator j = Σ_a u_a γ_a of a dissipator on a mixed state, read on its polynomial
# Σ_a i u_a γ₁γ_a; refused when j is not linear in the Majoranas of the backend, as σ⁻ behind the parity of the
# modes before it in the layout, its dissipator not Gaussian on a mixed state
function mixed_jump(state::GaussianState{Mixed}, j, x)
    u = zeros(ComplexF64, 2 + state.layout.count)
    p = state_polynomial(state, lower(weaken(state.system, ()), j))
    linear(α) = length(α) == 2 && α[1] == 1 && α[2] >= 3
    if !all(linear, keys(p.terms))
        bad = maximum(abs(c) for (α, c) in p.terms if !linear(α))
        throw(ArgumentError("the Gaussian backend evolves a mixed state under dissipators of operators linear " *
                            "in its Majoranas (C and dag(C) of fermions placed before any qubit, Sm and Sp of " *
                            "the first qubit of a system of qubits along :Z), and $(obs_name(x)) is not one: " *
                            "draw trajectories, as State{Trajectory}(system, \"FullyMixed\")." *
                            rounding_hint(bad, maximum(abs, values(p.terms)))))
    end
    for (α, c) in p.terms
        u[α[2]] = -im * c
    end
    u
end

function evolve_step(::Union{Tdvp, ApproxW}, pre::LyapunovEvolver, state::GaussianState{Mixed},
                     time::Number, dt::Number, step::Int, limits::Limits)
    if iszero(dt)
        return state
    end
    if !(dt isa Real)
        throw(ArgumentError("the Gaussian backend evolves by real time steps, not $dt"))
    end
    E, G, m = get!(pre.steps, Float64(dt)) do
        m = lyapunov_pieces(pre.X, dt)
        (lyapunov_propagator(pre.X, pre.Y, dt / m)..., m)
    end
    M = state.cov
    for _ in 1:m
        M = E * M * transpose(E) + G
    end
    single(state, antisymmetric(M))
end

"""
    steady_state(lindbladian, state::GaussianState{Mixed}; algo, limits, mpo_limits)

`(residual, state)`: the steady state of a constant Lindbladian of quadratic H and linear jumps, of trace
one, its covariance M solving X M + M Xᵀ + Y = 0 (Bravyi and König 2012, eq. (33)), exactly, and the residual
‖X M + M Xᵀ + Y‖² / ‖X‖², the entries of M being at most 1; refused when it is not unique or too weakly
damped to be found, an eigenvalue of X of real part above -1e-10 ‖X‖ (a mode the dissipators do not reach).

- `algo`, `limits`, `mpo_limits`: without use

# Examples

    r = steady_state(-im * H + sum(0.3 * Dissipator(C)(j) for j in 1:L), GaussianState(State{Mixed}(sys, "FullyMixed")))
"""
function steady_state(lindbladian, state::GaussianState{Mixed}; algo::Union{Nothing, SearchAlgorithm} = nothing,
                      limits::Limits = Limits(), mpo_limits::Limits = Limits())
    ops, fs = split_terms(lindbladian isa AbstractVector ? lindbladian : [lindbladian])
    if !isnothing(fs)
        throw(ArgumentError("steady_state takes a constant Lindbladian, not one depending on time"))
    end
    X, Y = lindblad_form(state, ops)
    # the modes of the system, the ancilla held empty
    s = 3:size(X, 1)
    Xs, Ys = X[s, s], Y[s, s]
    # >= so that a Lindbladian of zero (X = 0) is refused too
    top = maximum(real, eigvals(Xs); init = -Inf)
    if top >= -1e-10 * opnorm(Xs)
        throw(ArgumentError("the Lindbladian has no unique steady state, or one too weakly damped to be found here " *
                            "(an eigenvalue of real part $(round(top; sigdigits = 2)) of the Lyapunov equation): " *
                            "evolve the state towards it instead"))
    end
    M = copy(state.cov)
    M[s, s] = antisymmetric(lyap(Xs, Ys))
    residual = (norm(Xs * M[s, s] + M[s, s] * transpose(Xs) + Ys) / norm(Xs))^2
    (residual = residual, state = single(state, M; trace = 1.0))
end

# the evolver `e` of the trajectory `state`; when its evolution without jumps is not unitary, refused if one of
# its jumps does not keep the state Gaussian (the superposition it would make could not be evolved)
function checked_jumps(state::GaussianState, e)
    if !(e isa JumpEvolver) || unitary_evolver(e.pre)
        return e
    end
    for d in e.decays, j in d.jumps
        if first(operator_kind(state, scalararg(j))) == :poly
            throw(ArgumentError("under an evolution without jumps that is not unitary, the Gaussian backend draws the " *
                                "jumps that keep a state Gaussian: gates exp(A) of quadratic A, projectors on " *
                                "occupations, linear operators (c, c†, σ⁻, …); $(obs_name(j)) is none of them"))
        end
    end
    e
end

# whether the evolver steps the states it was compiled for by unitaries only
unitary_evolver(pre::GaussianEvolver) =
    all(values(pre.forms)) do form
        isnothing(pre.fs) ? eltype(form[2]) <: Real :
                            all(maximum(abs, imag(k); init = 0.0) <= rounding_tol * (norm(k) + abs(c)) for (c, k) in form[1])
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
    even = even_states(state)
    pre(sectors) = GaussianEvolver(system, state.axes, ops, fs, duration, lossy, sectors, even, Dict{Any, Any}(),
                                   Dict{Any, Any}())
    quadratic = all(op -> degree(polynomial(system, state.layout, lower(system, op); even)) <= 2, ops)
    # the sectors of the parity when the whole evolver conserves it and needs them, on even states only
    whole = polynomial(system, state.layout, lower(system, reduce(+, ops)); even)
    p = pre(even && !quadratic && conserves_parity(whole))
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
function step_parities(pre::GaussianEvolver, state::GaussianState)
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
        if eltype(h) <: Real
            return lost(turned(state, rotation(pre, σ, h, dt), even), real(c) * dt)
        end
        m = nonunitary_pieces(h, dt)
        return dissipated(state, rotation(pre, σ, h, dt / m), m, even, real(c) * dt)
    end
    forms, even = compiled(pre, σ)
    if pre.lossy
        return lossy_magnus_step(pre, forms, state, time, dt, even)
    end
    R, loss = magnus_rotation(pre, forms, time, dt)
    lost(turned(state, R, even), loss)
end

# a step of a mixed state under -iH(t): the rotation of the step, as for a pure state
function evolve_step(::Union{Tdvp, ApproxW}, pre::GaussianEvolver, state::GaussianState{Mixed},
                     time::Number, dt::Number, step::Int, limits::Limits)
    if iszero(dt)
        return state
    end
    if !(dt isa Real)
        throw(ArgumentError("the Gaussian backend evolves by real time steps, not $dt"))
    end
    forms, even = compiled(pre, nothing)
    R, _ = magnus_rotation(pre, forms, time, dt)
    turned(state, R, even)
end

# a step of the ket of a trajectory under an evolver depending on time, in pieces whose rotation R grows vectors
# by at most 3: ‖R‖₁‖R‖∞ ≤ 9 bounds ‖R‖₂², the condition of a complex orthogonal R, which grows as e^{Γt}
# under losses Γ; a piece is halved until it does
function lossy_magnus_step(pre::GaussianEvolver, forms, state::GaussianState{<:Ket}, time::Real, dt::Real, even::Bool)
    s, h = 0.0, dt
    while dt - s > 4 * eps() * dt
        h = min(h, dt - s)
        R, loss = magnus_rotation(pre, forms, time + s, h)
        if eltype(R) <: Real
            state = lost(turned(state, R, even), loss)
        elseif opnorm(R, 1) * opnorm(R, Inf) > 9
            h /= 2
            if h < 1e-9 * dt
                throw(ErrorException("the losses of the evolver grow too fast for the steps of a GaussianState: " *
                                     "take smaller steps"))
            end
            continue
        else
            state = dissipated(state, R, 1, even, loss)
        end
        s += h
    end
    state
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
            throw(ArgumentError(isnothing(state.sum) ? nonunitary_parity_what : nonunitary_sum_what))
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
ground_state(::Any, ::GaussianState{Mixed}; kwargs...) =
    error("ground_state searches a pure state: search it pure and mix it")

# the trajectory ψ, normalised, taken to e^{-βH/2}ψ normalised in `nsteps` steps of `beta / nsteps` (or of
# `beta_step`) by the non-unitary steps of the ket of a trajectory, `log_trace` the logarithm of its squared
# norm, as thermalize of a State or an ExactState; from a later `step`, not normalised again; `algo`, `limits`,
# `mpo_limits` and `total_beta` are without use, its steps exact
function thermalize(hamiltonian, beta::Real, state::GaussianState{Trajectory};
                    algo::EvolutionAlgorithm = Tdvp(), nsteps::Union{Nothing, Int} = nothing,
                    beta_step::Union{Nothing, Real} = nothing, limits::Limits = Limits(),
                    mpo_limits::Limits = Limits(), step::Int = 1, total_beta::Union{Nothing, Real} = nothing)
    nsteps = evolution_steps("thermalize", beta, nsteps, beta_step, "beta_step")
    if step == 1
        nonzero_trace(state, "thermalize")
        state = scaled(state, 1 / sqrt(state.trace))
    end
    if nsteps == 0
        return (log_trace = 0.0, state)
    end
    if !(hamiltonian isa IndexedOp{IsOperator})
        throw(ArgumentError("thermalize of a GaussianState takes the Hamiltonian as an operator"))
    end
    pre = gaussian_evolver(state, [-0.5 * hamiltonian], nothing; lossy = true)
    log_trace = 0.0
    for k in 0:nsteps-1
        state = evolve_step(Tdvp(), pre, state, 0.0, beta / nsteps, step + k, limits)
        if !(state.trace > 0)
            throw(ErrorException("thermalize left a trajectory of norm zero at step $(step + k)"))
        end
        log_trace += log(state.trace)
        state = scaled(state, 1 / sqrt(state.trace))
    end
    (; log_trace, state)
end

thermalize(::Any, ::Real, ::GaussianState{Pure}; kwargs...) =
    error("a thermal state needs a mixed representation or trajectories: draw trajectories from a mixed one, " *
          "as State{Trajectory}(system, \"FullyMixed\")")

"""
    thermalize(hamiltonian, beta, state::GaussianState{Mixed}; nsteps | beta_step, step = 1, …)

`(log_trace, state)`: the mixed state ρ, normalised, taken to e^{-βH/2} ρ e^{-βH/2} normalised, exactly, and the
logarithm of its trace, as for a State; from a later `step`, its trace is not divided out but enters
`log_trace`, as for a State. ρ is held by a pure Gaussian purification whose first half follows e^{-βH/2} (the
product rule of the covariances of ρ and e^{-βH}, Surace and Tagliacozzo 2022, eqs. (86), (153), needs I - M
invertible, false for pure modes and near it at large β).

- `nsteps`, `beta_step`: the steps, of no effect on the result but by rounding
- `step`: the number of the first step
- `algo`, `limits`, `mpo_limits`, `total_beta`: without use

# Examples

    r = thermalize(H, 2.0, GaussianState(State{Mixed}(sys, "FullyMixed")); nsteps = 1)
"""
function thermalize(hamiltonian, beta::Real, state::GaussianState{Mixed};
                    algo::EvolutionAlgorithm = Tdvp(), nsteps::Union{Nothing, Int} = nothing,
                    beta_step::Union{Nothing, Real} = nothing, limits::Limits = Limits(),
                    mpo_limits::Limits = Limits(), step::Int = 1, total_beta::Union{Nothing, Real} = nothing)
    nsteps = evolution_steps("thermalize", beta, nsteps, beta_step, "beta_step")
    if !isfinite(beta)
        throw(ArgumentError("thermalize takes a finite beta, not $beta"))
    end
    t = nonzero_trace(state, "thermalize")
    if !(t > 0)
        throw(ArgumentError("thermalize takes a state of positive trace, not $t"))
    end
    if step == 1
        state = single(state, state.cov; trace = 1.0)
    end
    if nsteps == 0
        return (log_trace = 0.0, state)
    end
    if !(hamiltonian isa IndexedOp{IsOperator})
        throw(ArgumentError("thermalize of a GaussianState takes the Hamiltonian as an operator"))
    end
    c0, k0, _ = gaussian_generator(state.system, -0.5 * hamiltonian; axes = state.axes, even = false)
    c, k = unitary_form(c0, k0, ""; lossy = true)
    if couples_ancilla(k)
        throw(ArgumentError(mixed_odd_what))
    end
    n = size(state.cov, 1)
    db = beta / nsteps
    m = nonunitary_pieces(k, db)
    R = Matrix{ComplexF64}(I, 2n, 2n)
    R[1:n, 1:n] = step_operator(k, db / m)
    Γ = purified_covariance(state.cov)
    log_trace = log(state.trace)
    for s in 0:nsteps-1
        Γ, lognorm = nonunitary_step(Γ, R, m)
        log_trace += lognorm + 2 * real(c) * db
        if !isfinite(log_trace)
            throw(ErrorException("the trace of thermalize left the range of the floating-point numbers at step " *
                                 "$(step + s): take a smaller beta, or a Hamiltonian without its large constant"))
        end
    end
    (; log_trace, state = single(state, antisymmetric(Γ[1:n, 1:n]); trace = 1.0))
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
the weight an operator keeps of it, as for a State. A mixed state ρ becomes O ρ O† for O a gate exp(A)
without odd term, a product of Majoranas, a projector on occupations, `C` or `dag(C)` of a mode, or a
channel of one Kraus operator among them; other operators are refused.
"""
apply(op::IndexedOp{IsOperator}, state::GaussianState{Pure}; limits::Limits = Limits()) = gaussian_applied(op, state)

apply(a::IndexedOp{IsOperator}, state::GaussianState{Trajectory}; limits::Limits = Limits(),
      rng::AbstractRNG = Random.default_rng()) = trajectory_applied(a, state, limits)

apply(a::IndexedOp{IsChannel}, state::GaussianState{Trajectory}; limits::Limits = Limits(),
      rng::AbstractRNG = Random.default_rng()) = trajectory_applied(a, state, limits, rng)

ket_applied(a::IndexedOp{IsOperator}, state::GaussianState{Trajectory}; limits::Limits = Limits()) =
    gaussian_applied(a, state)

apply(op::IndexedOp{IsOperator}, state::GaussianState{Mixed}; limits::Limits = Limits()) = gaussian_applied(op, state)

function apply(a::IndexedOp{IsChannel}, state::GaussianState{Mixed}; limits::Limits = Limits())
    if !is_pure(a)
        throw(ArgumentError("the Gaussian backend does not apply to a mixed state a channel of several Kraus " *
                            "operators yet: draw it on trajectories, CreateState{Trajectory}"))
    end
    gaussian_applied(kraus(a), state)
end

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
    kind, data = operator_kind(state, s)
    if state isa GaussianState{Mixed}
        return mixed_applied(state, kind, data)
    end
    if kind == :gate
        c0, k, _ = data
        c, h = unitary_form(c0, k, ""; lossy = true)
        if eltype(h) <: Real
            return lost(turned(state, propagator(h, 1), !couples_ancilla(k)), real(c))
        end
        # exp(A) with A not anti-Hermitian, Gaussian
        m = nonunitary_pieces(h, 1)
        return dissipated(state, step_operator(h, 1 / m), m, !couples_ancilla(k), real(c))
    elseif kind == :null
        return scaled(state, 0)
    elseif kind == :monomial && isnothing(state.sum)
        α, c = data
        # c c(α) multiplies the norm by |c|, c(α) being unitary
        return single(state, flipped_covariance(state.cov, α); parities = flipped_parities(state.parities, [α]),
                      trace = abs2(c) * state.trace)
    elseif kind == :projector && isnothing(state.sum)
        occupations, λ = data
        return scaled(projected(state, occupations; null = true), λ)
    elseif kind == :mode && isnothing(state.sum)
        # c(β) z γ_a times the projector on mode j, a = 2j - 1: c_j = γ_a n_j and c_j† = γ_a (1 - n_j)
        β, j, occupied, z, αs = data
        x = projected(state, [j => occupied]; null = true)
        return single(x, flipped_covariance(flipped_covariance(x.cov, [2j - 1]), β);
                      parities = flipped_parities(x.parities, αs), trace = abs2(z) * x.trace)
    elseif kind == :linear && isnothing(state.sum)
        β, u, αs = data
        M, w = linear_applied(state.cov, vcat(u, zeros(size(state.cov, 1) - length(u))))
        if isnothing(M)
            return scaled(state, 0)
        end
        return single(state, flipped_covariance(M, β); parities = flipped_parities(state.parities, αs),
                      trace = w * state.trace)
    end
    applied_poly(state, state_polynomial(state, lower(weaken(state.system, ()), s)))
end

# O ρ O† for the operator O of `kind` and `data` (`operator_kind`): Gaussian unitaries, products of Majoranas,
# projectors on occupations and c(β) c_j, c(β) c_j†, their Majoranas on the ancilla dropped, the ancilla of
# a mixed state held empty; refused otherwise
function mixed_applied(state::GaussianState{Mixed}, kind::Symbol, data)
    physical(α) = filter(>(2), α)
    if kind == :gate
        c0, k, _ = data
        c, h = unitary_form(c0, k, ""; lossy = true)
        if !(eltype(h) <: Real)
            throw(ArgumentError("the Gaussian backend does not apply to a mixed state a gate exp(A) whose A is " *
                                "not anti-Hermitian yet: apply it to trajectories"))
        end
        return lost(turned(state, propagator(h, 1), !couples_ancilla(k)), real(c))
    elseif kind == :null
        return scaled(state, 0)
    elseif kind == :monomial
        α, c = data
        return single(state, flipped_covariance(state.cov, physical(α)); trace = abs2(c) * state.trace)
    elseif kind == :projector
        occupations, λ = data
        return scaled(projected(state, occupations; null = true), λ)
    elseif kind == :mode
        β, j, occupied, z, _ = data
        x = projected(state, [j => occupied]; null = true)
        return single(x, flipped_covariance(flipped_covariance(x.cov, [2j - 1]), physical(β));
                      trace = abs2(z) * x.trace)
    end
    throw(ArgumentError("the Gaussian backend applies to a mixed state Gaussian unitaries, products of " *
                        "Majoranas, projectors on occupations and c, c† of a mode, not this operator." *
                        mixed_rounding_hint(kind, data, size(state.cov, 1))))
end

# the rounding hint of a refusal of `mixed_applied`, for the message only: when the operator, without its
# monomials of 1e-12 of the largest or less, would be one monomial, a projector or c, c† of a mode
function mixed_rounding_hint(kind::Symbol, data, m::Int)
    if kind == :linear
        return rounding_hint(mode_deviation(data[2]), norm(data[2]))
    elseif kind != :poly
        return ""
    end
    scale = maximum(abs, values(data.terms))
    small(c) = abs(c) <= 1e-12 * scale
    kept = MajoranaPoly(Dict(α => c for (α, c) in data.terms if !small(c)))
    removed = maximum((abs(c) for c in values(data.terms) if small(c)); init = 0.0)
    linear = monomial_linear(kept, m)
    taken = length(kept) == 1 || !isnothing(occupation_projector(kept)) ||
            (!isnothing(linear) && mode_deviation(last(linear)) <= 1e-12 * norm(last(linear)))
    taken && removed > 0 ? rounding_hint(removed, scale) : ""
end

# how far the coefficients u of a linear operator are from those of c or c† of one mode, z(γ_a ± iγ_{a+1})
function mode_deviation(u::AbstractVector)
    best = Inf
    for a in 1:2:length(u)-1
        rest = norm(u[setdiff(eachindex(u), (a, a + 1))])
        best = min(best, rest + min(abs(u[a + 1] - im * u[a]), abs(u[a + 1] + im * u[a])))
    end
    best
end

# how the Gaussian backend applies the operator s exactly to a state of one component: (:gate, generator) for
# exp(A) with A quadratic, (:null, nothing) for zero, (:monomial, (α, c)), (:projector, (occupations, λ)) for
# λ times a product of projectors on occupations, (:mode, (β, j, occupied, z, αs)) for c(β) z c_j or c(β) z c_j†,
# (:linear, (β, u, αs)) for c(β) Σ_a u_a γ_a, αs the monomials, and (:poly, p) otherwise, the polynomial way
function operator_kind(state::GaussianState, s)
    generator = quadratic_gate(state, s)
    if !isnothing(generator)
        return :gate, generator
    end
    p = state_polynomial(state, lower(weaken(state.system, ()), s))
    if length(p) == 0
        return :null, nothing
    elseif length(p) == 1
        return :monomial, only(p.terms)
    end
    projector = occupation_projector(p)
    if !isnothing(projector)
        return :projector, projector
    end
    linear = monomial_linear(p, maximum((isempty(α) ? 0 : last(α) for α in keys(p.terms)); init = 0))
    if isnothing(linear)
        return :poly, p
    end
    β, u = linear
    αs = collect(keys(p.terms))
    nonzero = findall(!iszero, u)
    if length(nonzero) == 2 && isodd(nonzero[1]) && nonzero[2] == nonzero[1] + 1
        a, b = nonzero
        if u[b] == im * u[a] || u[b] == -im * u[a]
            return :mode, (β, (a + 1) ÷ 2, u[b] == im * u[a], 2 * u[a], αs)
        end
    end
    :linear, (β, u, αs)
end

# the generator of the gate s, exp(A) with A quadratic, as `gaussian_generator` gives it; nothing for another
# operator or a gate whose A is not quadratic
function quadratic_gate(state::GaussianState, s)
    A = gate_generator(s)
    if isnothing(A)
        return nothing
    end
    try
        gaussian_generator(state.system, A; axes = state.axes, even = even_states(state))
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
function scaled(state::GaussianState, c::Number)
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

# a density matrix times c, as for a State: its trace times c, real
function density_scaled(state::GaussianState{Mixed}, c::Number)
    if !isreal(c)
        throw(ArgumentError("a mixed Gaussian state is multiplied by a real number, not $c"))
    end
    single(state, state.cov; trace = real(c) * state.trace)
end

(c::Number * state::GaussianState{Mixed}) = density_scaled(state, c)
(state::GaussianState{Mixed} * c::Number) = density_scaled(state, c)
(state::GaussianState{Mixed} / c::Number) = density_scaled(state, inv(c))

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

applied_poly(::GaussianState{Mixed}, ::MajoranaPoly; cause::String = "the operator") =
    throw(ArgumentError("$cause is not Gaussian on a mixed state of the Gaussian backend"))

# the modes and occupations of an operator that is λ times a product of projectors on the occupations of
# modes, Π_m (1/2 + σ_m (i/2) γ_{2m-1}γ_{2m}), σ_m = 1 for occupied (-iγγ = 1 - 2n), as (m => occupied, λ);
# nothing for another operator. Decided by the exact equality of the polynomials (the factors of powers of
# two and ±i are exact): one whose expansion carries rounding takes the general path, exact too.
function occupation_projector(p::MajoranaPoly)
    singles = [α for α in keys(p.terms) if length(α) == 2 && isodd(α[1]) && α[2] == α[1] + 1]
    r = length(singles)
    if r == 0 || length(p) != 2^r || !haskey(p.terms, Int[])
        return nothing
    end
    λ = p.terms[Int[]] * 2^r
    if !all(p.terms[α] == λ * im / 2^r || p.terms[α] == -λ * im / 2^r for α in singles)
        return nothing
    end
    occupied = [(α[1] + 1) ÷ 2 => p.terms[α] == λ * im / 2^r for α in singles]
    product = λ * prod(MajoranaPoly(Int[] => 0.5, [2m - 1, 2m] => (o ? im : -im) / 2) for (m, o) in occupied)
    product == p ? (occupied, λ) : nothing
end

# the probability of finding mode m empty (s = 1) or occupied (s = -1) in a state of one component of
# covariance M, precise for a pure state (`occupation_probability`); for a mixed one (1 + s M_{2m-1,2m})/2
# (Bravyi and König 2012, eq. (20)), its relative precision eps/p, M_{2m-1,2m} known to eps only: below
# rounding_tol, as `measurement` takes it, zero
occupation_of(::GaussianState, M, m, s) = occupation_probability(M, m, s)

function occupation_of(::GaussianState{Mixed}, M, m, s)
    q = (1 + s * M[2m - 1, 2m]) / 2
    q <= rounding_tol ? 0.0 : q
end

# the state projected on these occupations of modes, exactly, not normalised, the covariance by Wick's
# theorem, for a mixed state too (Bravyi and König 2012, eq. (22)); when that has probability zero, an
# error, or with `null` on one component the null state
function projected(state::GaussianState, occupations; null::Bool = false)
    if isnothing(state.sum)
        M = copy(state.cov)
        t = state.trace
        for (m, o) in occupations
            s = o ? -1.0 : 1.0
            q = occupation_of(state, M, m, s)
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

function normalized(state::GaussianState)
    nonzero_trace(state, "project_on")
    isnothing(state.sum) ? single(state, state.cov; trace = 1.0) : scaled(state, 1 / sqrt(state.trace))
end

"""
    project_on(state::GaussianState, p; limits)

the state projected by `p` and normalised: exactly and without new components for a product of
projectors on the occupations of modes (N of a fermion, (1 ± σ)/2 along the axis of a qubit), as a
polynomial otherwise, see `apply`, which a mixed state refuses. On a mixed state, a result of
probability p is projected with a relative error up to about 1e-16/p, its covariance known to 1e-16
only; a probability below 1e-13 is taken as zero, as `collapse` takes it.
"""
function project_on(state::GaussianState, p::IndexedOp; limits::Limits = Limits())
    poly = state_polynomial(state, lower(weaken(state.system, ()), p))
    projector = occupation_projector(poly)
    if !isnothing(projector)
        return normalized(projected(state, first(projector)))
    end
    cause = "this measurement, outside the occupations of the modes,"
    normalized(applied_poly(state, poly; cause))
end

project_count(::GaussianState, sites, locals, T; limits::Limits = Limits()) =
    throw(ArgumentError("the Gaussian backend does not project on a total of operators of one site for now; " *
                        "measure the sites one by one"))

# drawn from the probabilities, without projecting the state
sample(state::GaussianState, pos::Int; rng = Random.default_rng()) =
    round(Int, sample(state, Level(pos); rng))

# the sites drawn one after the other, each collapsing the state
function sample(state::GaussianState; rng = Random.default_rng())
    out = Int[]
    for pos in eachindex(state.system.sites)
        r = collapse(state, pos; rng)
        push!(out, r.result)
        state = r.state
    end
    out
end
