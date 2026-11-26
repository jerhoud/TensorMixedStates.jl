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

# the physical parity of a state of one component and definite parity, nothing for another: that of the
# ancilla on the even states, ⟨-iγ₁γ₂⟩ = M₁₂ = ±1
state_parity(state::GaussianState{Pure}) =
    isnothing(state.sum) && state.definite ? Int(sign(state.cov[1, 2])) : nothing

"""
the evolver of a Gaussian state under L(t) = Σ_j f_j(t) L_j, L_j = c_j + (1/4) Σ k_j γγ with k_j real:
the time functions, nothing for a constant evolver, the duration of the evolution, and whether L keeps
the physical parity
"""
struct GaussianTimeEvolver
    ks::Vector{Matrix{Float64}}
    fs::Vector
    duration::Union{Nothing, Number}
    even::Bool
end

function prepared_evolver(state::GaussianState{Pure}, terms::AbstractVector, ::Limits = Limits();
                          duration::Union{Nothing, Number} = nothing)
    ops, fs = split_terms(terms)
    if !all(op -> op isa IndexedOp{IsOperator}, ops)
        throw(ArgumentError("the Gaussian backend evolves a pure state under an operator, -iH"))
    end
    what = "the Gaussian backend evolves a pure state under -iH with H Hermitian, for now"
    parity = state_parity(state)
    if isnothing(fs)
        c, k, _ = gaussian_generator(state.system, reduce(+, ops); axes = state.axes, parity)
        _, h = unitary_form(c, k, what)
        return GaussianEvolver(h, Dict{Float64, Matrix{Float64}}(), !couples_ancilla(k))
    end
    forms = [gaussian_generator(state.system, op; axes = state.axes, parity) for op in ops]
    ks = [last(unitary_form(c, k, what)) for (c, k, _) in forms]
    GaussianTimeEvolver(ks, fs, duration, !any(couples_ancilla(k) for (_, k, _) in forms))
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

function evolve_step(::Union{Tdvp, ApproxW}, pre::GaussianTimeEvolver, state::GaussianState{Pure},
                     time::Number, dt::Number, step::Int, limits::Limits)
    if iszero(dt)
        return state
    end
    if !(dt isa Real)
        throw(ArgumentError("the Gaussian backend evolves by real time steps, not $dt"))
    end
    turned(state, magnus_rotation(pre, time, dt), pre.even)
end

# the rotation R of the step from `time` to `time + dt`, dR/dt = (Σ_j f_j(t) k_j) R, by the
# commutator-free Magnus integrator of order 4 with substeps, as `magnus_evolved` for an ExactState:
# a substep against its two halves, its error kept within its share of 1000 rounding_tol over the
# duration
function magnus_rotation(pre::GaussianTimeEvolver, time::Number, dt::Number)
    duration = something(pre.duration, dt)
    function generator(cs)
        if !all(c -> isreal(c), cs)
            throw(ArgumentError("the Gaussian backend evolves under real time functions, a unitary evolution"))
        end
        sum(real(c) * k for (c, k) in zip(cs, pre.ks))
    end
    magnus(t, h, R) = foldl((R, cs) -> propagator(generator(cs), h) * R, step_coefs(pre.fs, t, h, 4); init = R)
    R = Matrix{Float64}(I, size(first(pre.ks)))
    s, f = 0.0, 1.0
    while 1 - s > 4 * eps()
        f = min(f, 1 - s)
        t, h = time + s * dt, f * dt
        whole = magnus(t, h, R)
        halves = magnus(t + h / 2, h / 2, magnus(t, h / 2, R))
        err = norm(halves - whole) / norm(halves) / 15
        if isnan(err)
            error("the evolver of a GaussianState gave NaN within a step")
        end
        allowed = max(800 * rounding_tol * abs(h / duration), 10 * eps())
        if err ≤ allowed
            R, s = halves, s + f
        end
        f *= clamp(0.9 * (allowed / max(err, eps()))^(1 / 5), 0.25, 2.0)
        if f < 1e-9
            error("an evolver depending on time varies too fast for the substeps of a GaussianState: " *
                  "take smaller steps")
        end
    end
    R
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
    degenerate_info(g)
    cov = zeros(M, M)
    cov[1, 2], cov[2, 1] = σ, -σ
    cov[3:end, 3:end] = g.cov
    found = GaussianState{Pure}(system, state.axes, state.layout, cov, nothing, true, 1.0, state.max_terms)
    (energy = real(c) + g.energy, width = sum(g.energies), state = found)
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
function projected(state::GaussianState{Pure}, occupations)
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
        return single(state, M; definite = false, trace = t)
    end
    s = as_sum(state)
    for (m, o) in occupations
        project!(s, m, o)
    end
    superposed(state, s)
end

normalized(state::GaussianState{Pure}) =
    isnothing(state.sum) ? single(state, state.cov; trace = 1.0) : scaled(state, 1 / sqrt(state.trace))

"""
    project_on(state::GaussianState, p; limits)

the state projected by `p` and normalised: exactly and without new components for a product of
projectors on the occupations of modes (N of a fermion, (1 ± σ)/2 along the axis of a qubit), as a
polynomial otherwise, see `apply`.
"""
function project_on(state::GaussianState{Pure}, p::IndexedOp; limits::Limits = Limits())
    poly = polynomial(state.system, state.layout, lower(weaken(state.system, ()), p))
    occupations = occupation_projector(poly)
    normalized(isnothing(occupations) ? apply(p, state; limits) : projected(state, occupations))
end

project_count(::GaussianState, sites, locals, T; limits::Limits = Limits()) =
    throw(ArgumentError("the Gaussian backend does not project on a total of operators of one site, which " *
                        "is not Gaussian; measure the sites one by one"))

sample(state::GaussianState{Pure}, pos::Int; rng = Random.default_rng()) =
    first(collapse(state, pos; rng))

# the sites drawn one after the other, each collapsing the state
function sample(state::GaussianState{Pure}; rng = Random.default_rng())
    out = Int[]
    for pos in eachindex(state.system.sites)
        r = collapse(state, pos; rng)
        push!(out, r.result)
        state = r.state
    end
    out
end
