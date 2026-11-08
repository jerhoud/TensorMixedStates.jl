# The solvers of an exact state: its exponentials by Krylov, for its evolutions, those depending
# on time by Magnus, and its thermalization, its searches by Lanczos and Arnoldi, and its random
# draws.

"""
the operator `a`, or its superoperator, see `superoperator`, prepared to act on the exact states
of `system`, and its compaction, which bounds its norm and its width and checks its flux
"""
function exact_acting(::Type{Pure}, system::System, a::IndexedOp{IsOperator})
    c = compacted_op(system, a)
    plain = weaken(system, ())
    return prepared_acted(plain, gathered(plain, a)), c
end

function exact_acting(::Type{Mixed}, system::System, a::IndexedOp)
    s = superoperator(a)
    c = compacted_super(system, s)
    return prepared_channel(weaken(system, ()), s), c
end

"""
    ExactEvolver

the evolver of an exact state: the linear map of each part and a bound of its norm, see
`exact_acting`, the time functions, see `split_terms`, and the duration of the whole evolution
"""
struct ExactEvolver
    actions::Vector{Function}
    bounds::Vector{Float64}
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
        return exact_acting(R, system, R === Pure ? followed_operator(a) : a)
    end
    return ExactEvolver(first.(parts), [ c.norm for (_, c) in parts ], fs, duration)
end

prepared_evolver(::ExactState, ::Union{PreMPO, MPO}, ::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    error("an ExactState is given the evolver rather than its PreMPO or MPO: " *
          "evolve(evolver, duration, state)")

prepared_evolver(::ExactState, x, ::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) = error("$(repr(x)) is not an evolver")

"""
the exponential ``e^{tA}x`` of the linear map `A` of norm at most `bound`, by Krylov to `tol` per
unit of time, see `KrylovKit.exponentiate`, as an array of norm one and the logarithm of its norm
"""
function exp_scaled(A, bound::Real, t::Number, x::AbstractArray, tol::Real)
    n = norm(x)
    if iszero(n)
        return (x, -Inf)
    end
    # KrylovKit holds its precision as absolute, so that it must stay relative: in pieces,
    # renormalised in between, of `A` shifted by its mean value `μ`, along which the norm of a
    # normal map, convex in log and stationary at the start, only grows
    m = max(1, ceil(Int, abs(t) * bound / 10))
    y, scale = x / n, log(n)
    for _ in 1:m
        μ = dot(y, A(y))
        z, info = KrylovKit.exponentiate(v -> A(v) - μ * v, t / m, y; tol, ishermitian = false,
                                         eager = true, verbosity = 0)
        if info.converged == 0
            error("an exponential of an ExactState did not converge in Krylov: take more steps")
        end
        nz = norm(z)
        if iszero(nz)
            return (z, -Inf)
        end
        y, scale = z * cis(imag(t / m * μ)) / nz, scale + log(nz) + real(t / m * μ)
    end
    return (y, scale)
end

"""
the exponential ``e^{tA}x``, see `exp_scaled`, refused when its norm leaves the floats
"""
function exact_exp(A, bound::Real, t::Number, x::AbstractArray, tol::Real)
    if iszero(norm(x))
        return x
    end
    y, scale = exp_scaled(A, bound, t, x, tol)
    if !(log(floatmin(Float64)) < scale < log(floatmax(Float64)))
        error("the evolution takes the norm of the state to e^$(round(scale; digits = 1)), " *
              "beyond the floats: evolve it in parts, normalizing it in between")
    end
    return exp(scale) * y
end

"""
the array `x` evolved from `time` for `dt` under `cs -> Σ_j c_j A_j` that `act` gives for the
values `cs` of the time functions `fs`, of norm at most `bound(cs)`, by the commutator-free
Magnus integrator of order 4, see `step_coefs`, its error over `duration` below
`1000 rounding_tol`
"""
function magnus_evolved(act, bound, fs, time::Number, dt::Number, duration::Number,
                        x::AbstractArray)
    if iszero(norm(x))
        return x
    end
    # the exponentials to a tenth of the error the substeps are allowed, so that their difference
    # measures that of the integrator
    tol = 100 * rounding_tol / abs(duration)
    magnus(t, h, y) = foldl((y, cs) -> exact_exp(act(cs), bound(cs), h, y, tol),
                            step_coefs(fs, t, h, 4); init = y)
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
        # its share of the error of the whole evolution, above the rounding of the operations
        allowed = max(800 * rounding_tol * abs(h / duration), 10 * eps())
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
    if iszero(dt)
        return state
    end
    system = state.system
    act(cs) = x -> sum(c * f(x) for (c, f) in zip(cs, pre.actions))
    bound(cs) = sum(abs(c) * b for (c, b) in zip(cs, pre.bounds))
    # a step taken alone, outside evolve, is the whole evolution
    duration = something(pre.duration, dt)
    unit = ones(length(pre.actions))
    x = isnothing(pre.functions) ?
        exact_exp(act(unit), bound(unit), dt, state.state, 1000 * rounding_tol / abs(duration)) :
        magnus_evolved(act, bound, pre.functions, time, dt, duration, state.state)
    return ExactState{R}(system, x, Val(:raw))
end

"""
the eigenvalue of the linear map `A` of the smallest real part, `which = :SR`, or of the largest,
`:LR`, and its eigenvector, searched by Krylov from `x` to `rounding_tol`, `why` saying what keeps
it from converging
"""
function exact_eigen(A, x::AbstractArray, which::Symbol, hermitian::Bool, why::String)
    # wider than the defaults of KrylovKit, which a spectrum crowded at its end needs
    vals, vecs, info = KrylovKit.eigsolve(A, x / norm(x), 1, which; ishermitian = hermitian,
                                          tol = rounding_tol, krylovdim = min(length(x), 100),
                                          maxiter = 1000, verbosity = 0)
    if info.converged == 0
        error(why)
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

"""
a random vector of norm one of the shape of the vector `x` of `system`, in the sector of its
charges declared, drawn by `rng` of Gaussian elements of type `eltype`, so uniformly
"""
function random_in_sector(system::System, x::AbstractVector, rng::AbstractRNG,
                          eltype::Type{<:Number} = ComplexF64)
    r = randn(rng, eltype, length(x))
    if is_charged(system)
        a, idx, _ = dense_legs(system, x)
        table = charge_table(idx)
        q = common_charge(a, table, rounding_tol * norm(a))
        if isnothing(q)
            error("the state has no definite charge, which a system conserving it cannot hold")
        end
        r = vec(of_charge(reshape(r, size(a)), table, q))
    end
    return r / norm(r)
end

# as random_state of a State, but uniform whatever `linkdims`: a pure state of the uniform measure
# on the vectors of its sector, a mixed one of the measure its purification on the system doubled
# induces, as that of the MPS, the ancillas on the right taking no fermionic sign
function random_state(::Type{ExactState{Pure}}, system::System, linkdims::Int;
                      rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64)
    if is_charged(system)
        error("cannot draw a random pure state on a system that conserves something, there " *
              "being no sector to draw it in: give the local states of the sector wanted")
    end
    d = prod(site_dims(system))
    return ExactState{Pure}(system, random_in_sector(system, zeros(d), rng, eltype), Val(:raw))
end

random_state(::Type{ExactState{Pure}}, system::System, states, linkdims::Int;
             rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64) =
    random_state(ExactState{Pure}(system, states), linkdims; rng, eltype)

random_state(state::ExactState{Pure}, linkdims::Int; rng::AbstractRNG = Random.default_rng(),
             eltype::Type{<:Number} = ComplexF64) =
    ExactState{Pure}(state.system, random_in_sector(state.system, state.state, rng, eltype),
                     Val(:raw))

random_state(::ExactState{Mixed}, ::Int; kwargs...) =
    error("randomising a mixed state is not implemented")

function random_state(::Type{ExactState{Mixed}}, system::System, linkdims::Int;
                      rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64)
    refuse_strong(system)
    if is_charged(system)
        error("cannot draw a random mixed state on a system that conserves something without a " *
              "sector to start from: name the local states its purification starts from")
    end
    d = prod(site_dims(system))
    return traced_half(system, random_in_sector(system ⊗ system, zeros(d^2), rng, eltype))
end

function random_state(::Type{ExactState{Mixed}}, system::System, states, linkdims::Int;
                      rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64)
    refuse_strong(system)
    doubled = system ⊗ system
    start = ExactState{Pure}(doubled, double(states, length(system)))
    return traced_half(system, random_in_sector(doubled, start.state, rng, eltype))
end

random_state(::Type{ExactState{Mixed}}, ::System, ::Union{Mixture, FullyMixed}, ::Int;
             kwargs...) =
    error("a random mixed state is drawn from the purification of a product state: give its " *
          "local states")

"""
the density matrix of `system` the pure vector `ψ` of `system ⊗ system` leaves once the second
half, on the right, is traced out
"""
function traced_half(system::System, ψ::AbstractVector)
    d = prod(site_dims(system))
    # ψ in the order of kron, the first half the slowest: m[a, i] = ψ(i, a)
    m = reshape(ψ, d, d)
    return ExactState{Mixed}(system, transpose(m) * conj(m), Val(:raw))
end

# the method and the limits of the MPS without use, `algo` not needed; the search of a State
# solved exactly, see `orthogonal_args`: the lowest state of the Hamiltonian plus `weight` times
# the projectors of the states of `orthogonal_to`, its energy that of the Hamiltonian alone
function ground_state(hamiltonian, state::ExactState{Pure};
                      algo::Union{Nothing, SearchAlgorithm} = nothing, limits::Limits = Limits(),
                      mpo_limits::Limits = Limits(), orthogonal_to::AbstractVector = [],
                      weight::Union{Nothing, Real} = nothing)
    system = state.system
    H, c = exact_acting(Pure, system, searched_operator(system, hamiltonian, "ground_state"))
    refs = map(orthogonal_to) do m
        if !(m isa Union{ExactState{Pure}, State{Pure}})
            error("ground_state searches a state orthogonal to pure states, not to a $(typeof(m))")
        end
        e = m isa ExactState ? m : ExactState(m)
        same_sites(state, e)
        return normalize(e)
    end
    w = something(weight, 1.1 * c.width)
    A = isempty(refs) ? H : v -> H(v) + w * sum(r.state * dot(r.state, v) for r in refs)
    # Lanczos started from the state given perturbed by a random state of its sector of the
    # charges declared, so that it reaches the sectors of the symmetries not declared, as DMRG
    # does, the state given still picking among degenerate states
    x = state.state + 0.01 * norm(state.state) * random_in_sector(system, state.state, Xoshiro(1))
    _, v = exact_eigen(A, x, :SR, true, "ground_state of an ExactState did not " *
                       "converge in Krylov, its spectrum crowded at its bottom: start from a " *
                       "state closer to the ground state")
    found = normalize(ExactState{Pure}(system, v, Val(:raw)))
    check_orthogonal(found, (; orthogonal_to = refs, weight = w), weight)
    return (energy = real(dot(v, H(v))) / real(dot(v, v)), state = found)
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
    L, _ = exact_acting(Mixed, system, a)
    _, x = exact_eigen(L, state.state, :LR, false,
                       "steady_state of an ExactState did not converge in Krylov, the spectrum " *
                       "of the Lindbladian crowded near zero: start from a state closer to " *
                       "the steady state, one evolved for a while for instance")
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
# MPS without use; `total_beta`, that of the whole thermalization this call is a part of, bounds
# its error
function thermalize(hamiltonian, beta::Real, state::ExactState{Mixed};
                    algo::EvolutionAlgorithm = Tdvp(), nsteps::Union{Nothing, Int} = nothing,
                    beta_step::Union{Nothing, Real} = nothing, limits::Limits = Limits(),
                    mpo_limits::Limits = Limits(), step::Int = 1,
                    total_beta::Union{Nothing, Real} = nothing)
    nsteps = evolution_steps("thermalize", beta, nsteps, beta_step, "beta_step")
    if nsteps == 0
        return (log_trace = 0., state = step == 1 ? normalize(state) : state)
    end
    system = state.system
    a = searched_operator(system, hamiltonian, "thermalize")
    # ρ ↦ Hρ + ρH, whose exponential at -β/2 is e^{-βH/2} ρ e^{-βH/2}
    K, c = exact_acting(Mixed, system, a)
    bound = c.norm
    x = (step == 1 ? normalize(state) : state).state
    log_trace = 0.
    for k in 0:nsteps - 1
        # the trace taken apart from the norm, which the exponential may take far beyond floats
        y, scale = exp_scaled(K, bound, -beta / nsteps / 2, x,
                              2000 * rounding_tol / abs(something(total_beta, beta)))
        t = real(tr(y))
        if !(t > 0)
            error("the trace of the state reached $t at step $(step + k) of thermalize")
        end
        x = y / t
        log_trace += scale + log(t)
    end
    return (; log_trace, state = ExactState{Mixed}(system, x, Val(:raw)))
end

thermalize(::Any, ::Real, ::ExactState{Pure}; kwargs...) =
    error("a thermal state needs a mixed representation: mix the state first")
