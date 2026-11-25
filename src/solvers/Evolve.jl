# The evolution of a state in time: evolve, a loop over the steps that evolve_step makes for each
# method, Tdvp and ApproxW, with the evolver prepared once and its MPOs kept, see cached_mpos.

export evolve

"""
    truncated_once(pre, limits)

`pre` truncated by `limits`, see `truncated`, kept in `pre` for the calls that follow, see
`cached_mpos`
"""
truncated_once(pre::PreMPO, limits::Limits) =
    limits == Limits() ? pre : cached_mpos(() -> truncated(pre, limits), pre, :truncated, limits)

"""
    check_pre(pre, state)

refuse a `PreMPO` prepared on another system than that of `state`, or for the matrices of
another representation, see `algebra`
"""
function check_pre(pre::PreMPO{S}, state::State{R}) where {S, R}
    if S !== algebra(R)
        error("the PreMPO was prepared on a state of the other representation: prepare it with " *
              "the state, PreMPO(state, op)")
    elseif pre.system ≠ state.system
        error("the PreMPO was prepared on another system than that of the state, " *
              "$(system_difference(pre.system, state.system)): prepare it with the state, " *
              "PreMPO(state, op)")
    end
    return pre
end

"""
    prepared_operator(state, op, mpo_limits = Limits())

the `PreMPO` of the operator `op` of a solver on `state`, or `op` itself when it is one,
truncated by `mpo_limits`, see `truncated_once`, refused when it depends on time, which an
evolution alone follows, see `prepared_evolver`
"""
prepared_operator(state::State, a::IndexedOp, mpo_limits::Limits = Limits()) =
    PreMPO(state, a; limits = mpo_limits)
function prepared_operator(state::State, pre::PreMPO, mpo_limits::Limits = Limits())
    check_pre(pre, state)
    if !isnothing(pre.functions)
        error("the PreMPO depends on time, which an evolution alone follows")
    end
    return truncated_once(pre, mpo_limits)
end

prepared_operator(::State, ::MPO, ::Limits = Limits()) =
    error("the solvers take an operator or its PreMPO, which keeps the MPOs they build, not an MPO")

prepared_operator(::State, x, ::Limits = Limits()) =
    error("$(repr(x)) is not an operator or its PreMPO")

"""
    prepared_evolver(state, evolver, mpo_limits = Limits(); duration = nothing)

the `PreMPO` of the evolver of `state`, see `prepared_operator`, with its time functions: an
operator, a vector of terms `[A, B => g]`, a term `A => g`, or the `PreMPO` of one of them.
`duration` is that of the whole evolution, which an MPS does not use. A state of one's own adds
a method returning what its `evolve_step` takes as `pre`.
"""
prepared_evolver(state::State, a::Union{IndexedOp, PreMPO, MPO}, mpo_limits::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    prepared_operator(state, a, mpo_limits)
prepared_evolver(state::State, pre::PreMPO, mpo_limits::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    truncated_once(check_pre(pre, state), mpo_limits)
prepared_evolver(state::State, terms::Union{AbstractVector, Pair}, mpo_limits::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    PreMPO(state, terms; limits = mpo_limits)

prepared_evolver(::State, x, ::Limits = Limits(); duration::Union{Nothing, Number} = nothing) =
    error("$(repr(x)) is not an evolver")

"""
    restored(e, st)

the MPS `st` written in the basis of the MPS `e`, made orthonormal from the right: the bond
spaces of `e` holding `st` as it was, which `expand` alters when the vectors it adds to a basis
are not quite orthogonal to it
"""
function restored(e::MPS, st::MPS)
    b = orthogonalize(e, 1)
    env = ITensor(1.)
    for j in length(st):-1:2
        env = env * st[j] * dag(b[j])
    end
    return MPS([ st[1] * env ; [ b[j] for j in 2:length(st) ] ])
end

"""
    expanded(state, mpo, maxdim)

the state enlarged by a global Krylov expansion under `mpo`, its Krylov vectors computed as `apply`
computes a product, "naive", truncated at rounding and to twice the bond dimension `maxdim` of
the step, which bounds the cost of the second, the product of `mpo` by the first, and leaves new
directions even to a state just below `maxdim`: truncated to one more than the bond dimension
of the state, the default of ITensorMPS, they keep nothing of some charged states, as
FullyMixed. The other options are those of ITensorMPS. The state is then written again in the
bonds the expansion enlarged, see `restored`. An expansion only enlarges the basis of the step,
which goes on without it, with a warning, when ITensorMPS fails on charged indices, as on a state
that `mpo` maps onto itself, and when the state written again is off by more than the precision
ITensorMPS gives the expansion, the square root of that of the numbers
"""
function expanded(state::State, mpo::MPO, maxdim::Int)
    st = state.state
    e = try
        expand(st, mpo; alg = "global_krylov",
               apply_kwargs = (; alg = "naive", truncation(Limits())...,
                               maxdim = twice_maxdim(maxdim)))
    catch err
        if !any(hasqns, st)
            rethrow()
        end
        @warn "the expansion before a tdvp step failed on charged indices " *
              "($(nameof(typeof(err)))), as on a state the evolver maps onto itself: the step " *
              "goes on without it"
        return state
    end
    new = State(state, restored(e, st))
    δ = norm(new - state) / norm(state)
    if δ > sqrt(eps())
        δ = round(δ; sigdigits = 2)
        @warn "the expansion before a tdvp step changed the state by $δ, beyond its precision: " *
              "the step goes on without it"
        return state
    end
    return new
end

# every bond of `st` as large as the space on its smaller side, which an expansion cannot
# enlarge; charges make the largest bonds smaller, and such states are still expanded
function bonds_full(st::MPS)
    d = [ dim(siteinds(st, i)) for i in eachindex(st) ]
    cap = typemax(Int) ÷ maximum(d)
    right = ones(Int, length(st))
    for i in length(st)-1:-1:1
        right[i] = min(right[i + 1] * d[i + 1], cap)
    end
    left = 1
    for b in 1:length(st)-1
        left = min(left * d[b], cap)
        if linkdim(st, b) < min(left, right[b])
            return false
        end
    end
    return true
end

"""
    evolve_step(algo, pre, state, time, dt, step, limits)

the state after one step of `dt` from `time` by the method `algo`, under the evolver `pre`,
whose time functions are `pre.functions`, one per term, `nothing` when it does not depend on
time, see `make_mpo`, `step` being the number of the step, from 1, and `limits` its limits, see
`step_limits`. An `EvolutionAlgorithm` of one's own adds a method, and so does a state of
one's own, `pre` being what its `prepared_evolver` returned.
"""
function evolve_step(algo::Tdvp, pre::PreMPO, state::State, time::Number, dt::Number, step::Int,
                     limits::Limits)
    # an evolution leaves the null state null, which ITensors fails on
    if vanished(state.state)
        return state
    end
    # KrylovKit takes its precision as absolute, and a vector below it for a fixed point: the
    # step is made on the state of norm one, its norm put back after, the evolution being linear.
    # The norm is the full product, not that of an orthogonality centre, which a state read back
    # from a file has lost, so that a resumed run takes the same steps
    n = sqrt(real(dot(state.state, state.state)))
    if iszero(n)
        return state
    end
    state = State(state, state.state / n)
    fs = pre.functions
    mpo = isnothing(fs) ? cached_mpos(() -> make_mpo(pre), pre, :mpo, nothing) :
                          make_mpo(pre, only(step_coefs(fs, time, dt, 2)))
    tr = truncation(limits)
    # before the step: from a state of small bond dimension, a product state, an expansion after
    # the step comes too late, which leaves an error of order dt; a single site has no bond
    small = maxlinkdim(state) < min(5, tr.maxdim)
    if length(state) > 1 && (small || due(algo.expand_period, step - 1)) && !bonds_full(state.state)
        state = expanded(state, mpo, tr.maxdim)
    end
    st = if length(state) == 1
        # tdvp of ITensorMPS takes two sites at least
        w = mpo[1]
        MPS([ noprime(exp(dt * w, inds(w; plev = 1), inds(w; plev = 0)) * state.state[1]) ])
    else
        # eager: KrylovKit tests convergence after every vector rather than once its space is
        # full, which the short local steps of tdvp reach with far fewer products
        tdvp(mpo, dt, state.state; nsweeps = 1, tr...,
             updater_kwargs = (; eager = true, krylov_kwargs(algo.krylov)...))
    end
    state = State(state, n * st)
    if due(algo.hermitianize_period, step)
        state = hermitianize(state; limits)
    end
    return state
end

function evolve_step(algo::ApproxW, pre::PreMPO, state::State, time::Number, dt::Number,
                     step::Int, limits::Limits)
    if vanished(state.state)
        return state
    end
    mpos = if isnothing(pre.functions)
        cached_mpos(() -> approx_W_mpos(pre, dt, algo), pre, :approx_W, (algo.order, algo.w, dt))
    else
        reduce(vcat, [ approx_W_mpos(pre, dt, algo, c)
                       for c in step_coefs(pre.functions, time, dt, algo.order) ])
    end
    st = state.state
    for mpo in mpos
        st = applied_mpo(mpo, st, limits; alg = algo.apply_algo)
    end
    state = State(state, st)
    if due(algo.hermitianize_period, step)
        state = hermitianize(state; limits)
    end
    return state
end

############### Jumps of the trajectories ###############

"""
    unravelled(a)

the operators and the dissipators of the generator `a`, `(ops, dissipators)`: the operator of
each part with no dissipator, see `pure_part`, and `rate => L` for `rate * Dissipator(L)`, `L` an
operator or a channel, placed on sites as in `a`
"""
function unravelled(a::Op{IsGenerator})
    if a isa SumOp
        parts = map(unravelled, a.subs)
        return (reduce(vcat, first.(parts)), reduce(vcat, last.(parts)))
    elseif a isa ScalarOp
        ops, ds = unravelled(a.arg)
        return (Op{IsOperator}[ scale(a.coef, x) for x in ops ],
                Pair{Float64, Op}[ a.coef * r => x for (r, x) in ds ])
    elseif a isa AtIndex
        ops, ds = unravelled(a.op)
        return (Op{IsOperator}[ AtIndex(x, a.index) for x in ops ],
                Pair{Float64, Op}[ r => AtIndex(x, a.index) for (r, x) in ds ])
    elseif a isa Dissipator
        return (Op{IsOperator}[], Pair{Float64, Op}[1.0 => a.arg])
    end
    return (Op{IsOperator}[pure_part(a)], Pair{Float64, Op}[])
end

"""
    Decay(f, rate, jumps, labels)

the jumps of a term of an evolver, which a trajectory draws: `rate`, the sum of their
``L^\\dagger L``, whose expectation value times `f(t)` (1 for `f = nothing`) is the rate of a
jump at time `t`, the operators `jumps`, and for each the dissipator it comes from, as written
and placed, and its index among the Kraus operators of that dissipator, 1 for an operator
"""
struct Decay
    f::Any
    rate::IndexedOp{IsOperator}
    jumps::Vector{IndexedOp{IsOperator}}
    labels::Vector{Tuple{IndexedOp, Int}}
end

"""
    JumpEvolver(pre, decays, keeps, exact, rng, jumps, stop)

the evolver of a trajectory under dissipators, or under operators keeping its trace: the evolver
of its ket under the effective operator ``K - \\frac{1}{2}\\sum_j L_j^\\dagger L_j``, see
`ket_evolver`, the `Decay` of each term
with jumps, none without dissipators, whether its operators `K` keep the trace, which brings the
trajectory back to it at each step unless `exact`, its jumps then weighed by the exact norm of
its ket, see `exact_norm`, the generator drawing them, see `with_rng`, `nothing` for the
default one, the list the jumps drawn are pushed to, see `evolve`, `nothing` for none, and the
end of the evolution, `nothing` when unknown, past which rounding may put the end of a step
"""
struct JumpEvolver
    pre::Any
    decays::Vector{Decay}
    keeps::Bool
    exact::Bool
    rng::Union{Nothing, AbstractRNG}
    jumps::Union{Nothing, AbstractVector}
    stop::Union{Nothing, Float64}
end

"""
    ket_evolver(state, terms, mpo_limits; duration = nothing)

the evolver of the ket of the trajectory `state` under the operators `terms`, as
`prepared_evolver` takes them, which `evolve_step` takes as `pre`, over the whole `duration` of
the evolution when it is known: a method for each type of states holding trajectories
"""
ket_evolver(state::State{Trajectory}, terms::AbstractVector, mpo_limits::Limits;
            duration::Union{Nothing, Number} = nothing) =
    PreMPO(state, terms; limits = mpo_limits)

"""
    exact_norm(state)

whether the backend of the trajectory `state` keeps the norm of its ket exactly, without
truncation: the jumps of its trajectory `:state` are then weighed by the norm the evolution
leaves, which corrects the rate taken linear over a step where it is not zero; `false` by
default
"""
exact_norm(::AbstractState) = false

function JumpEvolver(state::AbstractState{Trajectory}, ops::Vector, fs, mpo_limits::Limits;
                     duration::Union{Nothing, Number} = nothing)
    sys = state.system
    terms = []
    pures = []
    decays = Decay[]
    for (k, a) in enumerate(ops)
        pure, ds = a isa IndexedOp{IsGenerator} ? unravelled(a) :
                   ([followed_operator(a, Trajectory)], [])
        jumps = IndexedOp{IsOperator}[]
        labels = Tuple{IndexedOp, Int}[]
        for (r, x) in ds
            if r < 0
                error("$(obs_name(a)) holds a dissipator of negative rate, $r, whose jumps no " *
                      "trajectory draws: use State{Mixed}")
            end
            ks = x isa IndexedOp{IsChannel} ? kraus_operators(sys, x) : [x]
            for (k, j) in enumerate(ks)
                # a Kraus operator proportional to the identity has a null dissipator
                if !(scalararg(j) isa IdentityOp)
                    push!(jumps, sqrt(r) * j)
                    push!(labels, (x, k))
                end
            end
        end
        f = isnothing(fs) ? nothing : fs[k]
        op = isempty(pure) ? nothing : sum(pure)
        push!(pures, op)
        if !isempty(jumps)
            for j in jumps
                # refused here, off the system or of no definite flux, as apply refuses it,
                # rather than at the first jump drawn
                check_indices(sys, j)
                pure_steps(sys, j)
                check_jump(sys, j, a)
            end
            rate = sum(j -> dag(j) * j, jumps)
            push!(decays, Decay(f, rate, jumps, labels))
            op = isnothing(op) ? -0.5 * rate : op - 0.5 * rate
        end
        if !isnothing(op)
            push!(terms, isnothing(f) ? op : op => f)
        end
    end
    # the operators of the terms of one function of time anti-Hermitian together keep the trace
    keeps = all(unique(isnothing(fs) ? [nothing] : fs)) do f
        ks = [ k for k in eachindex(ops) if isequal(isnothing(fs) ? nothing : fs[k], f) &&
                                            !isnothing(pures[k]) ]
        isempty(ks) || prepare_op(sys, im * sum(pures[ks])).hermitian
    end
    return JumpEvolver(ket_evolver(state, terms, mpo_limits; duration), decays, keeps,
                       exact_norm(state), nothing, nothing, nothing)
end

# refuse the jump `j` of the generator `a` when it changes a charge `system` conserves strongly,
# as State{Mixed} refuses the dissipator
function check_jump(system::System, j::IndexedOp{IsOperator}, a)
    names = strong_names(system)
    if !isempty(names)
        flux = compacted_op(system, j).flux
        if !isnothing(flux) && !isempty(strong_part(flux, names))
            strong_refusal(a, names)
        end
    end
end

prepared_evolver(state::State{Trajectory}, a::IndexedOp, mpo_limits::Limits = Limits();
                 duration::Union{Nothing, Number} = nothing) =
    prepared_evolver(state, [a], mpo_limits; duration)

prepared_evolver(state::State{Trajectory}, terms::Union{AbstractVector, Pair},
                 mpo_limits::Limits = Limits(); duration::Union{Nothing, Number} = nothing) =
    trajectory_evolver(state, terms, mpo_limits; duration)

# the evolver of the trajectory `state` under `terms`, see `prepared_evolver`
function trajectory_evolver(state::AbstractState{Trajectory}, terms::Union{AbstractVector, Pair},
                            mpo_limits::Limits; duration::Union{Nothing, Number} = nothing)
    ops, fs = split_terms(terms isa Pair ? [terms] : terms)
    if all(a -> !(a isa IndexedOp{IsGenerator}) || is_pure(a), ops)
        e = JumpEvolver(state, ops, fs, mpo_limits; duration)
        # an evolver keeping the trace brings the trajectory back to it at each step, so that what
        # truncation loses is no weight
        return e.keeps && !e.exact ? e : e.pre
    end
    return JumpEvolver(state, ops, fs, mpo_limits; duration)
end

"""
    with_rng(pre, rng, jumps = nothing; weighed = true, stop = nothing)

the evolver `pre`, see `prepared_evolver`, drawing with `rng` and pushing its jumps to the list
`jumps` when it is given, see `evolve`, over an evolution ending at `stop` when it is given, a
trajectory that carries no weight, `weighed = false`, brought back to its trace as on a backend
whose norm truncation lowers: itself when it draws nothing
"""
with_rng(pre, ::AbstractRNG, ::Union{Nothing, AbstractVector} = nothing; weighed::Bool = true,
         stop::Union{Nothing, Number} = nothing) = pre
with_rng(e::JumpEvolver, rng::AbstractRNG, jumps::Union{Nothing, AbstractVector} = nothing;
         weighed::Bool = true, stop::Union{Nothing, Number} = nothing) =
    JumpEvolver(e.pre, e.decays, e.keeps, e.exact && weighed, rng, jumps,
                isnothing(stop) ? nothing : Float64(stop))

# the rate of the jumps of `d` on `state` at `time`, taken at the end `stop` of the evolution when
# rounding puts the end of a step past it, a negative rate function refused
function decay_rate(d::Decay, state::AbstractState{Trajectory}, time::Number,
                    stop::Union{Nothing, Float64})
    t = isnothing(stop) ? time : min(time, stop)
    c = isnothing(d.f) ? 1. : d.f(t)
    if c < 0
        error("the function of a term holding dissipators is $c at time $t, which gives " *
              "their jumps a negative rate, which no trajectory draws: use State{Mixed}")
    end
    return c * real(expect(state, d.rate))
end

# the rate of all the jumps of `e` on `state` at `time`, rounding below zero taken as zero
total_rate(e::JumpEvolver, state::AbstractState{Trajectory}, time::Number) =
    max(sum(d -> decay_rate(d, state, time, e.stop), e.decays), 0.0)

# One step of a trajectory: the time of the next jump drawn from the probability of no jump up to
# it, as QuTiP's mcsolve and QuantumOptics.jl's mcwf do, the jump drawn by its rate at that time
# and the rest of the step drawn again. That probability is taken from the rates at both ends of
# the step, linear in between, not from the norm, which truncation lowers. The trajectory is
# brought back to its trace or, when the operators do not keep it or the norm is exact, see
# `exact_norm`, weighed by the norm left over the probability drawn with, which corrects the rate
# taken linear where it is not zero.
# An algorithm of one's own calling evolve_step of Tdvp or ApproxW with a JumpEvolver draws the
# jumps with them
function evolve_step(algo::Union{Tdvp, ApproxW}, e::JumpEvolver, state::AbstractState{Trajectory},
                     time::Number, dt::Number, step::Int, limits::Limits)
    if isempty(e.decays)
        return followed(evolve_step(algo, e.pre, state, time, dt, step, limits), state, 1)
    elseif !(isreal(dt) && real(dt) > 0)
        error("a trajectory draws its jumps forward in real time, and the step is $dt")
    end
    return jump_step(algo, e, state, time, real(dt), step, limits)
end

# operators keeping the trace leave no trajectory of norm zero but by underflow
no_norm_error(dt, time) =
    error("a step of $dt from time $time leaves a trajectory of norm zero by underflow, its " *
          "jumps too frequent for it: take smaller steps")

# the step of `evolve_step` for a JumpEvolver, or its rest after a jump
function jump_step(algo::Union{Tdvp, ApproxW}, e::JumpEvolver, state::AbstractState{Trajectory},
                   time::Number, dt::Real, step::Int, limits::Limits)
    if is_null(state)
        return state
    end
    rng = something(e.rng, Random.default_rng())
    brought = e.keeps && !e.exact
    a0 = total_rate(e, state, time)
    stepped = evolve_step(algo, e.pre, state, time, dt, step, limits)
    # a step losing all the norm, by the operators or by underflow, gives no rate at its end
    a1 = is_null(stepped) ? a0 : total_rate(e, stepped, time + dt)
    # no jump up to s with the probability exp(-a0 s - c s²)
    c = (a1 - a0) / (2dt)
    lq = -dt * (a0 + c * dt)
    lu = log(rand(rng))
    if lu < lq
        if is_null(stepped)
            return brought ? no_norm_error(dt, time) : stepped
        end
        return brought ? followed(stepped, state, 1) : stepped / exp(lq / 2)
    end
    s = clamp(-2lu / (a0 + sqrt(max(a0^2 - 4c * lu, 0.0))), 0.0, dt)
    before = evolve_step(algo, e.pre, state, time, s, step, limits)
    if is_null(before)
        return brought ? no_norm_error(dt, time) : before
    end
    rates = [ decay_rate(d, before, time + s, e.stop) for d in e.decays ]
    jumped = if sum(rates) ≤ 0
        # the linear rate drew a jump where there is none: weight zero when losing trace
        brought ? followed(before, state, 1) : 0 * before
    else
        d = e.decays[draw_index(rates, rng)]
        i = draw_index([ real(expect(before, dag(j) * j)) for j in d.jumps ], rng)
        x = ket_applied(d.jumps[i], before; limits)
        if !isnothing(e.jumps)
            dissipator, kraus = d.labels[i]
            push!(e.jumps, (; time = time + s, dissipator, kraus))
        end
        brought ? followed(x, state, 1) :
                  followed(x, before, sum(rates) / ((a0 + 2c * s) * exp(-s * (a0 + c * s))))
    end
    return s < dt ? jump_step(algo, e, jumped, time + s, dt - s, step, limits) : jumped
end

"""
    evolution_steps(what, span, nsteps, step, step_name)

the number of steps of `what` over `span`, given as `nsteps`, at least one, or by its step, the
keyword `step_name`: the nearest whole number of steps covering `span`, none for a span shorter
than half a step. Exactly one of the two is given, the other being `nothing`.
"""
function evolution_steps(what::String, span::Number, nsteps, step, step_name::String)
    if isnothing(nsteps) && isnothing(step)
        error("$what needs its number of steps, nsteps = 100 for instance, or its step, " *
              "$step_name = 0.01")
    elseif !isnothing(nsteps) && !isnothing(step)
        error("$what takes nsteps or $step_name, not both")
    elseif !isnothing(step) && iszero(step)
        error("the $step_name of $what cannot be zero")
    elseif !isnothing(nsteps) && nsteps < 1
        error("$what takes at least one step, and nsteps is $nsteps")
    end
    return isnothing(nsteps) ? round(Int, abs(span / step)) : nsteps
end

"""
    evolve(evolver, duration, state; nsteps | time_step, algo = Tdvp(), limits = Limits(),
           mpo_limits = Limits(), time = 0, step = 1, rng, jumps = nothing)

the state evolved for `duration` under `evolver`, in `nsteps` steps of `duration / nsteps`, or
in steps of `time_step` adjusted to the nearest one dividing `duration` into a whole number of
them (none when `duration` is shorter than half a step), one of the two being given, by
the method `algo`, see `EvolutionAlgorithm`. The evolver is `-im * H` for a Hamiltonian `H`,
plus dissipators for a mixed state or a trajectory, or its `PreMPO(state, evolver)`, prepared
once for several calls and keeping the MPOs they build. One that depends on time is a vector of terms, each an
operator or an operator with its real function of time, `[A, B => g]` for ``A + g(t) B``,
``t`` starting at `time`.
`step` is the number of the first step, for the periods of the methods and for `limits`, see
`Limits`, which may give one value per step, step `k` taking the `k`-th. Its default cutoff is
`1e-26` for tdvp and `apply_algo = "naive"` or `"zipup"`, `eps()` for `"densitymatrix"`.
`mpo_limits` truncates the MPO of the evolver, see `PreMPO` (default none).
`rng` is the random number generator a trajectory draws the jumps of the dissipators with, the
global one by default, and `jumps` a vector each jump is pushed to as
`(time, dissipator, kraus)`: its time, the dissipator as written and placed, as `Sm(1)` or
`Dephase(Level)(3)`, and the index of the operator drawn among the Kraus operators of a channel,
`TensorMixedStates.kraus_operators(system, dissipator)`, 1 for an operator.

# Examples

    evolve(-im * H, 1., state; nsteps = 10, limits = Limits(cutoff = 1e-12, maxdim = 50))
    evolve([-im * H0, -im * H1 => t -> cos(t)], 1., state; time_step = 0.05,
           algo = ApproxW(order = 4))
    evolve(-im * H + Dissipator(Sm)(1), 1., trajectory; nsteps = 10, jumps = list)
"""
function evolve(evolver, duration::Number, state::AbstractState; algo::EvolutionAlgorithm = Tdvp(),
                nsteps::Union{Nothing, Int} = nothing, time_step::Union{Nothing, Number} = nothing,
                limits::Limits = Limits(), mpo_limits::Limits = Limits(), time::Number = 0.,
                step::Int = 1, rng::AbstractRNG = Random.default_rng(),
                jumps::Union{Nothing, AbstractVector} = nothing)
    nsteps = evolution_steps("an evolution", duration, nsteps, time_step, "time_step")
    if nsteps == 0
        return state
    end
    pre = with_rng(prepared_evolver(state, evolver, mpo_limits; duration), rng, jumps;
                   stop = time + duration)
    dt = duration / nsteps
    for k in 0:nsteps - 1
        state = evolve_step(algo, pre, state, time + k * dt, dt, step + k,
                            step_limits(limits, step + k))
    end
    return state
end
