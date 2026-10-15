# The evolution of a state in time: evolve, a loop over the steps that evolve_step makes for each
# method, Tdvp and ApproxW, with the evolver prepared once and its MPOs kept, see cached_mpos.

export evolve

"""
    check_pre(pre, state)

refuse a `PreMPO` prepared on another system than that of `state`
"""
function check_pre(pre::PreMPO, state::State)
    if pre.system ≠ state.system
        error("the PreMPO was prepared on another system than that of the state: prepare it " *
              "with the state, PreMPO(state, op)")
    end
    return pre
end

"""
    prepared_evolver(state, evolver)

the `PreMPO` of the evolver of `state` and its time functions, `nothing` for an evolver that
does not depend on time: an operator, its `PreMPO`, or a vector of operators, or its `PreMPO`,
paired with as many functions of time, `[A, B] => [f, g]`
"""
prepared_evolver(state::State, a::IndexedOp) = (PreMPO(state, a), nothing)
prepared_evolver(state::State{R}, pre::PreMPO{R}) where {R <: PM} = (check_pre(pre, state), nothing)

function prepared_evolver(state::State,
                          (ops, fs)::Pair{<:Union{AbstractVector, PreMPO}, <:AbstractVector})
    pre = ops isa PreMPO ? first(prepared_evolver(state, ops)) : PreMPO(state, [ o for o in ops ])
    if length(fs) ≠ pre.nterms
        error("an evolver of $(pre.nterms) terms takes as many time functions, " *
              "got $(length(fs))")
    end
    return (pre, collect(fs))
end

prepared_evolver(::State, ::MPO) =
    error("the solvers take an operator or its PreMPO, which keeps the MPOs they build, not an MPO")

prepared_evolver(::State, ::PreMPO) =
    error("the PreMPO was prepared on a state of the other representation")

prepared_evolver(::State, x) = error("$(repr(x)) is not an evolver")

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
    expanded(state, mpo)

the state enlarged by a global Krylov expansion under `mpo`, its Krylov vectors computed as `apply`
computes a product, "naive" and truncated at rounding alone: truncated to one more than the bond
dimension of the state, the default of ITensorMPS, they keep nothing of some charged states, as
FullyMixed. The other options are those of ITensorMPS. The state is then written again in the
bonds the expansion enlarged, see `restored`. An expansion only enlarges the basis of the step,
which goes on without it, with a warning, when ITensorMPS fails on charged indices, as on a state
that `mpo` maps onto itself, and when the state written again is off by more than the precision
ITensorMPS gives the expansion, the square root of that of the numbers
"""
function expanded(state::State, mpo::MPO)
    st = state.state
    e = try
        expand(st, mpo; alg = "global_krylov",
               apply_kwargs = (; alg = "naive", truncation(Limits())...))
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

"""
    evolve_step(algo, pre, coefs, state, time, dt, step, limits)

the state after one step of `dt` from `time` by the method `algo`, under the evolver `pre`,
its time functions `coefs`, `nothing` when it does not depend on time, `step` being the number
of the step, from 1, and `limits` its limits, see `step_limits`. An `EvolutionAlgorithm` of
one's own adds a method.
"""
function evolve_step(algo::Tdvp, pre::PreMPO, coefs, state::State, time::Number, dt::Number,
                     step::Int, limits::Limits)
    # an evolution leaves the null state null, which ITensors fails on
    if vanished(state.state)
        return state
    end
    mpo = isnothing(coefs) ? cached_mpos(() -> make_mpo(pre), pre, :mpo, nothing) :
                             make_mpo(pre, [ f(time + dt / 2) for f in coefs ])
    tr = truncation(limits)
    # before the step: from a state of small bond dimension, a product state, an expansion after
    # the step comes too late, which leaves an error of order dt; a single site has no bond
    small = maxlinkdim(state) < min(5, tr.maxdim)
    if length(state) > 1 && (small || due(algo.expand_period, step - 1))
        state = expanded(state, mpo)
    end
    # eager: KrylovKit tests convergence after every vector rather than once its space is full,
    # which the short local steps of tdvp reach with far fewer products
    st = tdvp(mpo, dt, state.state; nsweeps = 1, tr...,
              updater_kwargs = (; eager = true, krylov_kwargs(algo.krylov)...))
    state = State(state, st)
    if due(algo.hermitianize_period, step)
        state = hermitianize(state; limits)
    end
    return state
end

function evolve_step(algo::ApproxW, pre::PreMPO, coefs, state::State, time::Number, dt::Number,
                     step::Int, limits::Limits)
    if vanished(state.state)
        return state
    end
    mpos = if isnothing(coefs)
        cached_mpos(() -> approx_W_mpos(pre, dt, algo), pre, :approx_W, (algo.order, algo.w, dt))
    else
        reduce(vcat, [ approx_W_mpos(pre, dt, algo, c)
                       for c in step_coefs(coefs, time, dt, algo.order) ])
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

"""
    evolve(evolver, duration, state; algo = Tdvp(), nsteps = 1, limits = Limits(), time = 0,
           step = 1)

the state evolved for `duration` under `evolver`, in `nsteps` steps of `duration / nsteps` by
the method `algo`, see `EvolutionAlgorithm`. The evolver is `-im * H` for a Hamiltonian `H`,
plus dissipators for a mixed state, or its `PreMPO(state, evolver)`, prepared once for several
calls and keeping the MPOs they build; one that depends on time pairs a vector of them with as
many real functions of time, `[A, B] => [f, g]`, read at `time` plus the time of each step.
`step` is the number of the first step, for the periods of the methods and for `limits`, see
`Limits`, which may give one value per step, step `k` taking the `k`-th. Its default cutoff is
`1e-26` for tdvp and `apply_algo = "naive"`, `eps()` for `"densitymatrix"`.

# Examples

    evolve(-im * H, 1., state; nsteps = 10, limits = Limits(cutoff = 1e-12, maxdim = 50))
    evolve([-im * H0, -im * H1] => [t -> 1, t -> cos(t)], 1., state; algo = ApproxW(order = 4))
"""
function evolve(evolver, duration::Number, state::State; algo::EvolutionAlgorithm = Tdvp(),
                nsteps::Int = 1, limits::Limits = Limits(), time::Number = 0., step::Int = 1)
    if nsteps < 1
        error("an evolution takes at least one step, and nsteps is $nsteps")
    end
    pre, coefs = prepared_evolver(state, evolver)
    dt = duration / nsteps
    for k in 0:nsteps - 1
        state = evolve_step(algo, pre, coefs, state, time + k * dt, dt, step + k,
                            step_limits(limits, step + k))
    end
    return state
end
