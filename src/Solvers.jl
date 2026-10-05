# The solvers: tdvp and approx_W for time evolution, dmrg for ground states, steady_state for
# steady states and thermal_state for thermal states.

export Krylov, tdvp, dmrg, approx_W, steady_state, thermal_state

"""
    Krylov(; dim = nothing, maxiter = nothing, tol = nothing)

the parameters of the Krylov method of each local step: the exponentiation of `tdvp` and
`thermal_state`, the eigenvector search of `dmrg` and `steady_state`. A field left to
`nothing` keeps the default of KrylovKit.

- `dim`: the largest dimension of a Krylov space (default 30 for `tdvp`, 3 for `dmrg`, 8 for
  `steady_state`)
- `maxiter`: the number of Krylov spaces built in turn (default 100 for `tdvp`, 1 for `dmrg`,
  3 for `steady_state`)
- `tol`: the tolerance (default `1e-12` per unit of time for `tdvp`, `1e-14` for `dmrg`); for
  `tdvp`, it is what sets the number of vectors, `dim` only bounding it

# Examples

    Krylov(tol = 1e-10)            # tdvp: fewer vectors, a lower precision
    Krylov(dim = 8, maxiter = 3)   # dmrg: a more accurate local step
"""
@kwdef struct Krylov
    dim::Union{Nothing, Int} = nothing
    maxiter::Union{Nothing, Int} = nothing
    tol::Union{Nothing, Float64} = nothing
end

show(io::IO, k::Krylov) =
    print(io, "Krylov(", join((string(f, " = ", repr(getfield(k, f))) for f in fieldnames(Krylov)
                               if !isnothing(getfield(k, f))), ", "), ")")

"""
    krylov_kwargs(k, prefix = "")

the fields of `k` that are given, as keyword arguments of KrylovKit prefixed by `prefix`
"""
krylov_kwargs(k::Krylov, prefix = "") =
    NamedTuple(Symbol(prefix, f == :dim ? :krylovdim : f) => getfield(k, f)
               for f in fieldnames(Krylov) if !isnothing(getfield(k, f)))

"""
    check_nsteps(nsteps)

refuse an evolution in less than one step
"""
function check_nsteps(nsteps)
    if nsteps < 1
        error("an evolution takes at least one step, and nsteps is $nsteps")
    end
end

"""
    check_pre_system(pre, state)

refuse a `PreMPO` prepared on another system than that of `state`
"""
function check_pre_system(pre::PreMPO, state::State)
    if pre.system !== state.system
        error("the PreMPO was prepared on another System than that of the state: prepare it " *
              "with the state, PreMPO(state, op)")
    end
end

"""
    tdvp_step(mpo, dt, st, state, sweep, expand_period, hermitianize_period, limits,
              updater_kwargs)

the MPS `st` after step `sweep` of tdvp, of `dt` under `mpo`, expanded before the step and
hermitianized after it when they are due
"""
function tdvp_step(mpo, dt, st, state, sweep, expand_period, hermitianize_period, limits,
                   updater_kwargs)
    # before the step: from a state of small bond dimension, a product state, an expansion
    # after the first step comes too late, which leaves an error of order dt
    if sweep_due(expand_period, sweep - 1)
        st = expand(st, mpo; alg="global_krylov")
    end
    st = tdvp(mpo, dt, st; nsweeps = 1, limits.cutoff, limits.maxdim, limits.mindim,
              updater_kwargs)
    if sweep_due(hermitianize_period, sweep)
        st = hermitianize(State(state, st); limits).state
    end
    return st
end

"""
    tdvp(evolver, t, ::State; options...)
    tdvp(evolver, t, ::Simulation; options...)

evolve a state, or a simulation, for a time `t` by tdvp, in `nsteps` steps of `t / nsteps`.
`evolver` is `-im * H` for a Hamiltonian `H`, plus dissipators for a mixed state, or its
`PreMPO(state, evolver)`, to prepare it once for several calls. A simulation comes back with
its time advanced by `t`.

# Options

- `nsteps`: the number of steps (default 1)
- `first_step`: the step to start from (default 1), to continue an evolution: `t`, `nsteps`
  and `time_start` stay those of the whole evolution
- `time_start`: the simulation time the evolution starts from (default 0, and the time of
  the simulation for a `Simulation`)
- `coefs`: for a vector of evolvers, the real functions of time they are multiplied by, taken
  at the middle of each step, see [Time dependent evolvers](@ref)
- `expand_period`: enlarge the bond dimension by a global Krylov expansion before the first
  step and then every `expand_period` steps (default 0, never)
- `hermitianize_period`: make a mixed state hermitian every `hermitianize_period` steps
  (default 0, never)
- `limits`: constraints on the state, see `Limits`, which may give one value per step
  (default `Limits()`)
- `observer!`: an observer, see `TdvpObserver`
- `krylov`: the Krylov exponentiation of each local step, see `Krylov`

# Examples

    tdvp(-im * H, 1., state; nsteps = 10, limits = Limits(cutoff = 1e-10, maxdim = 50))
"""
function tdvp(pre::PreMPO{R}, t::Number, state::State{R};
    observer! = NoObserver(), coefs=nothing, expand_period = 0, hermitianize_period = 0,
    nsteps = 1, first_step = 1, time_start = zero(t), limits::Limits=Limits(),
    krylov::Krylov = Krylov()) where {R <: PM}
    check_pre_system(pre, state)
    check_nsteps(nsteps)
    time_dep = !isnothing(coefs)
    st = state.state
    dt = t / nsteps
    # eager: KrylovKit tests convergence after every vector rather than once its space is
    # full, which the short local steps of tdvp reach with far fewer products
    updater_kwargs = (; eager = true, krylov_kwargs(krylov)...)
    if !time_dep
        mpo = make_mpo(pre)
    end
    for sweep in first_step:nsteps
        current_time = time_start + sweep * dt
        if time_dep
            tf = current_time - dt / 2
            mpo = make_mpo(pre, map(f->f(tf), coefs))
        end
        st = tdvp_step(mpo, dt, st, state, sweep, expand_period, hermitianize_period,
                       sweep_limits(limits, sweep), updater_kwargs)
        # after the whole step, hermitianization included, so that a checkpoint holds the
        # state the measurements saw
        if sweep_done!(observer!; sweep, state = st, current_time, mpo)
            break
        end
    end
    return State(state, st)
end

tdvp(op, t::Number, state::State; kwargs...) =
    tdvp(PreMPO(state, op), t, state; kwargs...)

"""
    dmrg(hamiltonian, ::State; options...)
    dmrg(hamiltonian, ::Simulation; options...)

the ground state of a hermitian Hamiltonian by dmrg from the given state, returned as
`(energy, state)` or `(energy, simulation)`. The Hamiltonian may be given as its
`PreMPO(state, hamiltonian)` or its `make_mpo(state, hamiltonian)`, to build it once for
several searches. A mixed state is refused: search the ground
state of the pure state and `mix` it. The steady state of a Lindbladian is given by
`steady_state`.

# Options

- `nsweeps`: the number of sweeps of the whole search (default 1)
- `first_sweep`: the sweep to start from (default 1), to continue a search
- `limits`: constraints on the state, see `Limits`, which may give one value per sweep
  (default `Limits()`)
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `observer!`: an observer, see `DmrgObserver`
- `krylov`: the Krylov search of each local step, see `Krylov`

# Examples

    energy, state = dmrg(H, state; nsweeps = 10, limits = Limits(maxdim = [10, 20, 50]))
"""
function dmrg(mpo::MPO, state::State; nsweeps = 1, first_sweep = 1, observer! = NoObserver(),
              limits::Limits = Limits(), noise = 0., krylov::Krylov = Krylov())
    # ITensorMPS counts its sweeps from 1: a resumed search asks for those it has left, with
    # the tail of its schedules
    done = first_sweep - 1
    # ITensorMPS, asked for no sweep, gives an energy of 0
    if done ≥ nsweeps
        st = state.state
        return (real(inner(st', mpo, st) / inner(st, st)), state)
    end
    lim = resume_schedule(limits, done)
    e, st = dmrg(mpo, state.state; outputlevel = 0, nsweeps = nsweeps - done,
                 observer = observer!, lim.cutoff, lim.maxdim, lim.mindim,
                 noise = resume_schedule(noise, done), krylov_kwargs(krylov, "eigsolve_")...)
    return (e, State(state, st))
end

function dmrg(pre::PreMPO, state::State; kwargs...)
    check_pre_system(pre, state)
    return dmrg(make_mpo(pre), state; kwargs...)
end

function dmrg(op, state::State; kwargs...)
    if state isa State{Mixed} && op isa IndexedOp{Pure}
        error("dmrg finds the ground state of a pure state: search it pure and mix it")
    end
    return dmrg(make_mpo(state, op), state; kwargs...)
end

"""
    w_approx_coefs

for each order from 1 to 4, the fractions of the time step whose W approximations, applied in
turn, make up the approximation of that order
"""
const w_approx_coefs = Vector{ComplexF64}[
    [
        1.
    ],
    [
        0.5 + 0.5im,
        0.5 - 0.5im
    ],
    [
        0.10566243270259355 - 0.39433756729740643im,
        0.39433756729740643 + 0.10566243270259355im,
        0.39433756729740643 - 0.10566243270259355im,
        0.10566243270259355 + 0.39433756729740643im
    ],
    [
        0.2588533986109182 + 0.0447561340111419im,
        -0.03154685814880379 + 0.24911905427556322im,
        0.1908290521106672 - 0.23185374923210605im,
        0.16372881485443674,
        0.1908290521106672 + 0.23185374923210605im,
        -0.03154685814880379 - 0.24911905427556322im,
        0.2588533986109182 - 0.0447561340111419im,
    ]
]

"""
    check_w_approx(order, w)

refuse an `order` or a `w` that `approx_W` does not offer
"""
function check_w_approx(order::Int, w::Int)
    if order < 1 || order > length(w_approx_coefs)
        error("W approximation of order $order is not implemented")
    end
    if w ∉ (1, 2)
        error("W approximation is only defined for w=1 or 2 (not $w)")
    end
end

"""
    make_approx_W(pre, t; order, w, coefs = [1.])

the MPOs, to apply in turn, of the approximation of `order` of a step `t` by WI (`w = 1`) or
WII (`w = 2`), `coefs` being the values of the time functions
"""
function make_approx_W(pre::PreMPO, t::Number; order::Int, w::Int, coefs = [1.])
    check_w_approx(order, w)
    make = w == 1 ? make_approx_W1 : make_approx_W2
    return map(c -> make(pre, t * c, coefs), w_approx_coefs[order])
end

"""
    step_coefs(coefs, t, dt, order)

the values of the time functions `coefs` for each exponential of a step of `dt` from `t`, in
the order they are applied: at the middle of the step up to order 2, and from order 3 the
commutator-free Magnus integrator of order 4 of Blanes and Moan, two exponentials combining
the values at the two Gauss points of the step
"""
function step_coefs(coefs, t::Number, dt::Number, order::Int)
    if order ≤ 2
        return [ map(f -> f(t + dt / 2), coefs) ]
    end
    v1, v2 = (map(f -> f(t + (1/2 + s * √3/6) * dt), coefs) for s in (-1, 1))
    w1, w2 = (3 + 2√3) / 12, (3 - 2√3) / 12
    return [ w1 * v1 + w2 * v2, w2 * v1 + w1 * v2 ]
end

"""
    approx_W(evolver, t, ::State; order, options...)
    approx_W(evolver, t, ::Simulation; order, options...)

evolve a state, or a simulation, for a time `t` in `nsteps` steps of `t / nsteps`, each the
exponential of the evolver approximated at the given `order` by WI or WII approximations.
`evolver` is as for `tdvp`, and a simulation comes back with its time advanced by `t`.

# Options

- `order`: the order of the approximation, from 1 to 4, required
- `w`: 1 or 2 for WI or WII (default 2)
- `nsteps`: the number of steps (default 1)
- `first_step`: the step to start from (default 1), to continue an evolution: `t`, `nsteps`
  and `time_start` stay those of the whole evolution
- `time_start`: the simulation time the evolution starts from (default 0, and the time of
  the simulation for a `Simulation`)
- `coefs`: for a vector of evolvers, the real functions of time they are multiplied by, taken
  at the middle of each step up to order 2 and at two points from order 3, which keeps the
  order of the approximation, see [Time dependent evolvers](@ref)
- `hermitianize_period`: make a mixed state hermitian every `hermitianize_period` steps
  (default 0, never)
- `limits`: constraints on the state, see `Limits`, which may give one value per step
  (default `Limits()`)
- `observer!`: an observer, see `ApproxWObserver`
- `apply_algo`: the algorithm of the product of the state by each MPO, `"densitymatrix"`
  (default) or `"naive"`

# Examples

    approx_W(-im * H, 1., state; order = 4, nsteps = 10)
"""
function approx_W(pre::PreMPO{R}, t::Number, state::State{R}; coefs = nothing,
    hermitianize_period::Int = 0, nsteps::Int = 1, first_step::Int = 1, order::Int, w::Int = 2,
    observer! = NoObserver(), time_start = zero(t), limits::Limits=Limits(),
    apply_algo::String = "densitymatrix") where {R <: PM}
    check_pre_system(pre, state)
    check_apply_algo(apply_algo)
    check_nsteps(nsteps)
    st = state.state
    dt = t / nsteps
    time_dep = !isnothing(coefs)
    if !time_dep
        mpos = make_approx_W(pre, dt; order, w)
    end
    for sweep in first_step:nsteps
        current_time = time_start + sweep * dt
        if time_dep
            mpos = reduce(vcat, [ make_approx_W(pre, dt; order, w, coefs = c)
                                  for c in step_coefs(coefs, current_time - dt, dt, order) ])
        end
        lim = sweep_limits(limits, sweep)
        for mpo in mpos
            st = apply(mpo, st; alg = apply_algo, lim.cutoff, lim.maxdim, lim.mindim)
        end
        if sweep_due(hermitianize_period, sweep)
            st = hermitianize(State(state, st); limits = lim).state
        end
        if sweep_done!(observer!; sweep, state = st, current_time, mpos)
            break
        end
    end
    return State(state, st)
end

approx_W(op, t::Number, state::State; kwargs...) =
    approx_W(PreMPO(state, op), t, state; kwargs...)

"""
    check_mpo_algo(mpo_algo)

refuse an algorithm of ``L^\\dagger L`` that `steady_state` does not offer
"""
function check_mpo_algo(mpo_algo::String)
    if mpo_algo ∉ ("naive", "zipup")
        error("mpo_algo is \"naive\" or \"zipup\", not $(repr(mpo_algo))")
    end
end

"""
    steady_state(lindbladian, ::State; options...)
    steady_state(lindbladian, ::Simulation; options...)

the steady state of a Lindbladian ``L``, `-im * H` plus dissipators, or its
`PreMPO(state, lindbladian)`, by dmrg on ``L^\\dagger L`` from the given mixed state, returned
with trace one as `(value, state)` or `(value, simulation)`, `value` being the energy dmrg
reaches, zero for a steady state.

A Lindbladian may have several steady states, one in each sector of a quantity it conserves
for instance, and dmrg returns any combination of them, in general not a density matrix. The
one the system reaches from a given state is the limit of its evolution, which `tdvp` gives;
when a site can conserve the quantity, declaring it `strong` keeps the search in the sector of
the starting state, see [Conserving a quantity](@ref). A warning is given when the state found
is not hermitian, its `HermiticityError` above `1e-6`: the search has not converged, or the
steady state is not unique, see [Checking the accuracy](@ref).

# Options

- `nsweeps`: the number of sweeps of the whole search (default 1)
- `first_sweep`: the sweep to start from (default 1), to continue a search
- `limits`: constraints on the state, see `Limits`, which may give one value per sweep
  (default `Limits()`)
- `mpo_limits`: the truncation of the MPO of ``L^\\dagger L`` (default `Limits()`)
- `mpo_algo`: the algorithm computing ``L^\\dagger L``, `"naive"` (default) or `"zipup"`
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `observer!`: an observer, see `DmrgObserver`
- `krylov`: the Krylov search of each local step, see `Krylov` (default
  `Krylov(dim = 8, maxiter = 3)`, which the spectrum of ``L^\\dagger L``, crowded near zero,
  needs)

# Examples

    value, rho = steady_state(-im * H + D, rho; nsweeps = 10,
                              limits = Limits(maxdim = [10, 20, 50]))
"""
steady_state(op::IndexedOp{Mixed}, state::State{Mixed}; kwargs...) =
    steady_state(PreMPO(state, op), state; kwargs...)

function steady_state(pre::PreMPO{Mixed}, state::State{Mixed};
    limits::Limits = Limits(), nsweeps::Int = 1, first_sweep::Int = 1,
    observer! = NoObserver(), mpo_limits::Limits = Limits(), mpo_algo::String = "naive",
    noise = 0., krylov::Krylov = Krylov(dim = 8, maxiter = 3))
    check_pre_system(pre, state)
    check_mpo_algo(mpo_algo)
    l = make_mpo(pre)
    # the naive algorithm truncates only when asked to, the others take no such option
    extra = mpo_algo == "naive" ? (; truncate = mpo_limits != Limits()) : (;)
    l2 = apply(replaceprime(dag(l)', 2=>0), l;
               mpo_limits.cutoff, mpo_limits.maxdim, mpo_limits.mindim, alg = mpo_algo, extra...)
    e, st = dmrg(l2, state; nsweeps, first_sweep, limits, noise, observer!, krylov)
    # the trace, a complex number, takes away the phase of the eigenvector, which dividing by
    # its real part, as `normalize` does, would leave
    ρ = st / trace(st)
    δ = 1 - hermiticity(ρ)
    if δ > 1e-6
        @warn "the steady state found has a HermiticityError of $(short(δ)): the search has " *
              "not converged, or the steady state is not unique, see steady_state"
    end
    return (e, ρ)
end

"""
    thermal_state(hamiltonian, beta, ::State; options...)
    thermal_state(hamiltonian, beta, ::Simulation; options...)

the mixed state ``\\rho`` taken to ``e^{-\\beta H/2} \\rho \\, e^{-\\beta H/2}`` by tdvp in
imaginary time, in `nsteps` steps of `beta / nsteps`, and normalized to trace one. It is
returned as `(log_trace, state)` or `(log_trace, simulation)`, the time of the simulation
unchanged, `log_trace` being ``\\log \\mathrm{tr}(e^{-\\beta H/2} \\rho \\, e^{-\\beta H/2})`` for
``\\rho`` of trace one.

From `"FullyMixed"`, it gives the thermal state ``e^{-\\beta H}/Z``, and `log_trace` is
``\\log Z - \\sum_i \\log d_i``, ``d_i`` being the dimension of site `i`. A state commuting
with ``H`` gives the thermal state restricted to what it describes: `fully_mixed(system, N => m)`
the canonical state of `m` particles, a product state of populations ``e^{\\beta\\mu n_i}``
the grand canonical state at chemical potential ``\\mu``. A pure state is refused.

# Options

- `nsteps`: the number of steps (default 1)
- `first_step`: the step to start from (default 1), to continue a computation from the state
  it had reached: `beta` and `nsteps` stay those of the whole computation
- `log_trace`: the `log_trace` the computation starts from (default 0), to continue it
- `expand_period`: enlarge the bond dimension by a global Krylov expansion before the first
  step and then every `expand_period` steps (default 0, never)
- `hermitianize_period`: make the state hermitian every `hermitianize_period` steps (default 0,
  never)
- `limits`: constraints on the state, see `Limits`, which may give one value per step
  (default `Limits()`)
- `observer!`: an observer, see `ThermalObserver`
- `krylov`: the Krylov exponentiation of each local step, see `Krylov`

# Examples

    log_trace, rho = thermal_state(H, 2.0, State{Mixed}(system, "FullyMixed"); nsteps = 20,
                                   limits = Limits(cutoff = 1e-12, maxdim = 64))
"""
function thermal_state(hamiltonian::IndexedOp{Pure}, beta::Real, state::State{Mixed};
    observer! = NoObserver(), nsteps::Int = 1, first_step::Int = 1, log_trace::Real = 0.,
    expand_period::Int = 0, hermitianize_period::Int = 0, limits::Limits = Limits(),
    krylov::Krylov = Krylov())
    check_nsteps(nsteps)
    dbeta = beta / nsteps
    # on a mixed state, -H/2 becomes its Evolver, ρ ↦ -(Hρ + ρH)/2, the generator of
    # e^{-βH/2} ρ e^{-βH/2}
    mpo = make_mpo(state, -hamiltonian / 2)
    updater_kwargs = (; eager = true, krylov_kwargs(krylov)...)
    # a continued computation is not normalized again, which would change its rounding
    st = first_step == 1 ? normalize(state).state : state.state
    for step in first_step:nsteps
        st = tdvp_step(mpo, dbeta, st, state, step, expand_period, hermitianize_period,
                       sweep_limits(limits, step), updater_kwargs)
        t = real(trace(State(state, st)))
        if !(t > 0)
            error("the trace of the state reached $t at step $step of thermal_state: take " *
                  "smaller steps or larger limits")
        end
        st /= t
        log_trace += log(t)
        if sweep_done!(observer!; sweep = step, state = st, beta = step * dbeta, log_trace)
            break
        end
    end
    return (log_trace, State(state, st))
end

thermal_state(::IndexedOp{Pure}, ::Real, ::State{Pure}; kwargs...) =
    error("a thermal state needs a mixed representation: mix the state first, see ToMixed")

tdvp(op, t::Number, sim::Simulation; kwargs...) =
    Simulation(sim, tdvp(op, t, sim.state; time_start = sim.time, kwargs...), sim.time + t)

function dmrg(op, sim::Simulation; kwargs...)
    e, st = dmrg(op, sim.state; kwargs...)
    return (e, Simulation(sim, st))
end

approx_W(op, t::Number, sim::Simulation; kwargs...) =
    Simulation(sim, approx_W(op, t, sim.state; time_start = sim.time, kwargs...), sim.time + t)

function steady_state(op, sim::Simulation; kwargs...)
    e, st = steady_state(op, sim.state; kwargs...)
    return (e, Simulation(sim, st))
end

function thermal_state(hamiltonian, beta, sim::Simulation; kwargs...)
    l, st = thermal_state(hamiltonian, beta, sim.state; kwargs...)
    return (l, Simulation(sim, st))
end
