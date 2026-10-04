# The algorithms on a state or a simulation: tdvp and approx_W for time evolution, dmrg for
# ground states, and steady_state for the steady state of an open system.

export Krylov, tdvp, dmrg, approx_W, steady_state, thermal_state

"""
    Krylov(; dim = nothing, maxiter = nothing, tol = nothing)

the parameters of the Krylov method that solves each local step: `KrylovKit.exponentiate`,
computing the exponential of `tdvp`, and `KrylovKit.eigsolve`, the lowest eigenvector of
`dmrg` and `steady_state`. A field left to `nothing` keeps the default of the method.

# Fields

- `dim`: the largest dimension of a Krylov space, the `krylovdim` of KrylovKit (default 30 for
  `tdvp`, 3 for `dmrg`, 8 for `steady_state`)
- `maxiter`: the number of Krylov spaces built one after the other (default 100 for `tdvp`, 1
  for `dmrg`, 3 for `steady_state`): `tdvp` covers in several parts a step that one space does
  not cover at the tolerance, `dmrg` restarts from its best vectors
- `tol`: the tolerance (default `1e-12` per unit of time for `tdvp`, `1e-14` for `dmrg`)

`tdvp` tests its convergence after every vector, so that `dim` bounds the number of vectors
rather than fixing it, and `tol` is what sets it.

# Examples

    Krylov(tol = 1e-10)            # tdvp: fewer vectors, a lower precision
    Krylov(dim = 8, maxiter = 3)   # dmrg: a more accurate local step
"""
@kwdef struct Krylov
    dim::Union{Nothing, Int} = nothing
    maxiter::Union{Nothing, Int} = nothing
    tol::Union{Nothing, Float64} = nothing
end

# printed as the call that builds it, the fields left to the default of the method omitted
show(io::IO, k::Krylov) =
    print(io, "Krylov(", join((string(f, " = ", repr(getfield(k, f))) for f in fieldnames(Krylov)
                               if !isnothing(getfield(k, f))), ", "), ")")

"""
    krylov_kwargs(::Krylov, prefix = "")

the fields of a `Krylov` that are given, as the keyword arguments of KrylovKit, `dim` being
its `krylovdim`, with `prefix` in front, those left to `nothing` omitted so that the method
keeps its default.
"""
krylov_kwargs(k::Krylov, prefix = "") =
    NamedTuple(Symbol(prefix, f == :dim ? :krylovdim : f) => getfield(k, f)
               for f in fieldnames(Krylov) if !isnothing(getfield(k, f)))

"""
    check_nsweeps(nsweeps, name = "nsweeps")

refuse an evolution in less than one step, its number of steps being called `name`: it made no
step, while a simulation took the time on, and `approx_W` divided its duration by zero
"""
function check_nsweeps(nsweeps, name = "nsweeps")
    if nsweeps < 1
        error("an evolution takes at least one step, and $name is $nsweeps")
    end
end

"""
    check_pre_system(pre, state)

refuse a `PreMPO` prepared on another system than that of the state, whose indices differ:
ITensorMPS refused the product with a message naming neither
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

the MPS `st` of a state of the representation of `state` after step `sweep` of tdvp, a step
of `dt` under `mpo`: expanded before the step and hermitianized after it when they are due,
see `tdvp`
"""
function tdvp_step(mpo, dt, st, state, sweep, expand_period, hermitianize_period, limits,
                   updater_kwargs)
    # before the step rather than after it, on the schedule shifted by one: the first
    # step from a product state, of bond dimension one, left the tangent space and kept
    # an error of order dt, the expansion after it coming too late
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

evolve a state, or a simulation, for a time `t` with the tdvp algorithm, in `nsweeps` steps of
`t / nsweeps`. `evolver` is `-im * H` for a hamiltonian `H`, plus dissipators for a mixed
state, or its `PreMPO`, prepared once for the calls that evolve under it, one step at a time
for instance. A simulation comes back with its time advanced by `t`.

# Options

- `nsweeps`: the number of steps (default 1)
- `first_sweep`: the step to start from (default 1), to continue an evolution left unfinished:
  `t`, `nsweeps` and `time_start` are still those of the whole evolution, and a simulation
  is given at the time that evolution started
- `time_start`: the simulation time the evolution starts from (default 0, and the time of
  the simulation for a `Simulation`)
- `coefs`: for a vector of evolvers, the real functions of time they are multiplied by, taken
  at the middle of each step
- `expand_period`: enlarge the bond dimension of the state by a global Krylov expansion before the
  first step and then every `expand_period` steps (default 0, never)
- `hermitianize_period`: make a mixed state hermitian every `hermitianize_period` steps (default 0,
  never)
- `limits`: constraints on the state, see `Limits`, which may give one value per step
  (default `Limits()`)
- `observer!`: an observer, see `TdvpObserver`
- `krylov`: the parameters of the Krylov exponentiation of each local step, see `Krylov`
  (default `Krylov()`, those of `KrylovKit.exponentiate`)

# Examples

    tdvp(-im * H, 1., state; nsweeps = 10, limits = Limits(cutoff = 1e-10, maxdim = 50))
"""
function tdvp(pre::PreMPO{R}, t::Number, state::State{R};
    observer! = NoObserver(), coefs=nothing, expand_period = 0, hermitianize_period = 0,
    nsweeps = 1, first_sweep = 1, time_start = zero(t), limits::Limits=Limits(),
    krylov::Krylov = Krylov()) where {R <: PM}
    check_pre_system(pre, state)
    check_nsweeps(nsweeps)
    time_dep = !isnothing(coefs)
    st = state.state
    dt = t / nsweeps
    # KrylovKit builds its whole Krylov space, `dim` vectors, before it tests
    # convergence: tested after every vector, as `eager` does, it stops at the same
    # tolerance, which the short local steps of tdvp reach with far fewer products
    updater_kwargs = (; eager = true, krylov_kwargs(krylov)...)
    if !time_dep
        mpo = make_mpo(pre)
    end
    for sweep in first_sweep:nsweeps
        current_time = time_start + sweep * dt
        if time_dep
            tf = current_time - dt / 2
            mpo = make_mpo(pre, map(f->f(tf), coefs))
        end
        st = tdvp_step(mpo, dt, st, state, sweep, expand_period, hermitianize_period,
                       sweep_limits(limits, sweep), updater_kwargs)
        # once the whole sweep is done: the sweep is committed right after its measurements,
        # and a checkpoint cannot fall between them
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

the ground state of a hamiltonian by dmrg, starting from the given state, returned as
`(energy, state)`, or `(energy, simulation)`. The hamiltonian may be given by its `PreMPO` or its
MPO, prepared once for several searches. A hamiltonian is refused on a mixed state, where
the lowest eigenvector of the superoperator it gives is neither the ground state nor a
density matrix: search the ground state of the pure state, then `mix` it.

# Options

- `nsweeps`: the last sweep to do, that is the number of sweeps of the whole run (default 1)
- `first_sweep`: the sweep to start from (default 1), to continue a search left unfinished;
  past `nsweeps`, no sweep is done and the energy is that of the state given
- `limits`: constraints on the state, see `Limits`, which may give one value per sweep
  (default `Limits()`)
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `observer!`: an observer, see `DmrgObserver`
- `krylov`: the parameters of the Krylov search of the lowest eigenvector at each local step,
  see `Krylov` (default `Krylov()`, those of `ITensorMPS.dmrg`)

# Examples

    energy, state = dmrg(H, state; nsweeps = 10, limits = Limits(maxdim = [10, 20, 50]))
"""
function dmrg(mpo::MPO, state::State; nsweeps = 1, first_sweep = 1, observer! = NoObserver(),
              limits::Limits = Limits(), noise = 0., krylov::Krylov = Krylov())
    # ITensorMPS counts its sweeps from 1 and offers no way to start elsewhere, so a run
    # resuming at `first_sweep` asks for the sweeps it has left and is handed the tail of
    # its per sweep schedules. `tdvp` and `approx_W` drive their own loop and keep the
    # sweep numbers of the run instead, which is why only this one has to adapt.
    done = first_sweep - 1
    # no sweep left, for a search resumed after its last one: the energy is that of the state
    # given, where ITensorMPS, asked for no sweep, gave 0
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
    # a hamiltonian on a mixed state becomes the superoperator ρ ↦ Hρ + ρH, whose lowest
    # eigenvector is neither the ground state nor a density matrix
    if state isa State{Mixed} && op isa IndexedOp{Pure}
        error("dmrg finds the ground state of a pure state: search it pure and mix it")
    end
    return dmrg(make_mpo(state, op), state; kwargs...)
end

"""
    w_approx_coefs

for each order from 1 to 4, the fractions of the time step whose W approximations, applied one
after the other, make up the approximation of that order.
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

refuse an approximation by WI and WII of an order or of a `w` not offered
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

the MPOs of the approximation of the given `order` of a step `t`, to apply one after the
other, built from WI (`w = 1`) or WII (`w = 2`) approximations; `coefs` are the real values of
the coefficients of a time dependent evolver.
"""
function make_approx_W(pre::PreMPO, t::Number; order::Int, w::Int, coefs = [1.])
    check_w_approx(order, w)
    make = w == 1 ? make_approx_W1 : make_approx_W2
    return map(c -> make(pre, t * c, coefs), w_approx_coefs[order])
end

"""
    approx_W(evolver, t, ::State; order, options...)
    approx_W(evolver, t, ::Simulation; order, options...)

evolve a state, or a simulation, for a time `t` in `nsweeps` steps of `t / nsweeps`, each
approximating the exponential of the evolver at the given `order` with WI or WII
approximations. `evolver` is as for `tdvp`, and a simulation comes back with its time
advanced by `t`.

# Options

- `order`: the order of the approximation, from 1 to 4, required
- `w`: 1 or 2 for WI or WII (default 2)
- `nsweeps`: the number of steps (default 1)
- `first_sweep`: the step to start from (default 1), to continue an evolution left unfinished:
  `t`, `nsweeps` and `time_start` are still those of the whole evolution, and a simulation
  is given at the time that evolution started
- `time_start`: the simulation time the evolution starts from (default 0, and the time of
  the simulation for a `Simulation`)
- `coefs`: for a vector of evolvers, the real functions of time they are multiplied by, taken
  at the middle of each step
- `hermitianize_period`: make a mixed state hermitian every `hermitianize_period` steps (default 0,
  never)
- `limits`: constraints on the state, see `Limits`, which may give one value per step
  (default `Limits()`)
- `observer!`: an observer, see `ApproxWObserver`
- `apply_algo`: the algorithm of the product of the state by each MPO, as `ITensorMPS.apply`
  takes it: `"densitymatrix"` (default) or `"naive"`. `"fit"` is not offered, since it needs a
  number of sweeps of its own, nor `"zipup"`, which versions of ITensorMPS before 0.3.45 lack

# Examples

    approx_W(-im * H, 1., state; order = 4, nsweeps = 10)
"""
function approx_W(pre::PreMPO{R}, t::Number, state::State{R}; coefs = nothing,
    hermitianize_period::Int = 0, nsweeps::Int = 1, first_sweep::Int = 1, order::Int, w::Int = 2, observer! = NoObserver(),
    time_start = zero(t), limits::Limits=Limits(), apply_algo::String = "densitymatrix") where {R <: PM}
    check_pre_system(pre, state)
    check_apply_algo(apply_algo)
    check_nsweeps(nsweeps)
    st = state.state
    dt = t / nsweeps
    time_dep = !isnothing(coefs)
    if !time_dep
        mpos = make_approx_W(pre, dt; order, w)
    end
    for sweep in first_sweep:nsweeps
        current_time = time_start + sweep * dt
        if time_dep
            tf = current_time - dt / 2
            mpos = make_approx_W(pre, dt; order, w, coefs = map(f->f(tf), coefs))
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

refuse an algorithm of the product of two MPOs, computing ``L^\\dagger L``, that is not one of
those `steady_state` offers, `"naive"` and `"zipup"`
"""
function check_mpo_algo(mpo_algo::String)
    if mpo_algo ∉ ("naive", "zipup")
        error("mpo_algo is \"naive\" or \"zipup\", not $(repr(mpo_algo))")
    end
end

"""
    steady_state(lindbladian, ::State; options...)
    steady_state(lindbladian, ::Simulation; options...)

the steady state of a Lindbladian ``L``, of the form `-im * H` plus dissipators, or its
`PreMPO`, by dmrg on ``L^\\dagger L`` starting from the given mixed state. It is returned as `(value, state)`, or
`(value, simulation)`, where `value` is the "energy" dmrg reaches, zero for a steady state,
and the state is normalized to trace one.

# Options

- `nsweeps`: the last sweep to do, that is the number of sweeps of the whole run (default 1)
- `first_sweep`: the sweep to start from (default 1), to continue a search left unfinished
- `limits`: constraints on the state, see `Limits`, which may give one value per sweep
  (default `Limits()`)
- `mpo_limits`: the truncation of the MPO of ``L^\\dagger L`` (default `Limits()`)
- `mpo_algo`: the algorithm computing ``L^\\dagger L``, `"naive"` (default) or `"zipup"`
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `observer!`: an observer, see `DmrgObserver`
- `krylov`: the parameters of the Krylov search of each local step of dmrg, see `Krylov`
  (default `Krylov(dim = 8, maxiter = 3)`). The spectrum of ``L^\\dagger L`` crowds near zero,
  which the default of `ITensorMPS.dmrg`, three vectors, does not resolve: the search then
  stalls at a residual well above rounding, on a state that is not steady

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
    # an eigenvector of (L+)L has norm one and a sign of its own, the trace set to one makes it
    # the density matrix it stands for
    e, st = dmrg(l2, state; nsweeps, first_sweep, limits, noise, observer!, krylov)
    return (e, normalize(st))
end

"""
    thermal_state(hamiltonian, beta, ::State; options...)
    thermal_state(hamiltonian, beta, ::Simulation; options...)

the mixed state ``\\rho`` taken to ``e^{-\\beta H/2} \\rho \\, e^{-\\beta H/2}``, normalized to
trace one, by tdvp in imaginary time, in `nsteps` steps of `beta / nsteps`. It is returned as
`(log_trace, state)`, or `(log_trace, simulation)`, the time of the simulation unchanged, where
`log_trace` is ``\\log \\mathrm{tr}(e^{-\\beta H/2} \\rho \\, e^{-\\beta H/2})``, ``\\rho`` being
the state given normalized to trace one: each step is normalized, so that the trace, which
grows or shrinks exponentially with `beta`, never leaves the range of the numbers, and its
logarithms are summed.

From `"FullyMixed"`, the state at infinite temperature, it gives the thermal state
``e^{-\\beta H}/Z``, and `log_trace` is ``\\log Z - \\sum_i \\log d_i``, ``d_i`` being the
dimension of site `i`. More generally a state that commutes with ``H`` gives the thermal state
restricted to what it describes: `fully_mixed(system, N => m)` the canonical thermal state
of `m` particles, `"MixedSpin"` on `Tj` sites, at one electron per site,
that of the Heisenberg model, and a product state of populations ``e^{\\beta\\mu n_i}``
the grand canonical state at chemical potential ``\\mu``, for a hamiltonian conserving the
number of particles. A pure state is refused: ``e^{-\\beta H/2}`` takes it towards the ground
state, see `dmrg`.

# Options

- `nsteps`: the number of steps (default 1)
- `first_step`: the step to start from (default 1), to continue a computation left unfinished
  from the state it had reached: `beta` and `nsteps` are still those of the whole computation
- `log_trace`: the value of `log_trace` the computation starts from (default 0), that of the
  step before `first_step` to continue it
- `expand_period`: enlarge the bond dimension of the state by a global Krylov expansion before the
  first step and then every `expand_period` steps (default 0, never)
- `hermitianize_period`: make the state hermitian every `hermitianize_period` steps (default 0,
  never)
- `limits`: constraints on the state, see `Limits`, which may give one value per step
  (default `Limits()`)
- `observer!`: an observer, see `ThermalObserver`
- `krylov`: the parameters of the Krylov exponentiation of each local step, see `Krylov`
  (default `Krylov()`, those of `KrylovKit.exponentiate`)

# Examples

    log_trace, rho = thermal_state(H, 2.0, State{Mixed}(system, "FullyMixed"); nsteps = 20,
                                   limits = Limits(cutoff = 1e-12, maxdim = 64))
"""
function thermal_state(hamiltonian::IndexedOp{Pure}, beta::Real, state::State{Mixed};
    observer! = NoObserver(), nsteps::Int = 1, first_step::Int = 1, log_trace::Real = 0.,
    expand_period::Int = 0, hermitianize_period::Int = 0, limits::Limits = Limits(),
    krylov::Krylov = Krylov())
    check_nsweeps(nsteps, "nsteps")
    dbeta = beta / nsteps
    # a placed operator on a mixed state becomes its Evolver, A ρ + ρ A†, which for A = -H/2
    # is the generator of e^{-βH/2} ρ e^{-βH/2}
    mpo = make_mpo(state, -hamiltonian / 2)
    updater_kwargs = (; eager = true, krylov_kwargs(krylov)...)
    # a computation continued goes on from the state it had reached, already normalized:
    # normalized again, its rounding would differ from that of the uninterrupted one
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
