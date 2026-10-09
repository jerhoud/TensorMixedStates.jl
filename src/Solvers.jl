# The algorithms on a state or a simulation: tdvp and approx_W for time evolution, dmrg for
# ground states, and steady_state for the steady state of an open system.

export tdvp, dmrg, approx_W, steady_state

"""
    check_nsweeps(nsweeps)

refuse an evolution in less than one step: it made no step, while a simulation took the time
on, and `approx_W` divided its duration by zero
"""
function check_nsweeps(nsweeps)
    if nsweeps < 1
        error("an evolution takes at least one step, and nsweeps is $nsweeps")
    end
end

"""
    tdvp(evolver, t, ::State; options...)
    tdvp(evolver, t, ::Simulation; options...)

evolve a state, or a simulation, for a time `t` with the tdvp algorithm, in `nsweeps` steps of
`t / nsweeps`. `evolver` is `-im * H` for a hamiltonian `H`, plus dissipators for a mixed
state. A simulation comes back with its time advanced by `t`.

# Options

- `nsweeps`: the number of steps (default 1)
- `first_sweep`: the step to start from (default 1), to continue an evolution left unfinished:
  `t`, `nsweeps` and `time_start` are still those of the whole evolution, and a simulation
  is given at the time that evolution started
- `time_start`: the simulation time the evolution starts from (default 0, and the time of
  the simulation for a `Simulation`)
- `coefs`: for a vector of evolvers, the functions of time they are multiplied by, taken at
  the middle of each step, with real values on a mixed state
- `n_expand`: enlarge the bond dimension of the state by a global Krylov expansion every
  `n_expand` steps (default 0, never)
- `n_hermitianize`: make a mixed state hermitian every `n_hermitianize` steps (default 0,
  never)
- `limits`: constraints on the state, see `Limits`, which may give one value per step
  (default `Limits()`, none)
- `observer!`: an observer, see `TdvpObserver`
- the other options are passed to `ITensorMPS.tdvp`

# Examples

    tdvp(-im * H, 1., state; nsweeps = 10, limits = Limits(cutoff = 1e-10, maxdim = 50))
"""
function tdvp(pre::PreMPO{R}, t::Number, state::State{R};
    observer! = NoObserver(), coefs=nothing, n_expand = 0, n_hermitianize = 0,
    nsweeps = 1, first_sweep = 1, time_start = zero(t), limits::Limits=Limits(), kwargs...) where {R <: PM}
    check_nsweeps(nsweeps)
    time_dep = !isnothing(coefs)
    st = state.state
    dt = t / nsweeps
    if !time_dep
        mpo = make_mpo(pre)
    end
    for sweep in first_sweep:nsweeps
        current_time = time_start + sweep * dt
        if time_dep
            tf = current_time - dt / 2
            mpo = make_mpo(pre, map(f->f(tf), coefs))
        end
        lim = sweep_limits(limits, sweep)
        st = tdvp(mpo, dt, st; nsweeps = 1, lim.cutoff, lim.maxdim, lim.mindim, kwargs...)
        if sweep_due(n_hermitianize, sweep)
            st = hermitianize(State(state, st); limits = lim).state
        end
        if sweep_due(n_expand, sweep)
            st = expand(st, mpo; alg="global_krylov")
        end
        # once the whole sweep is done, the expansion included: the sweep is committed right
        # after its measurements, and a checkpoint cannot fall between them
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
`(energy, state)`, or `(energy, simulation)`. A hamiltonian is refused on a mixed state, where
the lowest eigenvector of the superoperator it gives is neither the ground state nor a
density matrix: search the ground state of the pure state, then `mix` it.

# Options

- `nsweeps`: the last sweep to do, that is the number of sweeps of the whole run (default 1)
- `first_sweep`: the sweep to start from (default 1), to continue a search left unfinished;
  past `nsweeps`, no sweep is done and the energy is that of the state given
- `limits`: constraints on the state, see `Limits`, which may give one value per sweep
  (default `Limits()`, none)
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `observer!`: an observer, see `DmrgObserver`
- the other options are passed to `ITensorMPS.dmrg`, `outputlevel` defaulting here to 0

# Examples

    energy, state = dmrg(H, state; nsweeps = 10, limits = Limits(maxdim = [10, 20, 50]))
"""
function dmrg(mpo::MPO, state::State; nsweeps = 1, first_sweep = 1, observer! = NoObserver(),
              limits::Limits = Limits(), noise = 0., kwargs...)
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
                 noise = resume_schedule(noise, done), kwargs...)
    return (e, State(state, st))
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
    make_approx_W(pre, t; order, w, coefs = [1.])

the MPOs of the approximation of the given `order` of a step `t`, to apply one after the
other, built from WI (`w = 1`) or WII (`w = 2`) approximations; `coefs` are the values of the
coefficients of a time dependent evolver.
"""
function make_approx_W(pre::PreMPO, t::Number; order::Int, w::Int, coefs = [1.])
    if order < 1 || order > length(w_approx_coefs)
        error("W approximation of order $order is not implemented")
    end
    if w ∉ (1, 2)
        error("W approximation is only defined for w=1 or 2 (not $w)")
    end
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

- `order`: the order of the approximation, from 1 to 4, required; 2 at most with `coefs`, the
  time functions being taken at the middle of each step
- `w`: 1 or 2 for WI or WII (default 2)
- `nsweeps`: the number of steps (default 1)
- `first_sweep`: the step to start from (default 1), to continue an evolution left unfinished:
  `t`, `nsweeps` and `time_start` are still those of the whole evolution, and a simulation
  is given at the time that evolution started
- `time_start`: the simulation time the evolution starts from (default 0, and the time of
  the simulation for a `Simulation`)
- `coefs`: for a vector of evolvers, the functions of time they are multiplied by, taken at
  the middle of each step, with real values on a mixed state
- `n_hermitianize`: make a mixed state hermitian every `n_hermitianize` steps (default 0,
  never)
- `limits`: constraints on the state, see `Limits`, which may give one value per step
  (default `Limits()`, none)
- `observer!`: an observer, see `ApproxWObserver`
- the other options are passed to the MPO product of ITensorMPS, `ITensorMPS.apply`

# Examples

    approx_W(-im * H, 1., state; order = 4, nsweeps = 10)
"""
function approx_W(pre::PreMPO{R}, t::Number, state::State{R}; coefs = nothing, n_hermitianize::Int = 0,
    nsweeps::Int = 1, first_sweep::Int = 1, order::Int, w::Int = 2, observer! = NoObserver(),
    time_start = zero(t), limits::Limits=Limits(), kwargs...) where {R <: PM}
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
            st = apply(mpo, st; lim.cutoff, lim.maxdim, lim.mindim, kwargs...)
        end
        if sweep_due(n_hermitianize, sweep)
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
    steady_state(lindbladian, ::State; options...)
    steady_state(lindbladian, ::Simulation; options...)

the steady state of a Lindbladian ``L``, of the form `-im * H` plus dissipators, by dmrg on
``L^\\dagger L`` starting from the given mixed state. It is returned as `(value, state)`, or
`(value, simulation)`, where `value` is the "energy" dmrg reaches, zero for a steady state,
and the state is normalized to trace one.

# Options

- `nsweeps`: the last sweep to do, that is the number of sweeps of the whole run (default 1)
- `first_sweep`: the sweep to start from (default 1), to continue a search left unfinished
- `limits`: constraints on the state, see `Limits`, which may give one value per sweep
  (default `Limits()`, none)
- `mpo_limits`: the truncation of the MPO of ``L^\\dagger L`` (default `Limits()`, none)
- `mpo_algo`: the algorithm computing ``L^\\dagger L``, `"naive"` (default) or `"zipup"`
- `observer!`: an observer, see `DmrgObserver`
- the other options are passed to `ITensorMPS.dmrg`, `outputlevel` defaulting here to 0

# Examples

    value, rho = steady_state(-im * H + D, rho; nsweeps = 10,
                              limits = Limits(maxdim = [10, 20, 50]))
"""
function steady_state(op::IndexedOp{Mixed}, state::State{Mixed};
    limits::Limits = Limits(), nsweeps::Int = 1,
    observer! = NoObserver(), mpo_limits::Limits = Limits(), mpo_algo::String = "naive",
    alg = nothing, kwargs...)
    if !isnothing(alg)
        @warn "the `alg` keyword of steady_state is now `mpo_algo`, matching the field of " *
              "the SteadyState phase. The old name still works and will be removed." maxlog = 1
        mpo_algo = alg
    end
    l = make_mpo(state, op)
    # the naive algorithm truncates only when asked to, the others take no such option
    extra = mpo_algo == "naive" ? (; truncate = mpo_limits != Limits()) : (;)
    l2 = apply(replaceprime(dag(l)', 2=>0), l;
               mpo_limits.cutoff, mpo_limits.maxdim, mpo_limits.mindim, alg = mpo_algo, extra...)
    # an eigenvector of (L+)L has norm one and a sign of its own, the trace set to one makes it
    # the density matrix it stands for
    e, st = dmrg(l2, state; nsweeps, limits, observer!, kwargs...)
    return (e, normalize(st))
end

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
