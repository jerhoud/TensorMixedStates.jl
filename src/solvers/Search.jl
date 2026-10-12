# The searches of a state by dmrg: ground_state, the lowest eigenvector of a Hamiltonian, and
# steady_state, that of L†L for a Lindbladian L, through search, which each method implements,
# calling a function at the end of each sweep.

export ground_state, steady_state

"""
    SweepHook(tol, energy, done, state, sweep_done)

the observer of a dmrg search from `state`: at the end of each sweep, `sweep_done(sweep,
energy, state)`, unless it is `nothing`, the sweeps counted after the `done` ones a resumed
search had made, and a stop once the energy changes by less than `tol` from one sweep to the
next, `energy` being that of the last
"""
mutable struct SweepHook <: ITensorMPS.AbstractObserver
    tol::Float64
    energy::Union{Nothing, Float64}
    done::Int
    state::State
    sweep_done
end

function ITensorMPS.checkdone!(o::SweepHook; energy, sweep, psi, kwargs...)
    stop = !isnothing(o.energy) && abs(o.energy - energy) < o.tol
    o.energy = energy
    if !isnothing(o.sweep_done)
        # normalized as the search normalizes what it returns, an eigenvector having a sign of
        # its own; normalize makes a new MPS, leaving that of dmrg untouched
        o.sweep_done(sweep + o.done, energy, normalize(State(o.state, psi)))
    end
    return stop
end

"""
    search(algo, mpo, state, limits; sweep_done = nothing, krylov = Krylov(), done = 0,
           energy = nothing)

the lowest eigenvalue of `mpo` and its eigenvector, normalized, searched from `state` by the
method `algo`, `limits` the limits of the sweeps, see `Limits`; `sweep_done(sweep, energy,
state)` is called at the end of each sweep, and `krylov` gives the defaults of the method of
the caller. A search resumed after `done` sweeps takes the rest of the schedules, `energy` being
that of the last sweep. A `SearchAlgorithm` of one's own adds a method.
"""
function search(algo::Dmrg, mpo::MPO, state::State, limits::Limits; sweep_done = nothing,
                krylov::Krylov = Krylov(), done::Int = 0, energy = nothing)
    # ITensorMPS, asked for no sweep, gives an energy of 0
    if done ≥ algo.nsweeps
        st = normalize(state)
        return (real(inner(st.state', mpo, st.state)), st)
    end
    # ITensorMPS counts its sweeps from 1: a resumed search asks for those it has left, with the
    # tail of its schedules
    hook = SweepHook(algo.tol, energy, done, state, sweep_done)
    e, st = dmrg(mpo, state.state; outputlevel = 0, nsweeps = algo.nsweeps - done,
                 observer = hook, truncation(later_steps(limits, done))...,
                 noise = later_steps(algo.noise, done),
                 krylov_kwargs(algo.krylov, krylov, "eigsolve_")...)
    return (e, normalize(State(state, st)))
end

"""
    prepared_operator(state, op)

the `PreMPO` of the operator `op` of a search on `state`, or `op` itself when it is one
"""
prepared_operator(state::State, a::IndexedOp) = PreMPO(state, a)
prepared_operator(state::State{R}, pre::PreMPO{R}) where {R <: PM} = check_pre(pre, state)

prepared_operator(::State, ::PreMPO) =
    error("the PreMPO was prepared on a state of the other representation")

prepared_operator(::State, ::MPO) =
    error("the solvers take an operator or its PreMPO, which keeps the MPOs they build, not an MPO")

prepared_operator(::State, x) = error("$(repr(x)) is not an operator to search the state of")

"""
    ground_state(hamiltonian, state; algo, limits = Limits())

the ground state of the Hamiltonian, searched from the pure `state` by the method `algo`, see
`SearchAlgorithm`, as `(energy, state)`. The Hamiltonian may be given as its
`PreMPO(state, hamiltonian)`, prepared once for several searches and keeping its MPO. `limits`,
see `Limits`, may give one value per sweep, with a default cutoff of `1e-26`. A mixed state is
refused: search the pure state and mix it; the steady state of a Lindbladian is given by
`steady_state`.

# Examples

    energy, state = ground_state(H, state; algo = Dmrg(nsweeps = 10),
                                 limits = Limits(maxdim = [10, 20, 50]))
"""
function ground_state(hamiltonian, state::State{Pure}; algo::SearchAlgorithm,
                      limits::Limits = Limits())
    pre = prepared_operator(state, hamiltonian)
    return search(algo, cached_mpos(() -> make_mpo(pre), pre, :mpo, nothing), state, limits)
end

ground_state(_, ::State{Mixed}; kwargs...) =
    error("ground_state searches a pure state: search it pure and mix it")

"""
    squared_lindbladian(pre, mpo_limits, mpo_algo)

the MPO of ``L^\\dagger L`` for the Lindbladian of `pre`, computed by `mpo_algo` and truncated by
`mpo_limits`, kept in `pre`
"""
function squared_lindbladian(pre::PreMPO{Mixed}, mpo_limits::Limits, mpo_algo::String)
    if mpo_algo ∉ ("naive", "zipup")
        error("mpo_algo is \"naive\" or \"zipup\", not $(repr(mpo_algo))")
    end
    return cached_mpos(pre, :squared_lindbladian, (mpo_limits, mpo_algo)) do
        l = cached_mpos(() -> make_mpo(pre), pre, :mpo, nothing)
        # the naive algorithm truncates only when asked to, the others take no such option
        extra = mpo_algo == "naive" ? (; truncate = mpo_limits ≠ Limits()) : (;)
        apply(replaceprime(dag(l)', 2 => 0), l; truncation(mpo_limits)..., alg = mpo_algo, extra...)
    end
end

"""
    steady_krylov

the defaults of the Krylov search of `steady_state`, which the spectrum of ``L^\\dagger L``,
crowded near zero, needs
"""
const steady_krylov = Krylov(dim = 8, maxiter = 3)

"""
    steady_result(value, state)

what `steady_state` returns of the eigenvector found: the state divided by its trace, a complex
number, which takes away the phase that dividing by its real part, as normalize does, would
leave, with a warning when it is not hermitian
"""
function steady_result(value, state::State{Mixed})
    ρ = state / trace(state)
    δ = 1 - hermiticity(ρ)
    if δ > 1e-6
        @warn "the steady state found has a HermiticityDeviation of $(round(δ; sigdigits = 2)): " *
              "the search has not converged, or the steady state is not unique, see steady_state"
    end
    return (value, ρ)
end

"""
    steady_state(lindbladian, state; algo, limits = Limits(), mpo_limits = Limits(),
                 mpo_algo = "naive")

the steady state of the Lindbladian ``L``, `-im * H` plus dissipators, the lowest eigenvector
of ``L^\\dagger L`` searched from the mixed `state` by the method `algo`, see
`SearchAlgorithm`, as `(value, state)` with trace one, `value` being the eigenvalue reached,
zero for a steady state. The Lindbladian may be given as its `PreMPO(state, lindbladian)`,
prepared once and keeping the MPO of ``L^\\dagger L``.

- `limits`: the limits of the sweeps, see `Limits`, with a default cutoff of `1e-26`
- `mpo_limits`: the truncation of the MPO of ``L^\\dagger L`` (default none)
- `mpo_algo`: the algorithm of the product ``L^\\dagger L``, `"naive"` or `"zipup"`

The Krylov search of `Dmrg` takes `dim = 8` and `maxiter = 3` by default, which the spectrum of
``L^\\dagger L``, crowded near zero, needs. A Lindbladian may have several steady states, one in
each sector of a quantity it conserves for instance, and dmrg returns any combination of them,
in general not a density matrix; declaring the quantity `strong` keeps the search in the
sector of the starting state. A warning is given when the state found is not hermitian, its
`HermiticityDeviation` above `1e-6`: the search has not converged, or the steady state is not
unique.

# Examples

    value, rho = steady_state(-im * H + D, rho; algo = Dmrg(nsweeps = 10),
                              limits = Limits(maxdim = [10, 20, 50]))
"""
function steady_state(lindbladian, state::State{Mixed}; algo::SearchAlgorithm,
                      limits::Limits = Limits(), mpo_limits::Limits = Limits(),
                      mpo_algo::String = "naive")
    l2 = squared_lindbladian(prepared_operator(state, lindbladian), mpo_limits, mpo_algo)
    return steady_result(search(algo, l2, state, limits; krylov = steady_krylov)...)
end

steady_state(_, ::State{Pure}; kwargs...) =
    error("steady_state searches a density matrix: mix the state")
