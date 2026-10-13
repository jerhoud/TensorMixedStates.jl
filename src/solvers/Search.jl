# The searches of a state by dmrg: ground_state, the lowest eigenvector of a Hamiltonian, and
# steady_state, that of L†L for a Lindbladian L, through search, which each method implements,
# calling a function at the end of each sweep.

export ground_state, steady_state

"""
    SweepHook(tol, energy, done, nsweeps, state, mpo, sweep_done)

the observer of a dmrg search from `state` of `nsweeps` sweeps: at the end of each sweep,
`sweep_done(sweep, energy, state, ends)`, unless it is `nothing`, the sweeps counted after the
`done` ones a resumed search had made, `ends` telling whether the search ends there, and a stop
once the energy changes by less than `tol` from one sweep to the next, `energy` being that of
the last. Given `mpo`, the energy is that of `mpo` on the state, and not the one dmrg gives,
which holds the penalty of a search orthogonal to other states.
"""
mutable struct SweepHook <: ITensorMPS.AbstractObserver
    tol::Float64
    energy::Union{Nothing, Float64}
    done::Int
    nsweeps::Int
    state::State
    mpo::Union{Nothing, MPO}
    sweep_done
end

function ITensorMPS.checkdone!(o::SweepHook; energy, sweep, psi, kwargs...)
    # normalized as the search normalizes what it returns, an eigenvector having a sign of its
    # own; normalize makes a new MPS, leaving that of dmrg untouched
    st = isnothing(o.sweep_done) && isnothing(o.mpo) ? nothing : normalize(State(o.state, psi))
    e = isnothing(o.mpo) ? energy : real(inner(st.state', o.mpo, st.state))
    stop = !isnothing(o.energy) && abs(o.energy - e) < o.tol
    o.energy = e
    if !isnothing(o.sweep_done)
        o.sweep_done(sweep + o.done, e, st, stop || sweep + o.done == o.nsweeps)
    end
    return stop
end

"""
    search(algo, mpo, state, limits; sweep_done = nothing, krylov = Krylov(), done = 0,
           energy = nothing[, orthogonal_to, weight])

the lowest eigenvalue of `mpo` and its eigenvector, normalized, searched from `state` by the
method `algo`, `limits` the limits of the sweeps, see `Limits`; `sweep_done(sweep, energy,
state, ends)` is called at the end of each sweep, `ends` telling whether the search ends there,
on its tolerance or its last sweep, with the state it returns, and `krylov` gives the defaults
of the method of the caller. A search
resumed after `done` sweeps takes the rest of the schedules, `energy` being that of the last
sweep, which it returns when no sweep is left. Given `orthogonal_to`, pure states on the system
of `state`, the search is that of the lowest eigenvector orthogonal to them, which dmrg makes
by adding `weight` times their projectors to `mpo`, the energy being that of `mpo` alone. A
`SearchAlgorithm` of one's own adds a method, taking these two keywords only if it searches
excited states.
"""
function search(algo::Dmrg, mpo::MPO, state::State, limits::Limits; sweep_done = nothing,
                krylov::Krylov = Krylov(), done::Int = 0, energy = nothing,
                orthogonal_to::Vector{State{Pure}} = State{Pure}[], weight::Real = 0.)
    # ITensorMPS, asked for no sweep, gives an energy of 0; a resumed search gives back that of its
    # last sweep, as the search it resumes
    if done ≥ algo.nsweeps
        st = normalize(state)
        return (something(energy, real(inner(st.state', mpo, st.state))), st)
    end
    # ITensorMPS counts its sweeps from 1: a resumed search asks for those it has left, with the
    # tail of its schedules
    excited = !isempty(orthogonal_to)
    hook = SweepHook(algo.tol, energy, done, algo.nsweeps, state, excited ? mpo : nothing,
                     sweep_done)
    kwargs = (; outputlevel = 0, nsweeps = algo.nsweeps - done, observer = hook,
              truncation(later_steps(limits, done))..., noise = later_steps(algo.noise, done),
              krylov_kwargs(algo.krylov, krylov, "eigsolve_")...)
    e, st = excited ? dmrg(mpo, [ m.state for m in orthogonal_to ], state.state; weight, kwargs...) :
                      dmrg(mpo, state.state; kwargs...)
    return (excited ? hook.energy : e, normalize(State(state, st)))
end

"""
    search_algo(algo, what)

the search method `algo` of `what`, which has no default, refused with what to write when it is
not given
"""
function search_algo(algo, what::String)
    if isnothing(algo)
        error("$what needs its search method: algo = Dmrg(nsweeps = 10), for instance")
    end
    return algo
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
    orthogonal_args(state, pre, orthogonal_to, weight)

the keywords of `search` for a state orthogonal to the pure states `orthogonal_to`, on the system
of `state`, none when there are none: the states of another sector of charge than `state` are
left out, orthogonal to it already, and the weight is by default the bound of the width of the
spectrum of the Hamiltonian of `pre`, see `term_width`, above any gap between its eigenvalues
"""
function orthogonal_args(state::State{Pure}, pre::PreMPO{Pure}, orthogonal_to::AbstractVector,
                         weight::Union{Nothing, Real})
    ms = State{Pure}[]
    for m in orthogonal_to
        if !(m isa State{Pure})
            error("ground_state searches a state orthogonal to pure states, not to a $(typeof(m))")
        elseif m.system ≠ state.system
            error("a state of orthogonal_to is on another system than that of the state")
        end
        if ITensorMPS.flux(m.state) == ITensorMPS.flux(state.state)
            push!(ms, m)
        end
    end
    return isempty(ms) ? (;) : (; orthogonal_to = ms, weight = something(weight, pre.width[]))
end

"""
    check_orthogonal(state, args)

warn when `state`, found by a search orthogonal to other states, see `orthogonal_args`, keeps
an overlap with one of them beyond the precision of the numbers: a weight below the gap, or a
search that has not converged
"""
function check_orthogonal(state::State{Pure}, args::NamedTuple)
    for m in get(args, :orthogonal_to, State{Pure}[])
        p = abs2(inner(m, state))
        if p > sqrt(eps())
            @warn "the state found keeps an overlap of $(round(p; sigdigits = 2)) with a state " *
                  "of orthogonal_to: the weight $(args.weight) is below the gap, or the search " *
                  "has not converged"
        end
    end
end

"""
    ground_state(hamiltonian, state; algo, limits = Limits(), orthogonal_to = [], weight = nothing)

the ground state of the Hamiltonian, searched from the pure `state` by the method `algo`, see
`SearchAlgorithm`, as `(energy, state)`. The Hamiltonian may be given as its
`PreMPO(state, hamiltonian)`, prepared once for several searches and keeping its MPO. `limits`,
see `Limits`, may give one value per sweep, with a default cutoff of `1e-26`. A mixed state is
refused: search the pure state and mix it; the steady state of a Lindbladian is given by
`steady_state`.

Given `orthogonal_to`, eigenstates of the Hamiltonian found before, the search gives the
lowest eigenstate orthogonal to them, an excited state, `weight` times their projectors being
added to the Hamiltonian; it must exceed the gap between the state searched and each of them,
and is by default a bound of the width of the spectrum, sure but slow to converge, a weight just
above the gap converging faster. Its precision is that of the states given, and the penalty is
exact only for eigenstates. A warning says when the state found keeps an overlap with one of
them. The energy is that of the Hamiltonian alone.

# Examples

    energy, state = ground_state(H, state; algo = Dmrg(nsweeps = 10),
                                 limits = Limits(maxdim = [10, 20, 50]))
    e1, excited = ground_state(H, state; algo = Dmrg(nsweeps = 20), orthogonal_to = [ground])
"""
function ground_state(hamiltonian, state::State{Pure}; algo::Union{Nothing, SearchAlgorithm} = nothing,
                      limits::Limits = Limits(), orthogonal_to::AbstractVector = State{Pure}[],
                      weight::Union{Nothing, Real} = nothing)
    algo = search_algo(algo, "ground_state")
    pre = prepared_operator(state, hamiltonian)
    args = orthogonal_args(state, pre, orthogonal_to, weight)
    e, st = search(algo, cached_mpos(() -> make_mpo(pre), pre, :mpo, nothing), state, limits; args...)
    check_orthogonal(st, args)
    return (e, st)
end

ground_state(_, ::State{Mixed}; kwargs...) =
    error("ground_state searches a pure state: search it pure and mix it")

"""
    check_mpo_algo(mpo_algo)

refuse an algorithm of the product ``L^\\dagger L`` other than `"naive"` and `"zipup"`
"""
function check_mpo_algo(mpo_algo::String)
    if mpo_algo ∉ ("naive", "zipup")
        error("mpo_algo is \"naive\" or \"zipup\", not $(repr(mpo_algo))")
    end
    return mpo_algo
end

"""
    squared_lindbladian(pre, mpo_limits, mpo_algo)

the MPO of ``L^\\dagger L`` for the Lindbladian of `pre`, computed by `mpo_algo` and truncated by
`mpo_limits`, kept in `pre`
"""
function squared_lindbladian(pre::PreMPO{Mixed}, mpo_limits::Limits, mpo_algo::String)
    check_mpo_algo(mpo_algo)
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
function steady_state(lindbladian, state::State{Mixed}; algo::Union{Nothing, SearchAlgorithm} = nothing,
                      limits::Limits = Limits(), mpo_limits::Limits = Limits(),
                      mpo_algo::String = "naive")
    algo = search_algo(algo, "steady_state")
    l2 = squared_lindbladian(prepared_operator(state, lindbladian), mpo_limits, mpo_algo)
    return steady_result(search(algo, l2, state, limits; krylov = steady_krylov)...)
end

steady_state(_, ::State{Pure}; kwargs...) =
    error("steady_state searches a density matrix: mix the state")
