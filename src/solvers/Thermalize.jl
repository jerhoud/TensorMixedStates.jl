# The thermal states: thermalize, an evolution in imaginary time of a density matrix by any method
# of evolution, normalized at each step.

export thermalize

"""
    thermal_step(algo, pre, state, dbeta, step, limits)

the density matrix `state`, of trace one, taken to ``e^{-d\\beta H/2} \\rho \\, e^{-d\\beta H/2}``
by one step of the method `algo`, `pre` being the `PreMPO` of `H` on the mixed state, that of
``\\rho \\mapsto H\\rho + \\rho H``, and the logarithm of its trace, by which it is divided
"""
function thermal_step(algo::EvolutionAlgorithm, pre::PreMPO{Mixed}, state::State{Mixed},
                      dbeta::Real, step::Int, limits::Limits)
    state = evolve_step(algo, pre, nothing, state, 0., -dbeta / 2, step, limits)
    t = nonzero_trace(state, "thermalize")
    if !(t > 0)
        error("the trace of the state reached $t at step $step of thermalize: take smaller steps " *
              "or larger limits")
    end
    return (state / t, log(t))
end

"""
    thermalize(hamiltonian, beta, state; algo = Tdvp(), nsteps = 1, limits = Limits(), step = 1)

the mixed state ``\\rho``, normalized to trace one, taken to
``e^{-\\beta H/2} \\rho \\, e^{-\\beta H/2}`` by an evolution in imaginary time in `nsteps` steps
of `beta / nsteps` by the method `algo`, see `EvolutionAlgorithm` (`ApproxW` has its limits
there), the state normalized to trace one at each step. It is returned as `(log_trace, state)`,
`log_trace` being ``\\log \\mathrm{tr}(e^{-\\beta H/2} \\rho \\, e^{-\\beta H/2})`` for ``\\rho``
of trace one, so that the `log_trace` of successive calls add up.

From `"FullyMixed"`, it gives the thermal state ``e^{-\\beta H}/Z``, and `log_trace` is
``\\log Z - \\sum_i \\log d_i``, ``d_i`` being the dimension of site `i`. A state commuting
with ``H`` gives the thermal state restricted to what it describes: `fully_mixed(system,
N => m)` the canonical state of `m` particles, a product state of populations
``e^{\\beta\\mu n_i}`` the grand canonical state at chemical potential ``\\mu``. The
Hamiltonian may be given as its `PreMPO` on the mixed state, prepared once and keeping its MPOs.
`step` is the number of the first step, for the periods of the method and for `limits`, which
may give one value per step; from a later step, the state is taken as a computation left it, of
trace one, and not normalized again. A pure state is refused.

# Examples

    log_trace, rho = thermalize(H, 2.0, State{Mixed}(system, "FullyMixed"); nsteps = 20,
                                limits = Limits(cutoff = 1e-12, maxdim = 64))
"""
function thermalize(hamiltonian, beta::Real, state::State{Mixed};
                    algo::EvolutionAlgorithm = Tdvp(), nsteps::Int = 1,
                    limits::Limits = Limits(), step::Int = 1)
    if nsteps < 1
        error("thermalize takes at least one step, and nsteps is $nsteps")
    end
    pre = prepared_operator(state, hamiltonian)
    # a continued computation is not normalized again, which would change its rounding, and the
    # basis an expansion picks with it
    if step == 1
        state = normalize(state)
    end
    log_trace = 0.
    for k in 0:nsteps - 1
        state, l = thermal_step(algo, pre, state, beta / nsteps, step + k,
                                step_limits(limits, step + k))
        log_trace += l
    end
    return (log_trace, state)
end

thermalize(_, ::Real, ::State{Pure}; kwargs...) =
    error("a thermal state needs a mixed representation: mix the state first")
