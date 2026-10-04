# The Thermalize phase, which takes a mixed state towards the thermal state of a hamiltonian.

export Thermalize

"""
    Thermalize(; hamiltonian, beta, beta_step, algo, limits, measurements, measurements_period,
                 options...)

a phase that takes a mixed state ``\\rho`` to ``e^{-\\beta H/2} \\rho \\, e^{-\\beta H/2}``,
normalized to trace one, by tdvp in imaginary time, see `thermal_state`. After
`CreateState{Mixed}(…, "FullyMixed")`, the state at infinite temperature, it gives the thermal
state ``e^{-\\beta H}/Z``; a state that commutes with ``H`` gives the thermal state restricted
to what it describes, as `fully_mixed` the canonical one. The time of the simulation does not
move.

# Fields

- `name`, `time_start`, `final_measurements`: the fields every phase has, see `AbstractPhase`
- `hamiltonian`: the hamiltonian ``H``
- `beta`: the inverse temperature ``\\beta``
- `beta_step`: the step of `beta`, adjusted to the nearest one that divides `beta` into a whole
  number of steps, and taken with the sign of `beta` (the phase is skipped when that number is
  zero)
- `algo`: the algorithm, `Tdvp(...)`, whose `expand_period`, `hermitianize_period` and `krylov`
  it takes (default `Tdvp()`)
- `limits`: constraints on the state, see `Limits` (default `Limits()`)
- `measurements`: the measurements to make during the phase, see `output` (default `[]`), after
  every `measurements_period` steps. The symbols `:beta` and `:log_trace` take the inverse
  temperature reached and the logarithm of the trace the state would have without being
  normalized, ``\\log Z - \\sum_i \\log d_i`` from `"FullyMixed"`, ``d_i`` being the dimension
  of site `i`, see `thermal_state`
- `measurements_period`: the number of steps between two measurements (default 1)

# Examples

    Thermalize(hamiltonian = -sum(Z(i)Z(i+1) for i in 1:9) - sum(X(i) for i in 1:10),
               beta = 2., beta_step = 0.1, limits = Limits(cutoff = 1e-12, maxdim = 64),
               measurements = "thermo" => [:beta, :log_trace, Z(5)])
"""
@kwdef struct Thermalize <: AbstractPhase
    name::String = "Thermalization"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    hamiltonian::IndexedOp{Pure}
    beta::Real
    beta_step::Real
    algo::Algo = Tdvp()
    limits::Limits = Limits()
    measurements = []
    measurements_period::Int = 1
    # refused here rather than when it runs, as the fields of Evolve
    function Thermalize(name, time_start, final_measurements, hamiltonian, beta, beta_step, algo,
                        limits, measurements, measurements_period)
        if iszero(beta_step)
            error("the beta_step of a Thermalize cannot be zero")
        end
        if !(algo isa Tdvp)
            error("Thermalize computes with Tdvp(...), not $algo")
        end
        return new(name, time_start, final_measurements, hamiltonian, beta, beta_step, algo, limits,
                   measurements, measurements_period)
    end
end

function run_phase(sim::Simulation, phase::Thermalize)
    # adjusted as the time step of Evolve, the step taking the sign of beta
    nsteps = round(Int, abs(phase.beta / phase.beta_step))
    if nsteps == 0
        log_message(sim, "Skipping a thermalization to beta $(phase.beta), shorter than half " *
                         "a step")
        return sim
    end
    beta = phase.beta
    if !(beta / nsteps ≈ phase.beta_step)
        log_message(sim, "Taking a step of $(beta / nsteps) rather than $(phase.beta_step) to " *
                         "reach beta $beta in $nsteps steps")
    end
    done, log_trace = resume_step(sim)
    log_message(sim, "Thermalizing state to beta $beta")
    algo = phase.algo
    l, sim = thermal_state(phase.hamiltonian, beta, sim; nsteps, first_step = done + 1,
                           log_trace = something(log_trace, 0.), algo.expand_period,
                           algo.hermitianize_period, phase.limits, algo.krylov,
                           observer! = ThermalObserver(sim, phase.measurements,
                                                       phase.measurements_period))
    log_message(sim, "Done, log_trace is $l")
    return sim
end
