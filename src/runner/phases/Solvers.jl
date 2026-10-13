# The phases of the solvers: Evolve and Thermalize, made of steps, and GroundState and
# SteadyState, a search each, committed at the end of each sweep.

export Evolve, Thermalize, GroundState, SteadyState

"""
    step_number(span, step, what, sim)

the number of steps of `step` covering `span`, the step adjusted to the nearest one that does so
and taken with the sign of `span`, which is logged, as a phase of no step
"""
function step_number(span::Number, step::Number, what::String, sim::Simulation)
    n = round(Int, abs(span / step))
    if n == 0
        log_message(sim, "Skipping $what of $span, shorter than half a step")
    elseif !(span / n ≈ step)
        log_message(sim, "Taking a step of $(span / n) rather than $step to cover $span in $n steps")
    end
    return n
end

"""
    Evolve(; evolver, duration, time_step, algo = Tdvp(), limits, measurements,
             measurements_period, name, time_start, final_measurements)

a phase of time evolution, see `evolve`.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `evolver`: `-im * H` for a Hamiltonian `H`, plus dissipators for a mixed state, or for an
  evolver depending on time `[A, B] => [f, g]`, or `A => f`, see `evolve`
- `duration`: the duration of the evolution
- `time_step`: the time step, adjusted to the nearest one dividing the duration into a whole
  number of steps, and taken with the sign of the duration (no step when that number is zero)
- `algo`: the method, see `EvolutionAlgorithm`
- `limits`: the limits of the state, see `Limits`
- `measurements`: the measurements made after every `measurements_period` steps (default 1),
  see `output`; the state the phase starts from is measured by the final measurements of the
  phase before

# Examples

    Evolve(evolver = -im * (Z(1) * Z(2) + Z(2) * Z(3)), duration = 2., time_step = 0.1,
           measurements = "data" => [X, Y, Z])
"""
@kwdef struct Evolve <: AbstractPhase
    name::String = "Time evolution"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    evolver
    duration::Number
    time_step::Number
    algo::EvolutionAlgorithm = Tdvp()
    limits::Limits = Limits()
    measurements = []
    measurements_period::Int = 1
    function Evolve(name, time_start, final_measurements, evolver, duration, time_step, algo,
                    limits, measurements, measurements_period)
        if iszero(time_step)
            error("the time step of an Evolve cannot be zero")
        end
        # A => f written as the vectors, A => [f] left to be refused
        if evolver isa Pair && !(first(evolver) isa Union{AbstractVector, PreMPO}) &&
           !(last(evolver) isa AbstractVector)
            evolver = [first(evolver)] => [last(evolver)]
        end
        return new(name, time_start, final_measurements, evolver, duration, time_step, algo,
                   limits, measurements, measurements_period)
    end
end

"""
    EvolveSteps(algo, pre, coefs, limits, nsteps, start, duration)

an `Evolve` prepared: its method, the `PreMPO` of its evolver and its time functions, its
limits, its number of steps, the time it starts from and its duration
"""
struct EvolveSteps{A <: EvolutionAlgorithm, P <: Union{Nothing, PreMPO}, C}
    algo::A
    pre::P
    coefs::C
    limits::Limits
    nsteps::Int
    start::Number
    duration::Number
end

function prepare_phase(phase::Evolve, sim::Simulation)
    n = step_number(phase.duration, phase.time_step, "an evolution", sim)
    if n == 0
        return EvolveSteps(phase.algo, nothing, nothing, phase.limits, 0, sim.time, phase.duration)
    end
    log_message(sim, "Evolving the state from time $(sim.time) to $(sim.time + phase.duration)")
    pre, coefs = prepared_evolver(sim.state, phase.evolver)
    return EvolveSteps(phase.algo, pre, coefs, phase.limits, n, sim.time, phase.duration)
end

step_count(e::EvolveSteps) = e.nsteps

# the times of the steps are those of evolve, the last one the end asked for
function run_step(e::EvolveSteps, sim::Simulation, k::Int)
    dt = e.duration / e.nsteps
    st = evolve_step(e.algo, e.pre, e.coefs, sim.state, e.start + (k - 1) * dt, dt, k,
                     step_limits(e.limits, k))
    return Simulation(sim, st, k == e.nsteps ? e.start + e.duration : e.start + k * dt)
end

"""
    Thermalize(; hamiltonian, beta, beta_step, algo = Tdvp(), limits, measurements,
                 measurements_period, name, time_start, final_measurements)

a phase that takes a mixed state ``\\rho`` to ``e^{-\\beta H/2} \\rho \\, e^{-\\beta H/2}``,
normalized to trace one, by an evolution in imaginary time, see `thermalize`: from
`"FullyMixed"`, the thermal state ``e^{-\\beta H}/Z``. The time does not move. It provides the
values `:beta`, the inverse temperature reached, and `:log_trace`, the logarithm of the trace
before normalization, see `thermalize`.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `hamiltonian`: the Hamiltonian ``H``
- `beta`: the inverse temperature ``\\beta``
- `beta_step`: its step, adjusted as the time step of `Evolve`
- `algo`: the method, see `EvolutionAlgorithm`
- `limits`: the limits of the state, see `Limits`
- `measurements`: the measurements made after every `measurements_period` steps (default 1)

# Examples

    Thermalize(hamiltonian = -sum(Z(i) * Z(i + 1) for i in 1:9) - sum(X(i) for i in 1:10),
               beta = 2., beta_step = 0.1, limits = Limits(maxdim = 64),
               measurements = "thermo" => [:beta, :log_trace, Z(5)])
"""
@kwdef struct Thermalize <: AbstractPhase
    name::String = "Thermalization"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    hamiltonian
    beta::Real
    beta_step::Real
    algo::EvolutionAlgorithm = Tdvp()
    limits::Limits = Limits()
    measurements = []
    measurements_period::Int = 1
    function Thermalize(name, time_start, final_measurements, hamiltonian, beta, beta_step, algo,
                        limits, measurements, measurements_period)
        if iszero(beta_step)
            error("the beta_step of a Thermalize cannot be zero")
        end
        return new(name, time_start, final_measurements, hamiltonian, beta, beta_step, algo,
                   limits, measurements, measurements_period)
    end
end

provided_values(::Thermalize) = (:beta, :log_trace)

"""
    ThermalSteps(algo, pre, limits, nsteps, beta)

a `Thermalize` prepared: its method, the `PreMPO` of its Hamiltonian on the mixed state, its
limits, its number of steps and the inverse temperature it reaches
"""
struct ThermalSteps{A <: EvolutionAlgorithm, P <: Union{Nothing, PreMPO}}
    algo::A
    pre::P
    limits::Limits
    nsteps::Int
    beta::Real
end

function prepare_phase(phase::Thermalize, sim::Simulation)
    if !(sim.state isa State{Mixed})
        error("a thermal state needs a mixed representation: mix the state first, see ToMixed")
    end
    n = step_number(phase.beta, phase.beta_step, "a thermalization", sim)
    if n == 0
        return ThermalSteps(phase.algo, nothing, phase.limits, 0, phase.beta)
    end
    log_message(sim, "Thermalizing the state to beta $(phase.beta)")
    return ThermalSteps(phase.algo, prepared_operator(sim.state, phase.hamiltonian), phase.limits,
                        n, phase.beta)
end

step_count(t::ThermalSteps) = t.nsteps

# as thermalize: normalized at the first step only, the logarithms of the traces added up
function run_step(t::ThermalSteps, sim::Simulation, k::Int)
    state = k == 1 ? normalize(sim.state) : sim.state
    st, l = thermal_step(t.algo, t.pre, state, t.beta / t.nsteps, k, step_limits(t.limits, k))
    return Simulation(sim, st; beta = k == t.nsteps ? t.beta : k * (t.beta / t.nsteps),
                      log_trace = (k == 1 ? 0. : sim.log_trace) + l)
end

"""
    run_search(phase, sim, done, mpo, result; args...)

run the search of `phase`, `GroundState` or `SteadyState`, of the lowest eigenvector of `mpo`
from the state of `sim`, its first `done` sweeps done, each sweep ending with `step_done` on
the state `result` makes of the eigenvector, and return the simulation with that state and the
value `:energy`; `args` are the other keywords of `search`
"""
function run_search(phase, sim::Simulation, done::Int, mpo::MPO, result; args...)
    energy = done == 0 ? nothing : sim.energy
    reached = Ref(sim)
    hook(k, e, st, ends) =
        reached[] = step_done(Simulation(sim, result(st); energy = e), phase, k; last = ends)
    e, st = search(phase.algo, mpo, sim.state, phase.limits; sweep_done = hook, done, energy,
                   args...)
    log_message(sim, "Done, the energy is $e")
    return Simulation(reached[], result(st); energy = e)
end

"""
    GroundState(; hamiltonian, algo, limits, orthogonal_to, weight, measurements,
                  measurements_period, name, time_start, final_measurements)

a phase that searches the ground state of a Hamiltonian from the pure state, see
`ground_state`, or an excited state, the value `:energy` taking the energy of each sweep.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `hamiltonian`: the Hamiltonian
- `algo`: the method, see `SearchAlgorithm`
- `limits`: the limits of the sweeps, see `Limits`, which may give one value per sweep
- `orthogonal_to`: the eigenstates the state searched is orthogonal to, by their names in the
  simulation, see `CopyState`, or given, and `weight`, see `ground_state` (default none)
- `measurements`: the measurements made after every `measurements_period` sweeps (default 1)

# Examples

    GroundState(hamiltonian = X(1) * X(2), algo = Dmrg(nsweeps = 10, tol = 1e-10),
                limits = Limits(maxdim = [10, 20, 30]))
    [GroundState(hamiltonian = H, algo = Dmrg(nsweeps = 10)), CopyState(:state => :ground),
     GroundState(hamiltonian = H, algo = Dmrg(nsweeps = 20), orthogonal_to = [:ground])]
"""
@kwdef struct GroundState <: AbstractPhase
    name::String = "Ground state search"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    hamiltonian
    algo::SearchAlgorithm = nothing
    limits::Limits = Limits()
    orthogonal_to::Vector = []
    weight::Union{Nothing, Real} = nothing
    measurements = []
    measurements_period::Int = 1
    # the method, which has no default, asked for with what to write
    GroundState(name, time_start, final_measurements, hamiltonian, algo, limits, orthogonal_to,
                weight, measurements, measurements_period) =
        new(name, time_start, final_measurements, hamiltonian, search_algo(algo, "GroundState"),
            limits, orthogonal_to, weight, measurements, measurements_period)
end

provided_values(::GroundState) = (:energy,)

function run_phase(phase::GroundState, sim::Simulation, done::Int)
    if !(sim.state isa State{Pure})
        error("GroundState searches a pure state: search it pure and mix it")
    end
    pre = prepared_operator(sim.state, phase.hamiltonian)
    ms = [ m isa Symbol ? state_named(sim, m) : m for m in phase.orthogonal_to ]
    args = orthogonal_args(sim.state, pre, ms, phase.weight)
    log_message(sim, isempty(args) ? "Searching the ground state" :
                     "Searching an excited state, of weight $(args.weight)")
    sim = run_search(phase, sim, done, cached_mpos(() -> make_mpo(pre), pre, :mpo, nothing),
                     identity; args...)
    check_orthogonal(sim.state, args)
    return sim
end

"""
    SteadyState(; lindbladian, algo, limits, mpo_limits, mpo_algo, measurements,
                  measurements_period, name, time_start, final_measurements)

a phase that searches the steady state of a Lindbladian from the mixed state, see
`steady_state`, the value `:energy` taking the eigenvalue of ``L^\\dagger L`` of each sweep, 0
for a steady state, and the state measured after a sweep having trace one.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `lindbladian`: the Lindbladian ``L``, `-im * H` plus dissipators
- `algo`: the method, see `SearchAlgorithm`
- `limits`: the limits of the sweeps, see `Limits`
- `mpo_limits`, `mpo_algo`: the truncation and the algorithm of the MPO of ``L^\\dagger L``,
  see `steady_state`
- `measurements`: the measurements made after every `measurements_period` sweeps (default 1)

# Examples

    SteadyState(lindbladian = -im * H + D, algo = Dmrg(nsweeps = 10, tol = 1e-5),
                limits = Limits(maxdim = [10, 20, 50]))
"""
@kwdef struct SteadyState <: AbstractPhase
    name::String = "Steady state search"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    lindbladian
    algo::SearchAlgorithm = nothing
    limits::Limits = Limits()
    mpo_limits::Limits = Limits()
    mpo_algo::String = "naive"
    measurements = []
    measurements_period::Int = 1
    # refused when the phase is written, as the fields of Evolve and ApproxW, the method asked for
    # with what to write
    SteadyState(name, time_start, final_measurements, lindbladian, algo, limits, mpo_limits,
                mpo_algo, measurements, measurements_period) =
        new(name, time_start, final_measurements, lindbladian, search_algo(algo, "SteadyState"),
            limits, mpo_limits, check_mpo_algo(mpo_algo), measurements, measurements_period)
end

provided_values(::SteadyState) = (:energy,)

function run_phase(phase::SteadyState, sim::Simulation, done::Int)
    if !(sim.state isa State{Mixed})
        error("SteadyState searches a density matrix: mix the state first, see ToMixed")
    end
    l2 = squared_lindbladian(prepared_operator(sim.state, phase.lindbladian), phase.mpo_limits,
                             phase.mpo_algo)
    log_message(sim, "Searching the steady state")
    sim = run_search(phase, sim, done, l2, st -> st / trace(st); krylov = steady_krylov)
    steady_result(sim.energy, sim.state)
    return sim
end
