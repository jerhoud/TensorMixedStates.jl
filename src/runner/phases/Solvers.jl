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
        log_message(sim, "Taking a step of $(span / n) rather than $step to cover $span in " *
                         "$n steps")
    end
    return n
end

"""
    Evolve(; evolver, duration, time_step | nsteps, algo = Tdvp(), limits, mpo_limits,
             measurements, measurements_period, on, name, time_start, final_measurements)

a phase of time evolution, see `evolve`.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `evolver`: `-im * H` for a Hamiltonian `H`, plus dissipators for a mixed state, or for an
  evolver depending on time `[A, B] => [f, g]`, see `evolve`, or `A => f` for a single term
- `duration`: the duration of the evolution
- `time_step`: the time step, adjusted to the nearest one dividing the duration into a whole
  number of steps, and taken with the sign of the duration (no step when that number is zero),
  or `nsteps`, the number of steps, one of the two being given
- `algo`: the method, see `EvolutionAlgorithm`
- `limits`: the limits of the state, see `Limits`
- `mpo_limits`: the truncation of the MPO of the evolver, see `PreMPO` (default none)
- `measurements`: the measurements made after every `measurements_period` steps (default 1),
  see `output`; the state the phase starts from is measured by the final measurements of the
  phase before
- `on`: the name of the state the phase evolves, see `CopyState`, or a vector of names, the
  states of one system and one representation evolved together under the evolver, the time
  advancing once, for two time correlations for instance (default `:state`)

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
    time_step::Union{Nothing, Number} = nothing
    nsteps::Union{Nothing, Int} = nothing
    algo::EvolutionAlgorithm = Tdvp()
    limits::Limits = Limits()
    mpo_limits::Limits = Limits()
    measurements = []
    measurements_period::Int = 1
    on::Union{Symbol, Vector{Symbol}} = :state
    function Evolve(name, time_start, final_measurements, evolver, duration, time_step, nsteps,
                    algo, limits, mpo_limits, measurements, measurements_period, on)
        evolution_steps("an Evolve", duration, nsteps, time_step, "time_step")
        # A => f written as the vectors, A => [f] left to be refused
        if evolver isa Pair && !(first(evolver) isa Union{AbstractVector, PreMPO}) &&
           !(last(evolver) isa AbstractVector)
            evolver = [first(evolver)] => [last(evolver)]
        end
        return new(name, time_start, final_measurements, evolver, duration, time_step, nsteps,
                   algo, limits, mpo_limits, measurements, measurements_period, on)
    end
end

"""
    EvolveSteps(algo, pre, coefs, limits, nsteps, start, duration, on)

an `Evolve` prepared: its method, the `PreMPO` of its evolver and its time functions, its
limits, its number of steps, the time it starts from, its duration and the names of the states
it evolves
"""
struct EvolveSteps{A <: EvolutionAlgorithm, P <: Union{Nothing, PreMPO}, C}
    algo::A
    pre::P
    coefs::C
    limits::Limits
    nsteps::Int
    start::Number
    duration::Number
    on::Vector{Symbol}
end

function prepare_phase(phase::Evolve, sim::Simulation)
    names = state_names(phase.on)
    n = if isnothing(phase.nsteps)
        step_number(phase.duration, phase.time_step, "an evolution", sim)
    else
        phase.nsteps
    end
    if n == 0
        return EvolveSteps(phase.algo, nothing, nothing, phase.limits, 0, sim.time, phase.duration,
                           names)
    end
    states = [ state_named(sim, name) for name in names ]
    # one PreMPO for all, which is that of a representation and of a system
    if any(st -> typeof(st) ≠ typeof(states[1]) || st.system ≠ states[1].system, states)
        error("the states an Evolve evolves together need one system and one representation")
    end
    what = names == [:state] ? "the state" : "the states " * join(names, ", ")
    log_message(sim, "Evolving $what from time $(sim.time) to $(sim.time + phase.duration)")
    pre, coefs = prepared_evolver(states[1], phase.evolver, phase.mpo_limits)
    return EvolveSteps(phase.algo, pre, coefs, phase.limits, n, sim.time, phase.duration, names)
end

step_count(e::EvolveSteps) = e.nsteps

# the times of the steps are those of evolve, the last one the end asked for
function run_step(e::EvolveSteps, sim::Simulation, k::Int)
    dt = e.duration / e.nsteps
    for name in e.on
        st = evolve_step(e.algo, e.pre, e.coefs, state_named(sim, name), e.start + (k - 1) * dt,
                         dt, k, step_limits(e.limits, k))
        sim = with_state(sim, name, st)
    end
    return Simulation(sim, sim.state, k == e.nsteps ? e.start + e.duration : e.start + k * dt)
end

"""
    Thermalize(; hamiltonian, beta, beta_step | nsteps, algo = Tdvp(), limits, mpo_limits,
                 measurements, measurements_period, on, name, time_start, final_measurements)

a phase that takes a mixed state ``\\rho`` to ``e^{-\\beta H/2} \\rho \\, e^{-\\beta H/2}``,
normalized to trace one, by an evolution in imaginary time, see `thermalize`: from
`"FullyMixed"`, the thermal state ``e^{-\\beta H}/Z``. The time does not move. It provides the
values `:beta`, the inverse temperature reached, and `:log_trace`, the logarithm of the trace
before normalization, see `thermalize`: 0 both for a `beta` shorter than half a step, which
leaves the state normalized only.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `hamiltonian`: the Hamiltonian ``H``
- `beta`: the inverse temperature ``\\beta``
- `beta_step`: its step, adjusted as the time step of `Evolve`, or `nsteps`, the number of steps,
  one of the two being given
- `algo`: the method, see `EvolutionAlgorithm`
- `limits`: the limits of the state, see `Limits`
- `mpo_limits`: the truncation of the MPO of the Hamiltonian, see `PreMPO` (default none)
- `measurements`: the measurements made after every `measurements_period` steps (default 1)
- `on`: the name of the state the phase acts on, see `CopyState` (default `:state`)

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
    beta_step::Union{Nothing, Real} = nothing
    nsteps::Union{Nothing, Int} = nothing
    algo::EvolutionAlgorithm = Tdvp()
    limits::Limits = Limits()
    mpo_limits::Limits = Limits()
    measurements = []
    measurements_period::Int = 1
    on::Symbol = :state
    function Thermalize(name, time_start, final_measurements, hamiltonian, beta, beta_step,
                        nsteps, algo, limits, mpo_limits, measurements, measurements_period, on)
        evolution_steps("a Thermalize", beta, nsteps, beta_step, "beta_step")
        return new(name, time_start, final_measurements, hamiltonian, beta, beta_step, nsteps,
                   algo, limits, mpo_limits, measurements, measurements_period, on)
    end
end

provided_values(::Thermalize) = (:beta, :log_trace)

"""
    ThermalSteps(algo, pre, limits, nsteps, beta, on)

a `Thermalize` prepared: its method, the `PreMPO` of its Hamiltonian on the mixed state,
`nothing` for one shorter than half a step, made of a step that normalizes only, its limits, its
number of steps, the inverse temperature it reaches and the name of its state
"""
struct ThermalSteps{A <: EvolutionAlgorithm, P <: Union{Nothing, PreMPO}}
    algo::A
    pre::P
    limits::Limits
    nsteps::Int
    beta::Real
    on::Symbol
end

function prepare_phase(phase::Thermalize, sim::Simulation)
    state = state_named(sim, phase.on)
    if !(state isa State{Mixed})
        error("a thermal state needs a mixed representation: mix the state first, see ToMixed")
    end
    n = if isnothing(phase.nsteps)
        step_number(phase.beta, phase.beta_step, "a thermalization", sim)
    else
        phase.nsteps
    end
    # a step that only normalizes, giving the values all the same
    if n == 0
        return ThermalSteps(phase.algo, nothing, phase.limits, 1, 0., phase.on)
    end
    log_message(sim, "Thermalizing the state to beta $(phase.beta)")
    pre = prepared_operator(state, phase.hamiltonian, phase.mpo_limits)
    return ThermalSteps(phase.algo, pre, phase.limits, n, phase.beta, phase.on)
end

step_count(t::ThermalSteps) = t.nsteps

# as thermalize: normalized at the first step only, the logarithms of the traces added up
function run_step(t::ThermalSteps, sim::Simulation, k::Int)
    state = k == 1 ? normalize(state_named(sim, t.on)) : state_named(sim, t.on)
    if isnothing(t.pre)
        return with_state(sim, t.on, state; beta = 0., log_trace = 0.)
    end
    st, l = thermal_step(t.algo, t.pre, state, t.beta / t.nsteps, k, step_limits(t.limits, k))
    return with_state(sim, t.on, st; beta = k == t.nsteps ? t.beta : k * (t.beta / t.nsteps),
                      log_trace = (k == 1 ? 0. : sim.log_trace) + l)
end

"""
    run_search(phase, sim, done, mpo, result, check; args...)

run the search of `phase`, `GroundState` or `SteadyState`, of the lowest eigenvector of `mpo`
from the state `phase.on` of `sim`, its first `done` sweeps done, each sweep ending with
`step_done` on the state `result` makes of the eigenvector, and return the simulation with that
state and the eigenvalue as the value `name`; `check(value, state)` warns about the state of the
last sweep, before its step, and `args` are the other keywords of `search`
"""
function run_search(phase, sim::Simulation, done::Int, mpo::MPO, result, check, name::Symbol;
                    args...)
    energy = done == 0 ? nothing : getproperty(sim, name)
    reached = Ref(sim)
    function hook(k, e, st, ends)
        r = result(st)
        # before the step, which may stop the run: a resumed run skips the phase ended there
        if ends
            check(e, r)
        end
        reached[] = step_done(with_state(sim, phase.on, r; (name => e,)...), phase, k; last = ends)
    end
    e, st = search(phase.algo, mpo, state_named(sim, phase.on), phase.limits; sweep_done = hook,
                   done, energy, args...)
    return with_state(reached[], phase.on, result(st); (name => e,)...)
end

"""
    GroundState(; hamiltonian, algo, limits, mpo_limits, orthogonal_to, weight, measurements,
                  measurements_period, on, name, time_start, final_measurements)

a phase that searches the ground state of a Hamiltonian from the pure state, see
`ground_state`, or an excited state, the value `:energy` taking the energy of each sweep.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `hamiltonian`: the Hamiltonian
- `algo`: the method, see `SearchAlgorithm`
- `limits`: the limits of the sweeps, see `Limits`, which may give one value per sweep
- `mpo_limits`: the truncation of the MPO of the Hamiltonian, see `PreMPO` (default none)
- `orthogonal_to`: the eigenstates the state searched is orthogonal to, by their names in the
  simulation, see `CopyState`, or given, and `weight`, see `ground_state` (default none)
- `measurements`: the measurements made after every `measurements_period` sweeps (default 1)
- `on`: the name of the state the phase acts on, see `CopyState` (default `:state`)

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
    mpo_limits::Limits = Limits()
    orthogonal_to::Vector = []
    weight::Union{Nothing, Real} = nothing
    measurements = []
    measurements_period::Int = 1
    on::Symbol = :state
    # the method, which has no default, asked for with what to write
    GroundState(name, time_start, final_measurements, hamiltonian, algo, limits, mpo_limits,
                orthogonal_to, weight, measurements, measurements_period, on) =
        new(name, time_start, final_measurements, hamiltonian, search_algo(algo, "GroundState"),
            limits, mpo_limits, orthogonal_to, weight, measurements, measurements_period, on)
end

provided_values(::GroundState) = (:energy,)

function run_phase(phase::GroundState, sim::Simulation, done::Int)
    state = state_named(sim, phase.on)
    if !(state isa State{Pure})
        error("GroundState searches a pure state: search it pure and mix it")
    end
    pre = prepared_operator(state, phase.hamiltonian, phase.mpo_limits)
    ms = [ m isa Symbol ? state_named(sim, m) : m for m in phase.orthogonal_to ]
    args = orthogonal_args(state, pre, ms, phase.weight)
    log_message(sim, isempty(args) ? "Searching the ground state" :
                     "Searching an excited state, of weight $(args.weight)")
    return run_search(phase, sim, done, cached_mpos(() -> make_mpo(pre), pre, :mpo, nothing),
                      identity, (_, st) -> check_orthogonal(st, args, phase.weight), :energy;
                      args...)
end

"""
    SteadyState(; lindbladian, algo, limits, mpo_limits, measurements, measurements_period, on,
                  name, time_start, final_measurements)

a phase that searches the steady state of a Lindbladian from the mixed state, see
`steady_state`, the value `:residual` taking the eigenvalue of ``L^\\dagger L`` of each sweep, 0
for a steady state, and the state measured after a sweep having trace one. Its lists of
measurements cannot read `:energy`, which would be that of another phase.

- `name`, `time_start`, `final_measurements`: the fields of every phase, see
  `TensorMixedStates.AbstractPhase`
- `lindbladian`: the Lindbladian ``L``, `-im * H` plus dissipators
- `algo`: the method, see `SearchAlgorithm`
- `limits`: the limits of the sweeps, see `Limits`
- `mpo_limits`: the truncation of the MPO of ``L^\\dagger L``, see `PreMPO` (default none)
- `measurements`: the measurements made after every `measurements_period` sweeps (default 1)
- `on`: the name of the state the phase acts on, see `CopyState` (default `:state`)

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
    measurements = []
    measurements_period::Int = 1
    on::Symbol = :state
    # refused when the phase is written, as the fields of Evolve and ApproxW, the method asked for
    # with what to write
    SteadyState(name, time_start, final_measurements, lindbladian, algo, limits, mpo_limits,
                measurements, measurements_period, on) =
        new(name, time_start, final_measurements, lindbladian, search_algo(algo, "SteadyState"),
            limits, mpo_limits, measurements, measurements_period, on)
end

provided_values(::SteadyState) = (:residual,)

refused_reads(::SteadyState) =
    ((:energy, "the value of a steady state search is residual, the eigenvalue of L†L"),)

function run_phase(phase::SteadyState, sim::Simulation, done::Int)
    state = state_named(sim, phase.on)
    if !(state isa State{Mixed})
        error("SteadyState searches a density matrix: mix the state first, see ToMixed")
    end
    l2 = squared_lindbladian(prepared_operator(state, phase.lindbladian), phase.mpo_limits)
    log_message(sim, "Searching the steady state")
    return run_search(phase, sim, done, l2, st -> st / trace(st), steady_result, :residual;
                      krylov = steady_krylov)
end
