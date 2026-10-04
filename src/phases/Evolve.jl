# The Evolve phase, which evolves the state in time, and the algorithms it takes, Tdvp and
# ApproxW, with the methods of evolve that run them.

export Algo, Tdvp, ApproxW, Evolve

"""
    abstract type Algo

the supertype of the time evolution algorithms the `algo` field of `Evolve` takes: `Tdvp` and
`ApproxW`, and those an extension defines with a method of `TensorMixedStates.evolve`.
"""
abstract type Algo end

"""
    Tdvp(; expand_period = 0, hermitianize_period = 0, krylov = Krylov())

the tdvp algorithm, for the `algo` field of `Evolve`, see `tdvp`.

- `expand_period`: enlarge the bond dimension of the state by a global Krylov expansion before the
  first step and then every `expand_period` steps (default 0, never)
- `hermitianize_period`: make a mixed state hermitian every `hermitianize_period` steps (default 0,
  never)
- `krylov`: the parameters of the Krylov exponentiation of each local step, see `Krylov`
  (default `Krylov()`)

# Examples

    Tdvp()
    Tdvp(expand_period = 5)               # tdvp with expansion steps every 5 steps
    Tdvp(hermitianize_period = 3)         # tdvp, make hermitian every 3 steps
    Tdvp(krylov = Krylov(tol = 1e-10))    # tdvp, local steps at a lower precision
"""
@kwdef struct Tdvp <: Algo
    expand_period::Int = 0
    hermitianize_period::Int = 0
    krylov::Krylov = Krylov()
end

"""
    ApproxW(; order, w = 2, hermitianize_period = 0, apply_algo = "densitymatrix")

time evolution by WI or WII approximations of the exponential, combined into an approximation
of the given order, for the `algo` field of `Evolve`, see `approx_W`.

- `order`: the order of the approximation, from 1 to 4, required
- `w`: 1 or 2 for WI or WII (default 2)
- `hermitianize_period`: make a mixed state hermitian every `hermitianize_period` steps (default 0,
  never)
- `apply_algo`: the algorithm of the product of the state by each MPO, `"densitymatrix"`
  (default) or `"naive"`, see `approx_W`

# Examples

    ApproxW(order = 2)                            # order 2, WII
    ApproxW(order = 4, w = 1)                     # order 4, WI
    ApproxW(order = 4, hermitianize_period = 3)   # order 4, make hermitian every 3 steps
    ApproxW(order = 2, apply_algo = "naive")      # order 2, naive products
"""
@kwdef struct ApproxW <: Algo
    order::Int
    w::Int = 2
    hermitianize_period::Int = 0
    apply_algo::String = "densitymatrix"
    # checked when the phase is written rather than when it runs: corrected then, it could not
    # resume its checkpoint, which belongs to a simulation of other phases
    function ApproxW(order, w, hermitianize_period, apply_algo)
        check_w_approx(order, w)
        check_apply_algo(apply_algo)
        return new(order, w, hermitianize_period, apply_algo)
    end
end

# an algorithm is printed on one line, field by field, read from its type as a phase is, so
# that a field added to it shows up in the log and in `prog.jl` with nothing else to change
show(io::IO, s::Algo) =
    print(io, nameof(typeof(s)), "(",
          join(("$f = $(repr(getfield(s, f)))" for f in fieldnames(typeof(s))), ", "), ")")

"""
    Evolve(; duration, time_step, algo, evolver, measurements, measurements_period, limits,
             options...)

a phase of time evolution.

# Fields

- `name`, `time_start`, `final_measurements`: the fields every phase has, see `AbstractPhase`
- `limits`: constraints on the state, see `Limits` (default `Limits()`)
- `duration`: the duration of the evolution
- `time_step`: the time step, adjusted to the nearest one that divides the duration into a
  whole number of steps, and taken with the sign of the duration (the phase is skipped when
  that number is zero)
- `algo`: the algorithm, `Tdvp(...)` or `ApproxW(...)`
- `evolver`: `-im * H` for a hamiltonian `H`, plus dissipators for a mixed state, or
  `evolvers => coefs` for a time dependent one, see the `coefs` option of `tdvp`
- `measurements`: the measurements to make during the evolution, see `output` (default `[]`),
  after every `measurements_period` time steps. The state the phase starts from is not measured
  here: the `final_measurements` of the phase before measure it
- `measurements_period`: the number of time steps between two measurements (default 1)

# Examples

    Evolve(duration = 2., time_step = 0.1, algo = Tdvp(),
           evolver = -im * (Z(1)Z(2) + Z(2)Z(3)), measurements = "data" => [X, Y, Z])
"""
@kwdef struct Evolve <: AbstractPhase
    name::String = "Time evolution"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    limits::Limits = Limits()
    duration::Number
    time_step::Number
    algo::Algo
    evolver::Union{IndexedOp, Pair}
    measurements_period::Int = 1
    measurements = []
    # refused here rather than when it runs, as the fields of ApproxW
    function Evolve(name, time_start, final_measurements, limits, duration, time_step, algo,
                    evolver, measurements_period, measurements)
        if iszero(time_step)
            error("the time step of an Evolve cannot be zero")
        end
        # a single term and its function of time, written without the vectors, which failed on
        # a MethodError at the first step, and a function per term, which only the first step
        # checked
        if evolver isa Pair
            ops, fs = evolver
            if !(ops isa AbstractVector) && !(fs isa AbstractVector)
                evolver = [ops] => [fs]
            elseif !(ops isa AbstractVector && fs isa AbstractVector) || length(ops) ≠ length(fs)
                error("a time dependent evolver gives one function of time per term, as " *
                      "[A, B] => [f, g] or A => f")
            end
        end
        return new(name, time_start, final_measurements, limits, duration, time_step, algo, evolver,
                   measurements_period, measurements)
    end
end

"""
    evolve(algo, state, sim, phase; evolver, coefs, nsteps, kwargs...)

the simulation `sim` once the `Evolve` phase `phase` has evolved its state `state` with the
algorithm `algo`, in `nsteps` steps covering `phase.duration`: `evolver` is the evolver of the
phase, and `coefs` the functions of time of a time dependent one, or `nothing`. The state is
given apart from the simulation so that a method is chosen by its type as well as by that of
the algorithm.

An algorithm of an extension, or an algorithm for the state of an extension, comes with a
method of its own. It is called before the phase has read its resume point, so that a method
can read it with `resume_step`, or run its steps with `run_steps`, which resumes, stops and
checkpoints them; `output(sim, phase.measurements; sweep)` writes the measurements of the phase,
every `phase.measurements_period` steps. Such a method takes `kwargs...` after the keywords it
uses, so that a keyword a later version of TMS passes does not break it.
"""
function evolve(algo::Tdvp, state::State, sim::Simulation, phase::Evolve; evolver, coefs,
                nsteps)
    done, _ = resume_step(sim)
    st = tdvp(evolver, phase.duration, state; coefs, algo.hermitianize_period, nsweeps = nsteps,
              time_start = sim.time, phase.limits, first_sweep = done + 1, algo.expand_period,
              algo.krylov,
              observer! = TdvpObserver(sim, phase.measurements, phase.measurements_period))
    return Simulation(sim, st)
end

function evolve(algo::ApproxW, state::State, sim::Simulation, phase::Evolve; evolver, coefs,
                nsteps)
    done, _ = resume_step(sim)
    st = approx_W(evolver, phase.duration, state; coefs, algo.hermitianize_period, nsweeps = nsteps,
                  time_start = sim.time, phase.limits, first_sweep = done + 1, algo.order,
                  algo.w, algo.apply_algo,
                  observer! = ApproxWObserver(sim, phase.measurements, phase.measurements_period))
    return Simulation(sim, st)
end

evolve(algo, state, ::Simulation, ::Evolve; kwargs...) =
    error("TensorMixedStates.evolve has no method for $(typeof(algo)) on a $(typeof(state))")

function run_phase(sim::Simulation, phase::Evolve)
    # the duration gives the direction, and a step of the other sign is adjusted as one that
    # does not divide it: of the opposite sign, it made no step while the time went on
    nsteps = round(Int, abs(phase.duration / phase.time_step))
    if nsteps == 0
        log_message(sim, "Skipping an evolution of $(phase.duration), shorter than half a " *
                         "time step")
        return sim
    end
    # the step is adjusted rather than the duration, so that the phase ends where it was asked
    # to: a duration of 1 in steps of 0.3 stopped at 0.9
    duration = phase.duration
    if !(duration / nsteps ≈ phase.time_step)
        log_message(sim, "Taking a time step of $(duration / nsteps) rather than " *
                         "$(phase.time_step) to cover the duration $duration in $nsteps steps")
    end
    time_stop = sim.time + duration
    # a phase resumed in its course goes on from the time of its checkpoint, not from its start
    log_message(sim, "Evolving state from simulation time $(resume_time(sim)) to $(time_stop)")
    evolver, coefs = phase.evolver isa Pair ? (first(phase.evolver), last(phase.evolver)) :
                                              (phase.evolver, nothing)
    sim = evolve(phase.algo, sim.state, sim, phase; evolver, coefs, nsteps)
    # a phase cut short by a checkpoint stops at the time it actually reached
    return Simulation(sim, sim.state, stopped(sim) ? committed_time(sim) : time_stop)
end
