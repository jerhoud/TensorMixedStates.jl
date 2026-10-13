# The Update phase, which changes the simulation by a function of one's own, its state, its time
# and its values, at once or step by step.

export Update

"""
    Update(f; nsteps = 1, provides = [], measurements, measurements_period, name, time_start,
           final_measurements)

a phase that changes the simulation by the function `f(c)`, which reads it as a measurement
reads its context, `c.state`, `c.time`, `c.step`, `c.rng`, `c.ground`, and changes it by
assignments:

- `c.state = st`, `c.time = t`: the state and the time
- `c.x = v`: a value of the simulation, measured and kept by the checkpoints, `c.x = nothing`
  removing it, which frees the memory of a state

The phase makes `nsteps` steps, the run measuring `measurements` every `measurements_period`
steps, committing and stopping between two, as for `Evolve`, `c.step` being the steps done
before. A stopped run resumes after the last step committed: `f` must draw with `c.rng` alone
and act only through `c`, `output(c, …)` and `log_message(c, …)`, what it changes elsewhere
being neither kept nor done again. `f` works on a copy, which a failure leaves the simulation
without.

- `provides`: the names of the values `f` sets, which the measurements may then read, see
  `TensorMixedStates.provided_values`

# Examples

    Update(c -> c.ground = c.state; provides = [:ground])
    Update(c -> c.ground = nothing)
    Update(nsteps = 1000, provides = [:m], measurements = "kicks" => [Z(1), :m],
           measurements_period = 10) do c
        c.m, _, c.state = collapse(apply(kick, c.state), Z(1); rng = c.rng)
    end
"""
struct Update <: AbstractPhase
    name::String
    time_start::Union{Nothing, Number}
    final_measurements
    update::Function
    nsteps::Int
    provides::Vector{Symbol}
    measurements
    measurements_period::Int
end

function Update(update::Function; name::AbstractString = "Update", time_start = nothing,
                final_measurements = [], nsteps::Int = 1, provides = Symbol[], measurements = [],
                measurements_period::Int = 1)
    for n in provides
        if n in reserved_values
            error("an Update cannot provide $n, a name the run reserves")
        end
    end
    return Update(name, time_start, final_measurements, update, nsteps, collect(Symbol, provides),
                  measurements, measurements_period)
end

provided_values(phase::Update) = Tuple(phase.provides)

prepare_phase(phase::Update, ::Simulation) = phase

step_count(phase::Update) = phase.nsteps

function run_step(phase::Update, sim::Simulation, ::Int)
    c = Editing(sim)
    phase.update(c)
    return Simulation(c)
end
