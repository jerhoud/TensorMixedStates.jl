# The samples of a simulation: Realizations, the last phase, whose function draws the disorder of a
# sample and returns its phases, and the loop running each sample from the branch point, the state
# the phases before Realizations reached, computed once.

export Realizations

"""
    Realizations(f; name = "Realizations")

the last phase of a simulation with `samples`, see `SimData`, where its samples branch off: each
sample starts from the state the phases before it reached, computed once, and runs the phases
`f(c)` returns, `c.sample` being the index of the sample and `c.rng` the generator `f` draws the
disorder of that sample with. The final measurements of the `SimData` are made at the end of each
sample, and the value `:sample` gives its index to the measurements.

`f` must return the same phases, by type and name, for every sample, and draw with `c.rng` alone,
so that a resumed sample gets the same phases. It is called once when the `SimData` is written,
with a generator of its own, so that its phases are checked before the run. The phase has no
`time_start` nor `final_measurements` of its own.

# Examples

    Realizations() do c
        J = [1 + 0.1 * randn(c.rng) for i in 1:9]
        H = -sum(J[i] * Z(i) * Z(i + 1) for i in 1:9) - sum(X(i) for i in 1:10)
        GroundState(hamiltonian = H, algo = Dmrg(nsweeps = 10),
                    final_measurements = "energies" => [:sample, :energy])
    end
"""
struct Realizations <: AbstractPhase
    name::String
    time_start::Nothing
    final_measurements::Vector{Any}
    realize::Function
end

Realizations(realize::Function; name::AbstractString = "Realizations") =
    Realizations(name, nothing, Any[], realize)

"""
    Realization(sample, rng)

what the function of `Realizations` is given: `c.sample`, the index of the sample, and `c.rng`,
the generator it draws the disorder of that sample with
"""
struct Realization
    sample::Int
    rng::Xoshiro
end

function getproperty(c::Realization, name::Symbol)
    if name ∉ (:sample, :rng)
        error("the function of Realizations is given c.sample and c.rng, not c.$name")
    end
    return getfield(c, name)
end

setproperty!(::Realization, ::Symbol, _) =
    error("the function of Realizations returns the phases of a sample and changes nothing: an " *
          "Update among those phases changes the simulation")

"""
    disorder_rng(seed, sample)
    sample_rng(seed, sample)

the generators of the sample `sample` of a run of seed `seed`: the one `Realizations` draws its
disorder with, and the one its phases draw with, `sim.rng`, independent of each other and of
those of the other samples, so that a sample is made again alone
"""
disorder_rng(seed::Integer, sample::Int) = Xoshiro(UInt64[seed, sample])

# the third number, 1, leaves room for several trajectories of one disorder
sample_rng(seed::Integer, sample::Int) = Xoshiro(UInt64[seed, sample, 1])

"""
    realized_phases(p, c)

the phases the function of `p` returns for `c`, flattened and checked as phases, refusing
another `Realizations` and a phase writing a file of its own, which each sample would overwrite
"""
function realized_phases(p::Realizations, c::Realization)
    ps = map(check_phase, flatten_phases(p.realize(c)))
    for q in ps
        if q isa Realizations
            error("the function of Realizations returns a Realizations, which only the last " *
                  "phase of a simulation may be")
        elseif !isempty(phase_files(q))
            error("$(type_name(typeof(q))) writes $(first(phase_files(q))), which every sample " *
                  "would overwrite: it cannot be among the phases of Realizations")
        end
    end
    return ps
end

"""
    branched_phases(phases, samples)

the phases the function of the last phase of `phases`, a `Realizations`, returns for the first
sample with a generator of its own, none without samples; a `Realizations` that is not the last
phase is refused, and so are samples without one and one without samples
"""
function branched_phases(phases::Vector, samples::Union{Nothing, Int})
    for (i, p) in enumerate(phases)
        if p isa Realizations && i < length(phases)
            error("Realizations is the last phase of a simulation, the phases of a sample being " *
                  "those its function returns")
        end
    end
    branched = !isempty(phases) && last(phases) isa Realizations
    if branched && isnothing(samples)
        error("a simulation ending with Realizations needs samples: SimData(samples = N, …)")
    elseif !branched && !isnothing(samples)
        error("a simulation with samples ends with the phase Realizations, which draws each sample")
    end
    return branched ? realized_phases(last(phases), Realization(1, Xoshiro(0))) : Any[]
end

"""
    check_realization(phases, reference, sample)

refuse the phases of the sample `sample` unless they are, by type and name, those of
`reference`, which `Realizations` returned when the `SimData` was written
"""
function check_realization(phases::Vector, reference::Vector, sample::Int)
    if length(phases) ≠ length(reference)
        error("Realizations returns $(length(phases)) phases for sample $sample, and " *
              "$(length(reference)) for the first: every sample needs the same phases")
    end
    for (k, (p, q)) in enumerate(zip(phases, reference))
        if typeof(p) ≠ typeof(q) || p.name ≠ q.name
            error("phase $k of sample $sample is $(type_name(typeof(p))) \"$(p.name)\", and " *
                  "$(type_name(typeof(q))) \"$(q.name)\" for the first: every sample needs the " *
                  "same phases")
        end
    end
    return nothing
end

"""
    run_simulation(sim, sd, samples, sample, start, done, ended, threading)

run the phases of the simulation `sd` and its final measurements, from the sample `sample`, the
phase `start` and `done` steps, the phase itself done when `ended`, see `Commit`: once without
samples, otherwise the phases before the branch point once, see `run_samples`
"""
function run_simulation(sim::Simulation, sd, samples::Union{Nothing, Int}, sample::Int,
                        start::Int, done::Int, ended::Bool, threading)
    if isnothing(samples)
        sim = run_phases(sim, sd.phases, start, done, ended, threading)
        phase_output(sim, sd.final_measurements)
        return sim
    end
    return run_samples(sim, sd, samples, sample, start, done, ended, threading)
end

"""
    run_samples(sim, sd, samples, sample, start, done, ended, threading)

run the simulation `sd`, which ends with `Realizations`, up to its sample `samples`, from where
`run_simulation` says: the phases before `Realizations` once, the branch point they reach kept
in `branch.h5`, then each sample from it, with its index `:sample`, its generator, the phases
`Realizations` returns and the final measurements of `sd`, the log writing the phases and steps
of the first only and a line for each. A run resumed beyond `samples` ends the sample it resumes,
whose final measurements the outputs reopened have lost. A run with a directory ends with a
checkpoint, which a run with more samples goes on from.
"""
function run_samples(sim::Simulation, sd, samples::Int, sample::Int, start::Int, done::Int,
                     ended::Bool, threading)
    run = getfield(sim, :run)
    c = run.checkpointer
    prefix = sd.phases[1:end - 1]
    if sample > samples
        log_line(run, "Sample $sample reached already, beyond the $samples asked for: only that " *
                      "one is ended")
    end
    final = max(samples, sample)
    fresh = sample == 0
    if fresh
        run.sample = 0
        sim = run_phases(sim, prefix, start, done, ended, threading; final = false)
        values = committed_values(getfield(sim, :values))
        run.branch = Commit(0, length(prefix) + 1, 0, false, run.phase_time, sim.time, sim.state,
                            values, Dict{Destination, Any}())
        if !isempty(c.dir)
            write_branch(c.dir, run.branch, c.seed)
        end
        sample = 1
    else
        run.branch = read_branch(c.dir)
    end
    b = run.branch
    for k in sample:final
        run.sample = k
        run.quiet = k > 1
        resumed = !fresh
        if fresh
            values = merge(b.values, Dict{Symbol, Any}(:sample => k, :step => 0,
                                                        :rng => sample_rng(c.seed, k)))
            sim = Simulation(b.state, b.time, values, run)
            start, done, ended = length(prefix) + 1, 0, false
            commit!(sim, start, 0)
        end
        fresh = true
        began = time()
        realized = realized_phases(last(sd.phases), Realization(k, disorder_rng(c.seed, k)))
        check_realization(realized, sd.realized, k)
        phases = [prefix; realized]
        sim = run_phases(sim, phases, start, done, ended, threading; final = false)
        # before the final measurements, as without samples: a new run makes them again; not for
        # a sample resumed after its phases, whose checkpoint is on the disk, where a run asked to
        # stop would stop again and again
        if k == final && !(resumed && start > length(phases))
            commit!(sim, length(phases) + 1, 0)
            checkpoint_or_stop(run; force = true)
        end
        phase_output(sim, sd.final_measurements)
        log_line(run, "Sample $k of $final done, $(round(time() - began; digits = 3)) seconds" *
                      (sim.state isa State ? ", maxdim $(maxlinkdim(sim.state))" : ""))
    end
    run.quiet = false
    return sim
end
