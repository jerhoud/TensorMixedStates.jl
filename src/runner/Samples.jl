# The samples of a simulation: Realizations, the last phase, whose function draws the disorder of a
# sample and returns its phases, the generators of a sample and the checks of its phases; runTMS
# runs the samples (Run.jl).

export Realizations

"""
    Realizations(f; name = "Realizations")

the last phase of a simulation with `samples` drawing a disorder, see `SimData`, where its
samples branch off: each sample starts from the state the phases before it reached, computed
once, and runs the phases `f(c)` returns, `c.sample` being the index of the sample and `c.rng`
the generator `f` draws the disorder of that sample with. The final measurements of the
`SimData` are made at the end of each sample, and the value `:sample` gives its index to the
measurements.

`f` must return the same phases, by type and name, for every sample, and draw with `c.rng` alone,
so that a resumed sample gets the same phases. It is called once when the `SimData` is written,
with a generator of its own, so that its phases are checked before the run, and they are checked
again for each sample before it runs. The phase has no `time_start` nor `final_measurements` of
its own.

# Examples

    Realizations() do c
        J = [1 + 0.1 * randn(c.rng) for i in 1:9]
        H = -sum(J[i] * Z(i) * Z(i + 1) for i in 1:9) - sum(X(i) for i in 1:10)
        GroundState(hamiltonian = H, algo = Dmrg(nsweeps = 10),
                    final_measurements = EachSample("energies") => :energy)
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
    splitmix(x)

the finalizer of SplitMix64 (Steele, Lea and Flood 2014), a bijection of `UInt64` mixing its bits
"""
function splitmix(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

"""
    disorder_rng(seed, sample)
    sample_rng(seed, sample)

the generators of a sample, for its disorder and for its phases (`sim.rng`), independent of each
other and of those of the other samples
"""
disorder_rng(seed::Integer, sample::Int) = Xoshiro(mixed_seed(seed, sample, 0))

# the stream 1 leaves room for several trajectories of one disorder
sample_rng(seed::Integer, sample::Int) = Xoshiro(mixed_seed(seed, sample, 1))

# an integer seed, Julia having changed how it hashes a seed of several numbers after 1.10
mixed_seed(seed::Integer, sample::Int, stream::Int) =
    splitmix(splitmix(splitmix(UInt64(seed)) ⊻ UInt64(sample)) ⊻ UInt64(stream))

"""
    realized_phases(p, c)

the phases the function of `p` returns for `c`, flattened and checked as phases, another
`Realizations` refused
"""
function realized_phases(p::Realizations, c::Realization)
    ps = map(check_phase, flatten_phases(p.realize(c)))
    if any(q -> q isa Realizations, ps)
        error("the function of Realizations returns a Realizations, which only the last phase of " *
              "a simulation may be")
    end
    return ps
end

"""
    trajectory_branch(phases, samples)

`phases`, those from the first that makes a trajectory run as the phases of each sample, as
`Realizations` runs them, when no `Realizations` draws the samples; a trajectory made without
samples, or before `Realizations`, refused
"""
function trajectory_branch(phases::Vector, samples::Union{Nothing, Int})
    i = findfirst(makes_trajectory, phases)
    j = findfirst(p -> p isa Realizations, phases)
    if isnothing(i)
        return phases
    elseif !isnothing(j) && i < j
        error("$(type_name(typeof(phases[i]))) \"$(phases[i].name)\" comes before Realizations, " *
              "so that every sample would go on from one trajectory: put it among the phases " *
              "Realizations returns")
    elseif !isnothing(j)
        return phases
    elseif isnothing(samples)
        error("a simulation making trajectories needs samples: SimData(samples = N, …), N = 1 " *
              "for a single trajectory")
    end
    sample = phases[i:end]
    return [phases[1:i - 1]; Realizations(_ -> sample; name = "Trajectories")]
end

"""
    branched_phases(phases, samples)

the phases of the dry call of the `Realizations` ending `phases`, none without samples; a
misplaced `Realizations`, samples without one and one without samples refused
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
        error("a simulation with samples ends with the phase Realizations, which draws each " *
              "sample, or makes trajectories, with CreateState{Trajectory} or " *
              "ConvertState{Trajectory}")
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
