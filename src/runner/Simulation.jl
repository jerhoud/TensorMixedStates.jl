# Simulation, a state, its time and its values as runTMS runs it, output, which measures it into
# the destinations of the run, log_message, which writes to its log, and the logger sending the
# warnings of the package there.

export Simulation, output, log_message, stopped, Backend

"""
    simulation_files

the files `runTMS` writes in the directory of a simulation, which no destination may be named
after, nor be in its directory `prog`, which holds the files the program includes, nor, with
samples, be one of `sample_files` or in `parts` and `pending`, see `is_simulation_file`
"""
const simulation_files = Set(["log", "stop", "restart", "error", "running", "stamp", "description",
                              "lock", "prog.jl", "prog_args.json", checkpoint_files...])

"""
    is_simulation_file(name, sampled)

whether `name`, relative to the directory of a simulation, `sampled` or not, is a file `runTMS`
writes there, see `simulation_files`
"""
function is_simulation_file(name::AbstractString, sampled::Bool)
    n = normpath(name)
    under(d) = n == d || startswith(n, joinpath(d, ""))
    return n in simulation_files || under("prog") ||
           sampled && (n in sample_files || under("parts") || under("pending"))
end

"""
    reserved_values

the names no value of a simulation may take: those the measurements read otherwise, `c.state`
and `c.time`, those the run gives, `c.step` and `c.rng`, and `sim.data`; `c.sample` too in a
simulation with samples, see `gives_sample`
"""
const reserved_values = (:state, :time, :step, :rng, :data)

"""
    check_backend(backend)

`backend`, refused unless it is a type of states named without its representation, as
`ExactState`, or a `Backend`
"""
function check_backend(backend)
    t = backend isa Type ? Base.unwrap_unionall(backend) : nothing
    # the type as given, its parameters bound: Julia 1.14 no longer takes State{R}, R free, as a
    # subtype of AbstractState
    if !(t isa DataType && backend <: AbstractState && !isabstracttype(t))
        error("a backend is a type of states, as State or ExactState, or Backend(T; options...), " *
              "not $(repr(backend))")
    elseif backend isa DataType && backend.name.wrapper !== backend
        error("a backend is named without its representation: backend = $(nameof(backend))")
    end
    return backend
end

check_backend(::Nothing) = nothing

"""
    Backend(T; options...)

the backend `T`, a type of states as `ExactState` or `GaussianState`, whose states the phases create
from a `State` as `T(state; options...)`; `backend = T` takes its defaults. Converting a state already
of type `T` with options is refused. Experimental: it may change in a version 2.x.

# Examples

    SimData(name = "xy", backend = Backend(GaussianState; axes = :X), phases = [...])
"""
struct Backend{T}
    options::NamedTuple
    function Backend(T::Type; options...)
        check_backend(T)
        o = NamedTuple(options)
        if haskey(o, :limits)
            error("the limits of a conversion are given to the phase, not to Backend")
        elseif !isempty(o) && !any(R -> hasmethod(T, Tuple{State{R}}, keys(o)), (Pure, Mixed, Trajectory))
            error("$(nameof(T)) takes no option $(join(keys(o), ", ")) when created from a State, see its " *
                  "documentation")
        end
        return new{T}(o)
    end
end

check_backend(backend::Backend) = backend

function show(io::IO, b::Backend{T}) where T
    print(io, "Backend(", nameof(T))
    if !isempty(b.options)
        print(io, "; ", join(("$k = $(repr(v))" for (k, v) in pairs(b.options)), ", "))
    end
    print(io, ")")
end

# the backend as the log names it
backend_name(T::Type) = string(nameof(T))
backend_name(b::Backend) = sprint(show, b)

"""
    RunData(dir, outputs, log, redirect, checkpointer, provided, backend)

what the simulations of a run share: its directory, empty for a run whose output is redirected
to the stream `redirect`, its outputs, its log, its checkpointing, the values its phases
provide, the type of states its phases create by default, see `SimData`, the phase being run
and the time it started from, when its step being run started, in seconds of `time()`, whether
it stopped, the warnings its log received, the imaginary parts of
Hermitian measurements dropped and not yet reported, how many and the largest, as
`(imaginary part, measurement, time)`, the lists of measurements of the phase being run made
into `Measure`s, see `phase_output`, the sample being run, 0 without one, the phases it runs,
whether the log leaves out the lines of its phases and steps, and the outputs the samples are
merged into in a run without a directory, while its outputs record those of a sample
"""
mutable struct RunData
    dir::String
    outputs::Outputs
    log::IO
    redirect::Union{Nothing, IO}
    checkpointer::Checkpointer
    provided::Vector{Symbol}
    backend::Union{Type, Backend}
    phase::Int
    phase_time::Number
    step_start::Float64
    stopped::Bool
    warnings::Int
    dropped::Int
    largest::Union{Nothing, Tuple}
    measures::IdDict{Any, Vector}
    sample::Int
    phases::Vector{Any}
    quiet::Bool
    merged::Union{Nothing, Outputs}
end

RunData(dir, outputs, log, redirect, checkpointer, provided, backend) =
    RunData(dir, outputs, log, redirect, checkpointer, provided, backend, 0, 0., time(), false, 0,
            0, nothing, IdDict{Any, Vector}(), 0, Any[], false, nothing)

# whether the run `run` gives the value `:sample`, a simulation with samples
gives_sample(run::RunData) = run.outputs.sampled

# the values the run gives a simulation with samples, when `sampled`, and one making trajectories,
# when `trajectories`, which no phase may name
given_values(sampled::Bool, trajectories::Bool) =
    ((sampled ? (:sample,) : ())..., (trajectories ? (:log_weight, :jumps, :jump_list) : ())...)

given_values(run::RunData) = given_values(gives_sample(run), run.outputs.weighted)

"""
    Simulation(sim, state[, time = sim.time]; values...)

the simulation `runTMS` runs and returns: a state, its time and its values, which the
measurements read, `sim.energy` as `c.energy`. The form above gives `sim` another state, and
possibly another time and other values, sharing its run. A value is a number, an array, a
state, anything a phase leaves for the measurements and the phases after it; a checkpoint keeps
it, so that it must not be changed in place.

- `sim.state`, `sim.time`: the state and the time
- `sim.step`: the steps the phase being run has done
- `sim.rng`: the random number generator of the run, for every draw of a phase, see `SimData`
- `sim.name`: the value `name`, a state named by `CopyState` for instance
- `sim.data`: the rows gathered in each `Data(name)`, under `name`, a table, see `Data`

# Examples

    sim = runTMS(sim_data)
    rows = sim.data["magnetization"]
"""
struct Simulation{S <: Union{Nothing, AbstractState}}
    state::S
    time::Number
    values::Dict{Symbol, Any}
    run::RunData
end

function Simulation(sim::Simulation, state::Union{Nothing, AbstractState}, time::Number = sim.time;
                    values...)
    for name in keys(values)
        if name in reserved_values || name in given_values(getfield(sim, :run))
            error("a simulation cannot be given a value named $name, a name the run reserves")
        end
    end
    return Simulation(state, time, merge(getfield(sim, :values), values), getfield(sim, :run))
end

"""
    run_backend(sim)

the type of states the phases of the run of `sim` create by default, see `SimData`
"""
run_backend(sim::Simulation) = getfield(sim, :run).backend

"""
    with_values(sim; values...)

the simulation with the values given, the reserved ones included, which the run gives
"""
with_values(sim::Simulation; values...) =
    Simulation(sim.state, sim.time, merge(getfield(sim, :values), values), getfield(sim, :run))

"""
    data_rows(sim)

the rows gathered in each `Data` of the run of `sim`, by its name, their averages over the
samples when they are averaged; in a sample, those of the samples merged
"""
function data_rows(sim::Simulation)
    run = getfield(sim, :run)
    o = run.outputs
    if o.recording
        o = !isnothing(run.merged) ? run.merged :
            data_outputs(run.dir, o.format, last(read_merged(run.dir)))
    end
    return lock(o.lock) do
        names = unique(String[ n for n in map(data_name, collect(keys(o.opened)))
                               if !isnothing(n) ])
        Dict{String, DataTable}(n => data_rows(o, n) for n in names)
    end
end

function getproperty(sim::Simulation, name::Symbol)
    if name === :state || name === :time
        return getfield(sim, name)
    elseif name === :data
        return data_rows(sim)
    end
    values = getfield(sim, :values)
    if !haskey(values, name)
        error("the simulation has no value $name")
    end
    return values[name]
end

propertynames(sim::Simulation) = (:state, :time, :data, keys(getfield(sim, :values))...)

show(io::IO, sim::Simulation) = print(io, "Simulation($(sim.state), $(sim.time), ...)")

"""
    Editing(sim)
    Simulation(c)

the simulation `sim` open to changes, as `Update` hands it to its function, read as a
`Simulation` and changed by assignments: `c.state = st`, `c.time = t`, `c.x = v` for a value,
`c.x = nothing` removing it, `c.state = nothing` the state, before a phase creating another;
and the simulation it makes once changed
"""
mutable struct Editing
    state::Union{Nothing, AbstractState}
    time::Number
    values::Dict{Symbol, Any}
    run::RunData
end

Editing(sim::Simulation) =
    Editing(sim.state, sim.time, copy(getfield(sim, :values)), getfield(sim, :run))

Simulation(c::Editing) =
    Simulation(getfield(c, :state), getfield(c, :time), getfield(c, :values), getfield(c, :run))

getproperty(c::Editing, name::Symbol) = getproperty(Simulation(c), name)

propertynames(c::Editing) = propertynames(Simulation(c))

function setproperty!(c::Editing, name::Symbol, x)
    if name === :state && !(x isa Union{Nothing, AbstractState})
        error("c.state takes a state, or nothing to remove it, not a $(typeof(x))")
    elseif name === :time && !(x isa Number)
        error("c.time takes a number, not a $(typeof(x))")
    elseif name === :state || name === :time
        setfield!(c, name, x)
    elseif name in reserved_values || name in given_values(getfield(c, :run))
        error("$name cannot be changed, a name the run reserves")
    elseif isnothing(x)
        delete!(getfield(c, :values), name)
    else
        getfield(c, :values)[name] = x
    end
    return x
end

"""
    stopped(sim)

whether the run `runTMS` returned stopped before the end of its phases, at `max_time`, on the
file `stop` or on an interrupt: running it again resumes it

# Examples

    sim = runTMS(sim_data)
    if stopped(sim)
        println("stopped at time ", sim.time, ", to be resumed")
    end
"""
stopped(sim::Simulation) = getfield(sim, :run).stopped

"""
    state_named(sim, name)

the state of the simulation named `name`, `:state` the current one
"""
function state_named(sim::Simulation, name::Symbol)
    st = name === :state ? sim.state : get(getfield(sim, :values), name, nothing)
    if !(st isa AbstractState)
        error("the simulation has no state named $name")
    end
    return st
end

"""
    with_state(sim, name, state; values...)

the simulation with `state` under the name `name`, `:state` the current one, and the `values`
"""
with_state(sim::Simulation, name::Symbol, state; values...) =
    name === :state ? Simulation(sim, state; values...) :
                      Simulation(sim, sim.state; (name => state,)..., values...)

"""
    state_names(on)

the names of the states a phase acts on, its field `on` a name or a vector of them
"""
state_names(on::Symbol) = [on]
state_names(on::Vector{Symbol}) = on

"""
    run_destination(run, d)

where the run writes what goes to `d`: the stream of a redirected run for a name, `d` otherwise,
a name of a file of the simulation refused
"""
function run_destination(run::RunData, name::AbstractString)
    if !isnothing(run.redirect)
        return Stream(run.redirect, run.outputs.format)
    elseif !isempty(run.dir) && is_simulation_file(name, gives_sample(run))
        error("cannot write to $name, a file of the simulation directory: choose another name")
    end
    return name
end

run_destination(run::RunData, d::EachSample) = EachSample(run_destination(run, d.inner))

run_destination(::RunData, d) = d

"""
    row_at!(run, phase, step)

the place where the rows measured from now on are written, see `put_keyed_row!`: the step `step`
of the phase `phase` of the sample being run
"""
row_at!(run::RunData, phase::Int, step::Int) = (run.outputs.at[] = (run.sample, phase, step))

"""
    measured_values(sim)

the values the measurements of `sim` read, `missing` for those its phases provide that are not
there yet
"""
function measured_values(sim::Simulation)
    values = getfield(sim, :values)
    return merge(Dict{Symbol, Any}(k => missing for k in getfield(sim, :run).provided), values)
end

"""
    output(sim, destination => measurements)
    output(sim, [destination1 => measurements1, ...])

measure each list on the state of the simulation, at its time and with its values, and write it
to its destination: a file name, `"stdout"`, `"stderr"`, `""` or `Data(name)`, see `Destination`
and `measure`. The files are taken in the directory of the simulation, whose own files are
refused.

# Examples

    output(sim, "magnetization" => [X, Z(1)])
    output(sim, [ "energy" => :energy, Data("purity") => Purity ])
"""
function output(sim::Simulation, list::AbstractVector)
    if isempty(list)
        return nothing
    elseif isnothing(sim.state)
        error("the simulation has no state to measure: a phase must create one first")
    end
    run = getfield(sim, :run)
    items = [ run_destination(run, first(p)) => last(p) for p in list ]
    values = measured_values(sim)
    # the weight of a trajectory, what its steps carried and what its trace holds, one for any
    # other state of a simulation making trajectories
    if sim.state isa AbstractState{Trajectory}
        values[:log_weight] = get(values, :log_weight, 0.) + log(real(trace(sim.state)))
    elseif run.outputs.weighted
        values[:log_weight] = get(values, :log_weight, 0.)
    end
    if run.outputs.weighted
        values[:jumps] = get(values, :jumps, 0)
        values[:jump_list] = get(values, :jump_list, [])
    end
    # called for the state as it is: inferred for any state, the measurements would take
    # minutes to compile
    invokelatest(output!, run.outputs, sim.state, items; time = sim.time, values...)
    return nothing
end

output(c::Editing, list::AbstractVector) = output(Simulation(c), list)
output(sim::Union{Simulation, Editing}, p::Pair) = output(sim, [p])

"""
    phase_output(sim, list)

`output(sim, list)` for a list of measurements of the phase being run, or of the simulation,
measured again and again: the run makes it once into `Measure`s, one for each destination, kept
until the next phase starts, see `run_phases`, so that each keeps what it prepares for the
system measured, see `Measure`
"""
function phase_output(sim::Simulation, list)
    run = getfield(sim, :run)
    items = get!(run.measures, list) do
        pairs = list isa Pair ? [list] : list
        by_destination(run.outputs, [ run_destination(run, first(p)) => last(p) for p in pairs ])
    end
    return output(sim, items)
end

function log_line(run::RunData, text)
    println(run.log, text)
    flush(run.log)
end

"""
    log_message(sim, text)

write the line `text` to the log of the simulation, or to the stream its output is redirected
to, at once
"""
function log_message(sim::Union{Simulation, Editing}, text)
    log_line(getfield(sim, :run), text)
    return nothing
end

"""
    log_ignored(sim, state, phase, names)

note in the log the options `names` of `phase` that the type of `state` ignores, see
`ignored_options`, those given only: an algorithm other than the default `Tdvp()`, limits other
than `Limits()`, a link dimension of `randomize` other than 0
"""
function log_ignored(sim::Simulation, state::AbstractState, phase, names)
    given(x::Limits) = x ≠ Limits()
    given(x::EvolutionAlgorithm) = x ≠ Tdvp()
    given(x::Integer) = !iszero(x)
    given(x) = !isnothing(x)
    label(n) = n === :randomize ? "the link dimension of randomize" : string(n)
    ignored = [ label(n) for n in names
                if n in ignored_options(state) && given(getfield(phase, n)) ]
    if !isempty(ignored)
        log_message(sim, "$(nameof(typeof(state))) ignores $(join(ignored, ", ", " and "))")
    end
end

"""
    time_text(run, t)
    time_text(sim, t)

the time `t` as the text tables of the run write it, in their format: `0.3` where Julia prints
`0.30000000000000004`
"""
time_text(run::RunData, t) = text_cell(t, Printf.Format(run.outputs.format))
time_text(sim::Union{Simulation, Editing}, t) = time_text(getfield(sim, :run), t)

"""
    report_dropped(run)

one warning in the log for the imaginary parts of Hermitian measurements dropped since the last
one, how many and the largest, which a run gathers rather than writing a line for each
"""
function report_dropped(run::RunData)
    if run.dropped > 0
        x, name, time = run.largest
        n = run.dropped
        at = isnothing(time) ? "" : " at time $(time_text(run, time))"
        log_line(run, "WARNING: $n imaginary part$(n == 1 ? "" : "s") of Hermitian measurements " *
                      "beyond what an approximation leaves dropped, the largest $x, of $name$at")
        run.warnings += 1
        run.dropped = 0
        run.largest = nothing
    end
    return nothing
end

"""
    SimLogger(run, parent)

the logger of a run: the warnings of this package go to its log, those of the imaginary parts
dropped by measurements gathered for `report_dropped`, everything else to `parent`, the logger
in place
"""
struct SimLogger <: Logging.AbstractLogger
    run::RunData
    parent::Logging.AbstractLogger
end

"""
    simulation_warning(level, _module)

whether a log message is a warning of this package, which `SimLogger` sends to the log of the
run
"""
simulation_warning(level, _module) =
    level >= Logging.Warn && !isnothing(_module) && Base.moduleroot(_module) === @__MODULE__

Logging.min_enabled_level(l::SimLogger) = min(Logging.Warn, Logging.min_enabled_level(l.parent))

Logging.shouldlog(l::SimLogger, level, _module, group, id) =
    simulation_warning(level, _module) || Logging.shouldlog(l.parent, level, _module, group, id)

Logging.catch_exceptions(l::SimLogger) = Logging.catch_exceptions(l.parent)

function Logging.handle_message(l::SimLogger, level, message, _module, group, id, file, line;
                                kwargs...)
    if simulation_warning(level, _module) && group === :imaginary_part
        run = l.run
        x = kwargs[:imaginary]
        if isnothing(run.largest) || abs(x) > abs(first(run.largest))
            run.largest = (x, kwargs[:measurement], kwargs[:time])
        end
        run.dropped += 1
    elseif simulation_warning(level, _module)
        l.run.warnings += 1
        log_line(l.run, "WARNING: $message")
    elseif level >= Logging.min_enabled_level(l.parent)
        # the level of this logger is the lower of the two: the parent's is checked again
        Logging.handle_message(l.parent, level, message, _module, group, id, file, line; kwargs...)
    end
    return nothing
end
