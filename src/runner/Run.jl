# runTMS, which runs a simulation described by a SimData phase after phase, or sample after sample
# from a branch point, in a directory of its own, with its log, its checkpoints and the resumption
# of an interrupted run, and cleanTMS, which removes that directory.

export runTMS, cleanTMS, SimData

"""
    measurement_lists(x)

the lists of measurements `output` takes, `destination => measurements` or a vector of them, as
a vector of pairs
"""
measurement_lists(p::Pair) = [p]
measurement_lists(v::AbstractVector) = v
measurement_lists(x) =
    error("$(repr(x)) is not a list of measurements, destination => measurements")

"""
    run_lists(phases, final_measurements)

the lists of measurements of a simulation, each the pairs `destination => measurements` that one
call writes, see `measurement_lists`: its final ones, those of the end of each phase and those
of their steps
"""
function run_lists(phases::Vector, final_measurements)
    lists = Any[final_measurements]
    append!(lists, [ q.final_measurements for q in phases ])
    append!(lists, map(step_measurements, phases))
    return [ measurement_lists(l) for l in lists ]
end

"""
    refused_reads(phase)

the values the lists of measurements of `phase` cannot read, each with the reason given when one
does: none by default
"""
refused_reads(::AbstractPhase) = ()

"""
    check_values(phases, final_measurements, sampled)

refuse a value measured by one of the lists of the phases, or by the final measurements, that no
phase provides, see `provided_values`, or that its phase refuses, see `refused_reads`, and
return the names of those they provide; `:sample` is read only when the simulation is `sampled`
"""
function check_values(phases::Vector, final_measurements, sampled::Bool)
    for q in phases, (value, why) in refused_reads(q)
        lists = Any[q.final_measurements, step_measurements(q)]
        for l in lists, p in measurement_lists(l), (name, m) in as_measure(last(p)).entries
            if value in measurement_reads(m)
                error("the measurement $name of a $(nameof(typeof(q))) reads $value: $why")
            end
        end
    end
    provided = unique(Symbol[ v for p in phases for v in provided_values(p) ])
    for l in run_lists(phases, final_measurements), p in l
        for (name, m) in as_measure(last(p)).entries
            for r in measurement_reads(m)
                if r === :sample && !sampled
                    error("the measurement $name reads sample, which only a simulation with " *
                          "samples has, see Realizations")
                elseif r ∉ provided && r ∉ (:step, :rng, :sample)
                    of_1x = r === :sweep ? ": the :sweep of 1.x is :step" : ""
                    error("the measurement $name reads $r, a value that no phase of the " *
                          "simulation provides$of_1x")
                end
            end
        end
    end
    return provided
end

"""
    check_names(phases)

refuse, when the `SimData` is written rather than when the run reaches the phase, a phase naming a
state by a name the run reserves, see `reserved_values`, `:state` aside, and one naming a state
twice in its `on`, which would act on it twice
"""
function check_names(phases::Vector)
    for p in phases
        names = acted_on(p)
        for n in [ names ; collect(provided_values(p)) ]
            if n ≠ :state && n in reserved_values
                error("$(type_name(typeof(p))) cannot name a state $n, a name the run reserves")
            end
        end
        for n in names
            if count(==(n), names) > 1
                error("$(type_name(typeof(p))) names the state $n twice in its on, which would " *
                      "act on it twice")
            end
        end
    end
end

"""
    check_destinations(phases, final_measurements, sampled; live = true)

refuse, before a run starts rather than when it writes there, in a `sampled` simulation a file or
`Data` both averaged and written for each sample, and in a run with a directory (`live`) a file
of the simulation directory as the destination of a list of measurements of the phases or of the
final measurements, or as a file a phase writes, see `phase_files` and `simulation_files`, and a
text file given two lists, as one holds one list, the lists one call sends it making one
"""
function check_destinations(phases::Vector, final_measurements, sampled::Bool; live::Bool = true)
    for p in phases, f in phase_files(p)
        if live && is_simulation_file(f)
            error("cannot save to $f, a file of the simulation directory: choose another name")
        end
    end
    held = Dict{String, Vector{String}}()
    averaged = Dict{Any, Bool}()
    for l in run_lists(phases, final_measurements), (a, m) in by_destination(Outputs(; sampled), l)
        d = a isa Union{Averaged, EachSample} ? a.inner : a
        # without a directory, a name is the stream the output goes to, and only a Data is kept
        both = (live ? d isa Union{TextFile, JsonFile, Data} : d isa Data) &&
               get!(averaged, d, a isa Averaged) ≠ (a isa Averaged)
        if both
            error("$d is written both averaged and for each sample: give EachSample(…) another " *
                  "destination")
        end
        if !live || !(d isa Union{TextFile, JsonFile})
            continue
        elseif is_simulation_file(d.path)
            error("cannot write to $(d.path), a file of the simulation directory: choose another " *
                  "name")
        elseif d isa TextFile && get!(held, d.path, first.(m.entries)) ≠ first.(m.entries)
            error("$(d.path) holds the list $(join(held[d.path], ", ")), and another needs " *
                  "another file")
        end
    end
end

"""
    check_format(format)

refuse a Printf format that does not write one number
"""
function check_format(format::AbstractString)
    Outputs(; format)
    return format
end

"""
    checked_phases(phases, realized, samples)
    checked_phases(sd)

the phases the checks of a simulation see: those of its dry call in place of `Realizations`
"""
checked_phases(phases::Vector, realized::Vector, samples) =
    isnothing(samples) ? phases : [phases[1:end - 1]; realized]

"""
    SimData(; name = "simulation", phases, options...)

the description of a simulation, which `runTMS` runs.

- `name`: the name of the simulation, and of the directory its results are written to
- `phases`: the phases, see `TensorMixedStates.AbstractPhase`, as a vector which may contain
  vectors to any depth, flattened
- `description`: the text of the file `description` of the simulation (default `""`)
- `time_start`: the time the simulation starts at (default 0)
- `final_measurements`: the measurements made at the end, see `output` (default `[]`)
- `format`: the Printf format of the numbers written to text files (default `"%.13g"`, thirteen
  significant digits, the precision below which a difference is taken for rounding)
- `seed`: the seed of the random number generator of the run, `sim.rng`, an integer of zero or
  above, drawn at random by default and written in the file `stamp`, so that a run can be made
  again
- `backend`: the type of states `CreateState` makes when its own `backend` is not given, `State`
  (MPS) by default or `ExactState`, so that the phases of a program need no change to run on
  another backend, `LoadState` keeping the type of its file and a phase the backend cannot run
  raising an error
- `samples`: the number of samples of a simulation whose last phase is `Realizations`, which
  draws the disorder of each, the final measurements being made at the end of each sample, the
  files and `Data` holding averages over the samples, see `EachSample` (default `nothing`, no
  samples); the environment variable `TMS_SAMPLES` overrides it, so that a finished run is
  prolonged by running it again with more, and the program run several times at once runs the
  samples on as many cores, see `runTMS`
- `tolerance`: the errors the samples stop at, for an averaged destination `d`, a file or a
  `Data` (a `Data` alone in a run without a directory), `d => atol` or
  `d => (atol = …, rtol = …)`, or a list of them: every error of `d`, `err(c)`, within
  `atol + rtol * |c|`, at every row, the samples stopping at the first that reaches it, the same
  whichever process runs them, those an earlier launch merged kept (default `[]`, all the
  `samples`)
- `min_samples`: the samples each row of `d` holds before its tolerance is looked at, so that it
  does not stop on an error underestimated by chance (default 10)

A measurement of a value that no phase provides is refused here, see
`TensorMixedStates.provided_values`, and so is a state named by a name the run reserves, or twice
in the `on` of a phase.

# Examples

    SimData(
        name = "ising",
        phases = [
            CreateState{Pure}(10, Qubit(), "Up"),
            GroundState(hamiltonian = -sum(Z(i) * Z(i + 1) for i in 1:9) - sum(X(i) for i in 1:10),
                        algo = Dmrg(nsweeps = 10), limits = Limits(maxdim = [10, 20, 50]),
                        final_measurements = "data" => [:energy, X, Z]),
        ],
    )
"""
struct SimData
    name::String
    description::String
    time_start::Number
    final_measurements
    format::String
    seed::Union{Nothing, Integer}
    backend::Type
    samples::Union{Nothing, Int}
    tolerance::Vector{Pair{Any, @NamedTuple{atol::Float64, rtol::Float64}}}
    min_samples::Int
    phases::Vector{Any}
    provided::Vector{Symbol}
    realized::Vector{Any}
    function SimData(name, description, time_start, final_measurements, format, seed, backend,
                     samples, tolerance, min_samples, phases)
        # which Julia 1.10 refuses
        if !isnothing(seed) && seed < 0
            error("the seed of a simulation is an integer of zero or above, got $seed")
        elseif !isnothing(samples) && samples < 1
            error("the number of samples of a simulation is 1 or more, got $samples")
        elseif !isnothing(samples) && !isnothing(seed) && seed > typemax(UInt64)
            error("the seed of a simulation with samples is at most typemax(UInt64), got $seed")
        end
        samples = isnothing(samples) ? nothing : Int(samples)
        ps = map(check_phase, flatten_phases(phases))
        realized = branched_phases(ps, samples)
        checked = checked_phases(ps, realized, samples)
        provided = check_values(checked, final_measurements, !isnothing(samples))
        check_names(checked)
        lists = run_lists(checked, final_measurements)
        if isnothing(samples) && any(first(p) isa EachSample for l in lists for p in l)
            error("EachSample writes a row for each sample, and the simulation has no samples")
        end
        tolerance = tolerances(tolerance)
        if !isempty(tolerance) && isnothing(samples)
            error("a tolerance stops the samples, and the simulation has no samples")
        elseif min_samples < 1
            error("min_samples is 1 or more, got $min_samples")
        end
        for (d, _) in tolerance
            if !(plain_destination(Outputs(), d) isa Union{TextFile, JsonFile, Data})
                error("the tolerance names $(repr(d)), which is not averaged: name a file or a " *
                      "Data, as \"energy\" or Data(\"x\")")
            elseif !any(first(p) == d for l in lists for p in l)
                error("the tolerance names $(repr(d)), which no list of measurements writes " *
                      "averaged: name a destination of the lists, as \"energy\" or Data(\"x\")")
            end
        end
        return new(name, description, time_start, final_measurements, check_format(format), seed,
                   check_backend(backend), samples, tolerance, Int(min_samples), ps, provided,
                   realized)
    end
end

# the fields of 2.0, before backend, and of 2.1, before samples
SimData(name, description, time_start, final_measurements, format, seed, phases) =
    SimData(name, description, time_start, final_measurements, format, seed, State, phases)

SimData(name, description, time_start, final_measurements, format, seed, backend, phases) =
    SimData(name, description, time_start, final_measurements, format, seed, backend, nothing, [],
            10, phases)

SimData(; name::AbstractString = "simulation", description::AbstractString = "",
        time_start::Number = 0., final_measurements = [], format::AbstractString = "%.13g",
        seed::Union{Nothing, Integer} = nothing, backend = State,
        samples::Union{Nothing, Integer} = nothing, tolerance = [], min_samples::Integer = 10,
        phases) =
    SimData(name, description, time_start, final_measurements, format, seed, backend, samples,
            tolerance, min_samples, phases)

"""
    tolerances(t)

the tolerance of a `SimData` as a vector of pairs `destination => (atol, rtol)`, from a pair
`d => atol` or `d => (atol = …, rtol = …)` or a vector of them
"""
tolerances(p::Pair) = tolerances([p])

function tolerances(v::AbstractVector)
    out = Pair{Any, @NamedTuple{atol::Float64, rtol::Float64}}[]
    for p in v
        tol = p isa Pair ? last(p) : nothing
        t = tol isa Real ? (atol = Float64(tol), rtol = 0.) :
            tol isa NamedTuple && issubset(keys(tol), (:atol, :rtol)) &&
            all(x -> x isa Real, values(tol)) ?
            (atol = Float64(get(tol, :atol, 0.)), rtol = Float64(get(tol, :rtol, 0.))) : nothing
        if isnothing(t) || !(t.atol ≥ 0 && t.rtol ≥ 0)
            error("a tolerance is destination => atol or destination => (atol = …, rtol = …), " *
                  "of zero or above, got $(repr(p))")
        end
        push!(out, first(p) => t)
    end
    return out
end

tolerances(x) =
    error("a tolerance is destination => atol or destination => (atol = …, rtol = …), or a " *
          "list of them, got $(repr(x))")

checked_phases(sd::SimData) = checked_phases(sd.phases, sd.realized, sd.samples)

# what the merge of the samples of `sd` needs, see `merge_samples!`
merging(sd::SimData) = (; sd.format, sd.tolerance, sd.min_samples)

show(io::IO, s::SimData) =
    print(io, """
        SimData(
            name = $(repr(s.name)),
            description = $(repr(s.description)),
            time_start = $(s.time_start),
            final_measurements = $(s.final_measurements),
            format = $(repr(s.format)),
            seed = $(repr(s.seed)),
            backend = $(s.backend),
            samples = $(repr(s.samples)),
            tolerance = $(s.tolerance),
            min_samples = $(s.min_samples),
            phases =
        $(s.phases))""")

"""
    includes

the files included since the package was loaded, each with the file that included it, `nothing`
for an include from the REPL, which `record_include` notes, see `__init__`. A file included
again keeps what it included before, which it may skip this time, as a file it includes once.
"""
const includes = (edges = Set{Pair{Union{Nothing, String}, String}}(), lock = ReentrantLock())

# called by include before it runs the file, the source path being still that of the file
# including it
function record_include(::Module, path::String)
    lock(() -> push!(includes.edges, Base.source_path(nothing) => path), includes.lock)
    return nothing
end

"""
    program_files(src_path)

the files the program `src_path` included, directly or not, after the package was loaded and
still there, each by its path relative to the directory of the program, a `..` taken as
`_parent`, as the directory `prog` of a simulation holds it
"""
function program_files(src_path::String)
    files = Dict{String, String}()
    todo = [src_path]
    edges = lock(() -> copy(includes.edges), includes.lock)
    while !isempty(todo)
        f = pop!(todo)
        for (from, path) in edges
            # one removed since cannot be compared
            if from == f && path ≠ src_path && path ∉ values(files) && isfile(path)
                rel = joinpath(replace(splitpath(relpath(path, dirname(src_path))),
                                       ".." => "_parent")...)
                if haskey(files, rel)
                    error("$(files[rel]) and $path would both be copied to " *
                          "$(joinpath("prog", rel)): rename one of them")
                end
                files[rel] = path
                push!(todo, path)
            end
        end
    end
    return files
end

"""
    stored_files(dir)

the files of the directory `prog` of the simulation directory `dir`, by their paths relative to
it
"""
function stored_files(dir::String)
    prog = joinpath(dir, "prog")
    if !isdir(prog)
        return String[]
    end
    return [ relpath(joinpath(root, f), prog) for (root, _, fs) in walkdir(prog) for f in fs ]
end

"""
    same_program(dir, src_path)

whether the checkpoint of the directory `dir` may be resumed by the program `src_path` given the
arguments `ARGS`: the same, byte for byte, as `prog.jl`, including the same files as those of
the directory `prog`, given the arguments `prog_args.json` holds, none when it is absent. A
program with no file, run from the REPL, is always accepted.
"""
function same_program(dir::String, src_path)
    if isnothing(src_path) || src_path == ""
        return true
    end
    prog, args_file = joinpath(dir, "prog.jl"), joinpath(dir, "prog_args.json")
    args = isfile(args_file) ? JSON.parsefile(args_file) : []
    if !(isfile(prog) && read(src_path) == read(prog) && args == ARGS)
        return false
    end
    files = program_files(src_path)
    return sort(collect(keys(files))) == sort(stored_files(dir)) &&
           all(read(path) == read(joinpath(dir, "prog", rel)) for (rel, path) in files)
end

"""
    copy_program_files(dir, src_path)

copy the files the program `src_path` included into the directory `prog` of the simulation
directory `dir`, see `program_files`, in place of those it held
"""
function copy_program_files(dir::String, src_path::String)
    files = program_files(src_path)
    prog = joinpath(dir, "prog")
    # read before the directory is emptied, which may hold them when the copy of the program is
    # run
    sources = Dict(rel => read(path) for (rel, path) in files)
    rm(prog; recursive = true, force = true)
    for (rel, data) in sources
        mkpath(dirname(joinpath(prog, rel)))
        write(joinpath(prog, rel), data)
    end
    return nothing
end

"""
    active_runs

the runs going on in the process, the directories they hold and the threading before the first
of them, under a lock: only the outermost run sets the threading and the behaviour of Ctrl-C,
which the last to end puts back, and a directory is held by one run of the process at a time,
see `claim_directory`
"""
const active_runs = (count = Ref(0), threading = Ref{Union{Nothing, ThreadingState}}(nothing),
                     dirs = Set{String}(), lock = ReentrantLock())

"""
    claim_directory(dir)

take the directory `dir` for a run, refused when another run of the process holds it, a nested
one or one on another thread, whose file `running`, naming this process, would be taken for one
a dead process left, see `live_marker`
"""
function claim_directory(dir::String)
    lock(active_runs.lock) do
        if dir in active_runs.dirs
            error("$dir is the directory of another run of this process, going on: give each run " *
                  "a name of its own")
        end
        push!(active_runs.dirs, dir)
    end
    return nothing
end

release_directory(dir::String) = lock(() -> delete!(active_runs.dirs, dir), active_runs.lock)

"""
    live_marker(text)

whether the file `running` holding `text`, the machine, the process and the time a run started,
marks a run that may still be going on: not when it names a process of this machine that no
longer exists, or this very process, whose number the dead one had, the runs of this process
holding their directories apart, see `claim_directory`. A file naming another machine is taken
as live, there being no way to know.
"""
function live_marker(text::AbstractString)
    parts = split(text)
    if length(parts) < 2 || parts[1] ≠ gethostname() || isnothing(tryparse(Int32, parts[2]))
        return true
    end
    pid = parse(Int32, parts[2])
    # libuv answers signal 0 on every system: 0 for a process that exists, a refusal of
    # permission for one of another user, and ESRCH for none
    return pid ≠ getpid() && ccall(:uv_kill, Cint, (Cint, Cint), pid, 0) ≠ Base.UV_ESRCH
end

"""
    check_removable(name)

refuse to remove the directory `name` when it contains the current directory: `rm` would empty
it, the program included, before failing
"""
function check_removable(name::AbstractString)
    if ispath(name) && !islink(name) &&
       startswith(joinpath(realpath(pwd()), ""), joinpath(realpath(name), ""))
        error("cannot remove \"$name\", which contains the current directory")
    end
end

"""
    check_unlocked(dir)

refuse a directory whose file `running` marks a run that may still be going on, see
`live_marker`, and return the parts of a stale one, `nothing` without one
"""
function check_unlocked(dir::String)
    file = joinpath(dir, "running")
    if !isfile(file)
        return nothing
    end
    marker = read(file, String)
    if live_marker(marker)
        parts = split(marker)
        since = length(parts) ≥ 3 ? ", started $(parts[3])" : ""
        error("$dir holds the file running of another run$since, which may still be going on: " *
              "if none is, remove $(joinpath(dir, "running"))")
    end
    return split(marker)
end

"""
    cleanTMS(sim_data)
    cleanTMS(name)

remove the directory of a simulation, refused while a run may be using it

# Examples

    cleanTMS(sim_data)
"""
cleanTMS(sd::SimData) = cleanTMS(sd.name)

function cleanTMS(name::AbstractString)
    if isempty(name)
        error("a simulation run in a directory needs a name, that of its directory")
    end
    check_removable(name)
    if isdir(name)
        if lock(() -> abspath(name) in active_runs.dirs, active_runs.lock)
            error("$(abspath(name)) is the directory of a run of this process, going on")
        end
        check_unlocked(abspath(name))
        for p in part_numbers(abspath(name))
            check_unlocked(part_dir(abspath(name), p))
        end
    end
    rm(name; recursive = true, force = true)
    return nothing
end

"""
    launch_setting(name, value, parse)

the setting of `runTMS` given by the environment variable `name`, read by `parse`, which gives
`nothing` for what it cannot read, `value` when the variable is not set
"""
function launch_setting(name::String, value, parse)
    s = get(ENV, name, "")
    if isempty(s)
        return value
    end
    x = try
        parse(s)
    catch
        nothing
    end
    if isnothing(x)
        error("the environment variable $name is $(repr(s)), which runTMS cannot read")
    end
    return x
end

"""
    threading_setting(s)

the threading named by the environment variable `TMS_THREADING`
"""
function threading_setting(s::AbstractString)
    t = lstrip(s, ':')
    if t == "nothing"
        return Some(nothing)
    elseif t in ("auto", "dense", "blocks")
        return Some(Symbol(t))
    end
    return nothing
end

"""
    check_threading(threading)

refuse a `threading` other than `:auto`, `:dense`, `:blocks` and `nothing`
"""
function check_threading(threading)
    if !(threading in (nothing, :auto, :dense, :blocks))
        error("unknown threading $(repr(threading)), expected :auto, :dense, :blocks or nothing")
    end
    return threading
end

"""
    threading_stamp(mode)

the lines of the `stamp` file saying how the run is threaded when it starts, see
`threading_settings`
"""
function threading_stamp(mode)
    s = threading_settings()
    return """
        Threading $(repr(mode))
        BLAS $(s.blas_library), $(s.blas) threads
        Julia threads $(s.julia)
        GC threads $(s.gc)
        Strided threads $(s.strided)
        Block sparse multithreading $(s.blocksparse ? "on" : "off")
        CPU threads $(Sys.CPU_THREADS)
        """
end

"""
    computing_packages

the packages doing the computations of a run, whose versions the `stamp` file gives
"""
const computing_packages = ["ITensors", "ITensorMPS", "NDTensors", "KrylovKit"]

"""
    packages_stamp()

the lines of the `stamp` file giving the versions of the computing packages, and of the other
packages the program loads into `Main`, but for the standard library
"""
function packages_stamp()
    loaded = Dict(id.name => m for (id, m) in Base.loaded_modules)
    version(m) = "$(nameof(m)) $(pkgversion(m))"
    stamp = "Packages " * join([ version(loaded[name]) for name in computing_packages
                                 if haskey(loaded, name) ], ", ") * "\n"
    others = [ getfield(Main, n) for n in names(Main; imported = true) if isdefined(Main, n) ]
    others = sort([ version(m) for m in others
                    if m isa Module && parentmodule(m) === m && !isnothing(pkgversion(m)) &&
                       !startswith(something(pkgdir(m), ""), Sys.STDLIB) &&
                       m !== TensorMixedStates && string(nameof(m)) ∉ computing_packages ])
    return isempty(others) ? stamp : stamp * "Program packages " * join(others, ", ") * "\n"
end

function enter_run()
    lock(active_runs.lock) do
        active_runs.count[] += 1
        if active_runs.count[] == 1
            active_runs.threading[] = ThreadingState()
            # an exception rather than an exit, so that the run writes its checkpoint
            Base.exit_on_sigint(false)
            return true
        end
        return false
    end
end

function leave_run()
    lock(active_runs.lock) do
        active_runs.count[] -= 1
        if active_runs.count[] == 0
            set_threading(active_runs.threading[])
            # impossible to read back: the default of Julia, on in a script, off in the REPL
            Base.exit_on_sigint(!isinteractive())
        end
    end
end

"""
    run_simulation(sim, sd, samples, sample, start, done, ended, threading)

run the phases and the final measurements of `sd` from where `Commit` says, or for a simulation
with samples run without a directory, each sample in turn, see `run_samples`
"""
function run_simulation(sim::Simulation, sd::SimData, samples::Union{Nothing, Int}, sample::Int,
                        start::Int, done::Int, ended::Bool, threading)
    if isnothing(samples)
        sim = run_phases(sim, sd.phases, start, done, ended, threading)
        row_at!(getfield(sim, :run), typemax(Int), 0)
        phase_output(sim, sd.final_measurements)
        return sim
    end
    return run_samples(sim, sd, samples, threading)
end

"""
    sample_threading

the threading of a sample, one thread, so that it is the same to the last bit whichever process
runs it: the threads of BLAS change the rounding of its products
"""
const sample_threading = ThreadingState(1, 1, false)

# the commit of the branch point of `sim`, the phases before `branch` done
branch_commit(sim::Simulation, branch::Int) =
    Commit(0, branch, 0, false, getfield(sim, :run).phase_time, sim.time, sim.state,
           committed_values(getfield(sim, :values)), Dict{Destination, Any}())

"""
    run_prefix(sim, sd, commit, threading)

run the phases of `sd` before its branch point, from `commit` or anew, then write the merge of no
sample, which holds their rows, and `branch.h5` last, which the parts wait for
"""
function run_prefix(sim::Simulation, sd::SimData, commit::Union{Nothing, Commit}, threading)
    run = getfield(sim, :run)
    c = run.checkpointer
    prefix = sd.phases[1:end - 1]
    start, done, ended = isnothing(commit) ? (1, 0, false) :
                         (commit.phase, commit.step, commit.last)
    sim = run_phases(sim, prefix, start, done, ended, threading; final = false)
    write_merged(c.dir, 0, resume_marks(run.outputs))
    # its files written now, before a merge writes them
    close!(run.outputs)
    empty!(run.outputs.opened)
    write_branch(c.dir, branch_commit(sim, length(prefix) + 1), c.seed)
    log_line(run, "Branch point written, the samples start from it")
    return sim
end

"""
    run_sample(sim, sd, b, sample, resumed, threading)

run the sample `sample` of `sd` from the branch point `b`, or from its commit `resumed`, which
`sim` holds: its phases, drawn by the function of `Realizations`, and the final measurements
"""
function run_sample(sim::Simulation, sd::SimData, b::Commit, sample::Int,
                    resumed::Union{Nothing, Commit}, threading)
    run = getfield(sim, :run)
    c = run.checkpointer
    prefix = sd.phases[1:end - 1]
    run.sample = sample
    if isnothing(resumed)
        values = merge(b.values, Dict{Symbol, Any}(:sample => sample, :step => 0,
                                                    :rng => sample_rng(c.seed, sample)))
        sim = Simulation(b.state, b.time, values, run)
        start, done, ended = length(prefix) + 1, 0, false
        commit!(sim, start, 0)
    else
        start, done, ended = resumed.phase, resumed.step, resumed.last
    end
    began = time()
    # no phase runs while Realizations draws: an error or a stop there names the sample only
    run.phase = 0
    run.phases = Any[]
    realized = realized_phases(last(sd.phases), Realization(sample, disorder_rng(c.seed, sample)))
    check_realization(realized, sd.realized, sample)
    phases = [prefix; realized]
    check_values(phases, sd.final_measurements, true)
    check_names(realized)
    check_destinations(phases, sd.final_measurements, true; live = !isempty(c.dir))
    sim = run_phases(sim, phases, start, done, ended, threading; final = false)
    # what follows is told of the sample, outside its phases
    run.phase = 0
    row_at!(run, typemax(Int), 0)
    phase_output(sim, sd.final_measurements)
    log_line(run, "Sample $sample done, $(round(time() - began; digits = 3)) seconds" *
                  (sim.state isa State ? ", maxdim $(maxlinkdim(sim.state))" : ""))
    return sim
end

"""
    run_samples(sim, sd, samples, threading)

in a run without a directory, the phases before the branch point, then each sample, its rows
recorded, then added to the outputs, as a part records them and the merge adds them, see
`run_part`
"""
function run_samples(sim::Simulation, sd::SimData, samples::Int, threading)
    run = getfield(sim, :run)
    if any(!(first(t) isa Data) for t in sd.tolerance)
        error("without a directory, the files are not written and a tolerance names a " *
              "Data(name) alone")
    end
    prefix = sd.phases[1:end - 1]
    sim = run_phases(sim, prefix, 1, 0, false, threading; final = false)
    b = branch_commit(sim, length(prefix) + 1)
    merged = run.outputs
    run.merged = merged
    run.outputs = Outputs(; merged.dir, merged.format, sampled = true, recording = true)
    each = isnothing(threading) ? nothing : sample_threading
    try
        for k in 1:samples
            run.quiet = k > 1
            sim = run_sample(sim, sd, b, k, nothing, each)
            replay!(merged, take_records!(run.outputs))
            if !isempty(sd.tolerance) &&
               within(merged.opened, merged.dir, sd.format, sd.tolerance, sd.min_samples)
                log_line(run, "Tolerance reached at sample $k")
                break
            elseif !isempty(sd.tolerance) && k == samples
                log_line(run, "Tolerance not reached at sample $k, the last asked for")
            end
        end
    finally
        close!(run.outputs)
        run.outputs = merged
        run.merged = nothing
        run.quiet = false
    end
    return sim
end

"""
    run_part(sim, sd, b, commit, dir, part, threading)

run the part `part` of the simulation `sd` of directory `dir`: samples taken from its queue, see
`take_sample!`, from the branch point `b`, or from `commit` for the one it holds, their rows
recorded, then written and merged, see `write_pending` and `merge_all!`, until none is left
"""
function run_part(sim::Simulation, sd::SimData, b::Commit, commit::Union{Nothing, Commit},
                  dir::String, part::Int, threading)
    run = getfield(sim, :run)
    c = run.checkpointer
    each = isnothing(threading) ? nothing : sample_threading
    if isnothing(commit)
        # what the part returns when it finds no sample left
        sim = Simulation(b.state, b.time, copy(b.values), run)
    end
    ran = false
    while true
        run.sample = 0
        run.phase = 0
        # a stop is looked for between two samples, where nothing is left to write, once one is
        # left; after a sample, so that a launch whose time is up still runs a step
        k = take_sample!(dir, part; stop = ran && stop_requested(c))
        if k === :stop
            throw(StopRun())
        elseif isnothing(k)
            break
        end
        resumed = !isnothing(commit) && commit.sample == k ? commit : nothing
        if isnothing(resumed)
            clear_records!(run.outputs)
        end
        commit = nothing
        run.quiet = ran
        ran = true
        try
            sim = run_sample(sim, sd, b, k, resumed, each)
        catch e
            if !(e isa StopRun || e isa InterruptException)
                sample_failed!(dir, k, "$(sprint(showerror, e)), see " *
                                       joinpath(part_dir(dir, part), "error"))
            end
            rethrow()
        end
        run.sample = 0
        write_pending(dir, k, take_records!(run.outputs))
        sample_ended!(dir, part)
        merged = merge_all!(dir, merging(sd))
        if !isnothing(merged)
            log_line(run, (merged.first == merged.last ? "Sample $(merged.last) merged" :
                           "Samples $(merged.first) to $(merged.last) merged") *
                          (merged.reached ? ", the tolerance reached" : ""))
        end
    end
    run.quiet = false
    return sim
end

"""
    stop_place(c, phases)

where a run stops, from its last commit
"""
function stop_place(c::Checkpointer, phases::Vector)
    k = c.last
    # no phases while the function of Realizations draws those of a sample
    name = k.phase ≤ length(phases) ? "phase \"$(phases[k.phase].name)\"" : "phase $(k.phase)"
    place = !isempty(phases) && k.phase > length(phases) ? "after the last phase" :
            k.step == 0 ? "before $name" : "after step $(k.step) of $name"
    return k.sample == 0 ? place : "$place of sample $(k.sample)"
end

"""
    stop_line(c, phases)

the line logging where a run stops
"""
function stop_line(c::Checkpointer, phases::Vector)
    then = isempty(c.dir) ? "with no directory to resume from" : "the simulation can be resumed"
    return "***** Stopping $(stop_place(c, phases)), $then *****"
end

"""
    summary_line(sd, run[, samples])

the line the console takes once a run returns: whether it ended or stopped, where its results
are, how many warnings its log received, and what `samples` says of the samples, see
`samples_line`
"""
function summary_line(sd::SimData, run::RunData, samples::AbstractString = "")
    at = stop_place(run.checkpointer, run.phases)
    state = !run.stopped ? "done" :
            isempty(run.dir) ? "stopped $at, with no directory to resume from" :
            "stopped $at: run it again to resume it"
    where = isempty(run.dir) ? "" : ", results in $(run.dir)"
    n = run.warnings
    warned = n == 0 ? "" : isempty(run.dir) ? ", $n warning$(n == 1 ? "" : "s") in its output" :
             ", $n warning$(n == 1 ? "" : "s") in its log"
    return "simulation \"$(sd.name)\" $state$where$warned" * (isempty(samples) ? "" : "; $samples")
end

"""
    error_line(run, sd, e)

the line the log and the file `error` take when the run fails with `e`: the phase it was in and
the time of its last commit
"""
function error_line(run, sd::SimData, e)
    k = isnothing(run) ? nothing : run.checkpointer.last
    where = isnothing(run) || run.phase == run.sample == 0 ? "" :
            run.phase == 0 ? " in sample $(run.sample), outside its phases" :
            " in phase $(run.phase) \"$(run.phases[run.phase].name)\"" *
            (run.sample == 0 ? "" : " of sample $(run.sample)") *
            ", last committed at time " * time_text(run, k.time)
    return "ERROR$where: " * sprint(showerror, e)
end

"""
    runTMS(sim_data; max_time = Inf, checkpoint_interval = 0, threading = :dense, output = nothing)

run the simulation, see `SimData`, in a directory named after it, and return the `Simulation`
it ends with. A checkpoint found in the directory is resumed from, if it was written by the same
program, `prog.jl` in the directory, including the same files, copied into its directory `prog`,
given the same arguments: copying an edited program there, or setting `TMS_ACCEPT_PROGRAM=1`,
resumes it with that one, the variable leaving a line in the log. A file `restart` in the
directory starts the simulation over, and `cleanTMS` removes the directory. A destination that is
a file of the directory, as `log` or `stamp`, or a text file given two lists, is refused before
the run starts.

- `max_time`: seconds after which the run stops, writing a checkpoint
- `checkpoint_interval`: seconds between two checkpoints, 0 for those of a stop only
- `threading`: how the contractions are threaded, see `set_threading`: `:dense`, `:blocks`,
  `:auto`, chosen before each phase from the system of the state, or `nothing`, leaving the
  threading as it is
- `output`: a stream to write everything to, `stdout` or `devnull` for instance, instead of a
  directory; nothing can then be resumed

The first three are overridden by the environment variables `TMS_MAX_TIME`,
`TMS_CHECKPOINT_INTERVAL` and `TMS_THREADING`, `TMS_MAX_TIME=3600 julia prog.jl` for instance.
The run also stops, writing a checkpoint, when the file `stop` appears in its directory, and on
an interrupt. Runs may be nested or run on several threads: the outermost one sets the threading
and the behaviour of Ctrl-C, which the last to end puts back, and an interrupt stops a nested run
and the runs it is in, each resuming where it was. A nested run, as in an `Update`, is done again
when the run it is in resumes, the same way with `output = devnull` and a `seed` only: in a
directory of its own, it resumes its checkpoint or starts over, a seed drawn anew.

The samples of a simulation, see `SimData`, are run by parts, each in a directory `parts/1`,
`parts/2`, … of its own, and merged in their order into its files and `Data`: running the
program several times at once runs them on as many cores, a process waiting while another runs
the phases before the branch point, and returning when no sample is left, with the `Data`
merged so far. Each sample runs on one thread, so that it is the same whichever process runs
it, but with `threading = nothing` and in a nested run, which leave the threading as it is. The
file `stop` of the directory stops all the parts.
"""
function runTMS(sd::SimData; max_time::Real = Inf, checkpoint_interval::Real = 0,
                threading::Union{Nothing, Symbol} = :dense, output::Union{Nothing, IO} = nothing)
    max_time = launch_setting("TMS_MAX_TIME", max_time, s -> parse(Float64, s))
    checkpoint_interval = launch_setting("TMS_CHECKPOINT_INTERVAL", checkpoint_interval,
                                         s -> parse(Float64, s))
    threading = check_threading(something(launch_setting("TMS_THREADING", Some(threading),
                                                         threading_setting)))
    live = isnothing(output)
    if live && isempty(sd.name)
        error("a simulation run in a directory needs a name, that of its directory")
    end
    positive(s) = (n = parse(Int, s); n ≥ 1 ? n : nothing)
    samples = isnothing(sd.samples) ? nothing : launch_setting("TMS_SAMPLES", sd.samples, positive)
    check_destinations(checked_phases(sd), sd.final_measurements, !isnothing(samples); live)
    settings = (; max_time, checkpoint_interval, threading, output, samples)
    if live && !isnothing(samples)
        return run_parts(sd, settings)
    end
    dir = live ? abspath(sd.name) : ""
    src_path = Base.source_path()
    claim() = with_lock(() -> claim_run(dir, src_path), dir)
    sim, run = run_in(sd, dir, settings; claim) do sim, commit, threading
        sample, start, done, ended = isnothing(commit) ? (0, 1, 0, false) :
                                     (commit.sample, commit.phase, commit.step, commit.last)
        return run_simulation(sim, sd, samples, sample, start, done, ended, threading)
    end
    # what a program, or an AI running it, reads without opening the log
    @info summary_line(sd, run)
    return sim
end

"""
    empty_directory(dir)

remove what the directory `dir` holds, but its lock, which the caller holds, see `with_lock`
"""
function empty_directory(dir::String)
    for f in readdir(dir)
        if f ≠ "lock"
            rm(joinpath(dir, f); recursive = true, force = true)
        end
    end
    return nothing
end

"""
    check_program(dir, src_path, what)

whether the program `src_path` is accepted in place of the one that wrote `what` in the
directory `dir`, see `same_program`: refused unless `TMS_ACCEPT_PROGRAM=1`
"""
function check_program(dir::String, src_path, what::AbstractString)
    if same_program(dir, src_path)
        return false
    elseif get(ENV, "TMS_ACCEPT_PROGRAM", "") ≠ "1"
        error("$what of \"$dir\" was written by another program, or one given other arguments: " *
              "to go on with this one, copy it onto prog.jl, and the files it includes into " *
              "prog, or set TMS_ACCEPT_PROGRAM=1")
    end
    return true
end

# the program `src_path` and its arguments kept in the directory `dir`, see `same_program`
function store_program(dir::String, src_path)
    if isnothing(src_path) || src_path == ""
        return nothing
    end
    prog = joinpath(dir, "prog.jl")
    # cp refuses to copy the program onto itself, when it is the copy that is run
    if !(isfile(prog) && samefile(src_path, prog))
        cp(src_path, prog; force = true)
    end
    copy_program_files(dir, src_path)
    # no arguments, no file, as `same_program` reads it
    if isempty(ARGS)
        rm(joinpath(dir, "prog_args.json"); force = true)
    else
        write(joinpath(dir, "prog_args.json"), JSON.json(ARGS))
    end
    return nothing
end

"""
    claim_run(dir, src_path; accepted = nothing)

take the directory `dir` for a run of the program `src_path`, under its lock: refused while
another run may be using it, emptied for a file `restart`, its program checked unless `accepted`
says, and its file `running` written; `(; stale, resume, accepted)`, see `check_unlocked`
"""
function claim_run(dir::String, src_path; accepted::Union{Nothing, Bool} = nothing)
    stale = check_unlocked(dir)
    if isfile(joinpath(dir, "restart"))
        check_removable(dir)
        empty_directory(dir)
        stale = nothing
    end
    resume = has_checkpoint(dir)
    if isnothing(accepted)
        accepted = resume && check_program(dir, src_path, "the checkpoint")
    end
    write(joinpath(dir, "running"), "$(gethostname()) $(getpid()) $(Dates.now())\n")
    return (; stale, resume, accepted)
end

"""
    run_in(work, sd, dir, settings; claim, outdir = dir, main = "", seed = sd.seed,
           recording = false, clear_stop = true)

the run of `sd` in the directory `dir`, empty for none, with the settings of `runTMS`, the
directory taken by `claim`, see `claim_run`: `work(sim, commit, threading)` called with the
commit resumed or `nothing`; the simulation it returns and the `RunData`. A part of a simulation
with samples gives the directory of the simulation as `outdir` and `main`, see `Checkpointer`,
its seed, and `recording`, see `Recorder`, `clear_stop` false for the phases before its branch
point
"""
function run_in(work::Function, sd::SimData, dir::String, settings; claim = nothing,
                outdir::String = dir, main::String = "", seed = sd.seed, recording::Bool = false,
                clear_stop::Bool = true)
    (; max_time, checkpoint_interval, output, samples) = settings
    live = !isempty(dir)
    file(name) = joinpath(dir, name)
    outermost = enter_run()
    # only the outermost run sets the threading of the process
    threading = outermost ? settings.threading : nothing
    # whether the run has taken its directory: only then does the way out write `error` and
    # remove `running`
    claimed = false
    started = false
    logio = devnull
    run = nothing
    try
        if !isnothing(threading)
            set_threading(threading == :auto ? :dense : threading)
        end
        stale = nothing
        resume = false
        accepted = false
        if live
            claim_directory(dir)
            claimed = true
            mkpath(dir)
            (; stale, resume, accepted) = claim()
            started = true
            # left over from the run before
            if clear_stop
                rm(file("stop"); force = true)
            end
            rm(file("error"); force = true)
            # a resumed run continues the log, keeping the history of every run, as a part does
            logio = open(file("log"), resume || !isempty(main) ? "a" : "w")
            if sd.description ≠ ""
                write(file("description"), sd.description)
            end
            store_program(dir, Base.source_path())
        else
            logio = output
        end
        commit, generation, seed = resume ? load_checkpoint(dir, isempty(main) ? dir : main) :
                                   (nothing, 0, something(seed, rand(RandomDevice(), UInt64)))
        c = Checkpointer(dir, seed; interval = checkpoint_interval, max_time, main)
        outputs = Outputs(; dir = outdir, sd.format, sampled = !isnothing(samples), recording)
        run = RunData(outdir, outputs, logio, output, c, sd.provided, sd.backend)
        if resume
            resume!(c, commit, generation)
            run.phase_time = commit.phase_time
        end
        if live
            write(file("stamp"), """
                Julia $VERSION
                TensorMixedStates $(pkgversion(TensorMixedStates))
                Date $(Dates.now())
                Host $(gethostname())
                Seed $seed
                Max time $max_time
                Checkpoint interval $checkpoint_interval
                """ * packages_stamp() * threading_stamp(threading))
        end
        sim = if resume
            Simulation(commit.state, commit.time, copy(commit.values), run)
        else
            values = Dict{Symbol, Any}(:rng => Xoshiro(seed), :step => 0)
            if !isnothing(samples)
                values[:sample] = 0
            end
            Simulation(nothing, sd.time_start, values, run)
        end
        # each launch, so that the log of a run resumed elsewhere says where each part ran
        log_line(run, (resume ? "\n" : "") * "Run on $(gethostname()), $(Dates.now()), Julia " *
                      "$VERSION, TensorMixedStates $(pkgversion(TensorMixedStates))")
        if sd.backend !== State
            log_line(run, "States created as $(nameof(sd.backend)) by default")
        end
        if !isnothing(stale)
            log_line(run, "Warning: the run before, process $(stale[2]) on $(stale[1]), ended " *
                          "without removing running")
        end
        if isnothing(samples) && !isempty(get(ENV, "TMS_SAMPLES", ""))
            log_line(run, "TMS_SAMPLES is not used: the simulation has no samples")
        end
        if accepted
            log_line(run, "Warning: " * (isempty(main) ? "resumed by a program other than the " *
                                                         "one that wrote the checkpoint" :
                                         "run by a program other than the one that started the " *
                                         "samples") * ", as TMS_ACCEPT_PROGRAM=1 allows")
        end
        if !resume
            commit!(sim, 1, 0)
        end
        sim = Logging.with_logger(SimLogger(run, Logging.current_logger())) do
            if resume
                log_line(run, "Resuming from the checkpoint: " *
                              (commit.sample == 0 ? "" : "sample $(commit.sample), ") *
                              "phase $(commit.phase), step $(commit.step), time " *
                              time_text(run, commit.time))
                # before anything is written: a destination opened anew would empty a file the
                # checkpoint continues
                reopen!(outputs, commit.marks)
            end
            launched, allocated = time(), Base.gc_bytes()
            try
                sim = work(sim, commit, threading)
                report_dropped(run)
            catch e
                if !(e isa StopRun || (e isa InterruptException && live))
                    rethrow()
                end
                report_dropped(run)
                if e isa InterruptException
                    log_line(run, "\n***** Interrupted, writing a checkpoint *****")
                    logged_checkpoint(run)
                end
                log_line(run, stop_line(c, run.phases))
                run.stopped = true
                # a nested run passes the interrupt on, so that the run it is in stops as well
                if e isa InterruptException && !outermost
                    rethrow()
                end
                # what the run resumes from, its last commit
                k = c.last
                sim = Simulation(k.state, k.time, k.values, run)
            end
            # this launch only, a resumed run counting from its own start
            log_line(run, "\n***** This launch took $(round(time() - launched; digits = 3)) " *
                          "seconds of wall clock, " *
                          "$(Base.format_bytes(Base.gc_bytes() - allocated)) " *
                          "allocated, peak memory $(Base.format_bytes(Sys.maxrss())) *****")
            return sim
        end
        return (sim, run)
    catch e
        if started
            if !isnothing(run)
                report_dropped(run)
            end
            write(file("error"), error_line(run, sd, e) * "\n")
            println(logio, "\n" * error_line(run, sd, e))
            flush(logio)
        end
        rethrow()
    finally
        if !isnothing(run)
            close!(run.outputs)
        end
        if live && logio !== devnull
            close(logio)
        end
        if started
            rm(file("running"); force = true)
        end
        if claimed
            release_directory(dir)
        end
        leave_run()
    end
end

"""
    live_run(dir)
    live_parts(dir)

whether the file `running` of the directory `dir` marks a run that may still be going on, see
`live_marker`; and the numbers of the parts of the simulation of directory `dir` so marked
"""
live_run(dir::String) = isfile(joinpath(dir, "running")) &&
                        live_marker(read(joinpath(dir, "running"), String))

function part_numbers(dir::String)
    parts = joinpath(dir, "parts")
    numbers = isdir(parts) ? [ tryparse(Int, f) for f in readdir(parts) ] : Union{Nothing, Int}[]
    return sort!(Int[ p for p in numbers if !isnothing(p) ])
end

live_parts(dir::String) = filter(p -> live_run(part_dir(dir, p)), part_numbers(dir))

"""
    launch_samples!(dir, samples, src_path; first = true)

under the lock of the directory `dir` of a simulation with samples, a process joining it: on its
first call, the directory emptied for a file `restart`, and its file `stop` and the samples failed
forgotten when no other process runs it; the samples asked for set in the queue, and the role of
the process, with what goes with it: `:part` once the branch point is written, with whether the
program was accepted, `:wait` while another process runs the phases before it, with that
process, `:stop` for a file `stop` found while waiting, `:prefix` for this one to run them, the
directory taken, see `claim_run`
"""
function launch_samples!(dir::String, samples::Int, src_path; first::Bool = true,
                         merging = (format = "%.13g", tolerance = [], min_samples = 10))
    prefix = live_run(dir)
    if !first && isfile(joinpath(dir, "stop"))
        return :stop => nothing
    elseif first && !prefix && isempty(live_parts(dir))
        if isfile(joinpath(dir, "restart"))
            check_removable(dir)
            empty_directory(dir)
        end
        rm(joinpath(dir, "stop"); force = true)
        if isfile(queue_file(dir))
            q = read_queue(dir)
            empty!(q.failed)
            write_queue(dir, q)
        end
    elseif first && isfile(joinpath(dir, "restart"))
        error("$dir holds the file restart while other processes run the simulation: let them " *
              "end, or stop them by the file stop, first")
    end
    branched = isfile(joinpath(dir, "branch.h5"))
    accepted = branched && check_program(dir, src_path, "the branch point")
    q = isfile(queue_file(dir)) ? read_queue(dir) :
        Queue(samples, nothing, 1, Int[], Dict{Int, Int}(), Dict{Int, String}())
    q.samples = samples
    if branched
        # the tolerance of this program, on the samples merged
        m, marks = read_merged(dir)
        (; format, tolerance, min_samples) = merging
        q.reached = !isempty(tolerance) && within(marks, dir, format, tolerance, min_samples) ?
                    m : nothing
    end
    write_queue(dir, q)
    if branched
        store_program(dir, src_path)
        if accepted
            open(joinpath(dir, "log"), "a") do io
                println(io, "Warning: samples run by a program other than the one that started " *
                            "them, as TMS_ACCEPT_PROGRAM=1 allows, from $(Dates.now())")
            end
        end
        return :part => accepted
    elseif prefix
        return :wait => split(read(joinpath(dir, "running"), String))
    end
    return :prefix => claim_run(dir, src_path)
end

"""
    claim_part(dir, src_path, accepted)

under the lock of the directory `dir` of a simulation with samples, the part a process runs, its
directory taken, see `claim_run`: an idle part holding a sample, else, while samples are left,
an idle part or a new one; `nothing` when none is left. The checkpoint of an idle part not of
the sample it holds, of a sample ended, is removed.
"""
function claim_part(dir::String, src_path, accepted::Bool)
    q = read_queue(dir)
    if !isempty(failed_samples(q))
        throw(failed_error(q))
    end
    numbers = part_numbers(dir)
    idle = filter(p -> !live_run(part_dir(dir, p)), numbers)
    holding = filter(p -> haskey(q.current, p), idle)
    part = !isempty(holding) ? first(holding) :
           isnothing(next_sample(q)) ? nothing :
           !isempty(idle) ? first(idle) : maximum(numbers; init = 0) + 1
    if isnothing(part)
        return nothing
    end
    pdir = part_dir(dir, part)
    mkpath(pdir)
    if has_checkpoint(pdir) &&
       JSON.parsefile(checkpoint_json(pdir))["sample"] ≠ get(q.current, part, 0)
        for f in checkpoint_files[1:4]
            rm(joinpath(pdir, f); force = true)
        end
    end
    return part => claim_run(pdir, src_path; accepted)
end

"""
    samples_line(dir, merged; tolerance = false)

what the console takes of the samples of the simulation of directory `dir` once a process
returns: the samples merged, out of those asked for, whether a `tolerance` was reached, and the
other processes running
"""
function samples_line(dir::String, merged::Int; tolerance::Bool = false)
    n = length(live_parts(dir))
    q = read_queue(dir)
    asked = q.samples
    count = merged ≤ asked ? "$merged of $asked samples merged" :
            "$merged samples merged, more than the $asked asked for"
    return (!isnothing(q.reached) ? "$merged samples merged, within the tolerance" :
            tolerance && merged ≥ asked ? "$count, the tolerance not reached" : count) *
           (n == 0 ? "" : ", $n other process$(n == 1 ? "" : "es") running")
end

"""
    run_parts(sd, settings)

run the simulation with samples `sd` in its directory alongside the other processes running it:
the phases before the branch point once, the others waiting, see `run_prefix`, then parts taken
one after the other, see `run_part`, until no sample is left. It returns the simulation of the
last part, its `Data` those merged, or the branch point when no part was taken.
"""
function run_parts(sd::SimData, settings)
    dir = abspath(sd.name)
    # another run of this process would take its marks for those a dead process left
    if lock(() -> any(d -> d == dir || startswith(d, joinpath(dir, "")), active_runs.dirs),
            active_runs.lock)
        error("$dir is the directory of another run of this process, going on: give each run " *
              "a name of its own")
    end
    src_path = Base.source_path()
    mkpath(dir)
    first = true
    role, claimed = :wait, nothing
    while role == :wait
        role, claimed = with_lock(dir) do
            launch_samples!(dir, settings.samples, src_path; first, merging = merging(sd))
        end
        if role == :wait
            if first
                host, pid = length(claimed) ≥ 2 ? claimed[1:2] : ("?", "?")
                @info "simulation \"$(sd.name)\": waiting for the phases before the branch " *
                      "point, which process $pid on $host runs; if it runs no more, remove " *
                      joinpath(dir, "running")
            end
            first = false
            sleep(1)
        end
    end
    b = nothing
    if role == :stop
        @info "simulation \"$(sd.name)\" stopped while waiting for the branch point, by the " *
              "file stop"
    elseif role == :prefix
        sim, run = run_in(sd, dir, settings; claim = () -> claimed,
                          clear_stop = false) do sim, commit, threading
            return run_prefix(sim, sd, commit, threading)
        end
        if run.stopped
            @info summary_line(sd, run)
            return sim
        end
        claimed = false
    end
    if role ≠ :stop
        # what a process that died while merging left
        merge_all!(dir, merging(sd))
        b, seed = read_branch(dir)
    end
    sim, run = nothing, nothing
    while role ≠ :stop
        taken = with_lock(() -> claim_part(dir, src_path, claimed), dir)
        if isnothing(taken)
            break
        end
        part, claim = taken
        sim, run = run_in(sd, part_dir(dir, part), settings; claim = () -> claim, outdir = dir,
                          main = dir, seed, recording = true) do sim, commit, threading
            return run_part(sim, sd, b, commit, dir, part, threading)
        end
        if run.stopped
            break
        end
    end
    ran = !isnothing(run)
    if !ran
        run = RunData(dir, Outputs(; dir, sd.format, sampled = true), devnull, nothing,
                      Checkpointer("", isnothing(b) ? something(sd.seed, 0) : seed), sd.provided,
                      sd.backend)
        run.stopped = role == :stop
        sim = isnothing(b) ? Simulation(nothing, sd.time_start, Dict{Symbol, Any}(), run) :
              Simulation(b.state, b.time, copy(b.values), run)
    end
    if !isnothing(b)
        merged, marks = read_merged(dir)
        run.outputs = data_outputs(dir, sd.format, marks)
        line = samples_line(dir, merged; tolerance = !isempty(sd.tolerance))
        @info ran ? summary_line(sd, run, line) :
              "simulation \"$(sd.name)\" has no sample left to run, results in $dir; $line"
    end
    return sim
end
