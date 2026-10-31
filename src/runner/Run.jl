# runTMS, which runs a simulation described by a SimData phase after phase, in a directory of its
# own, with its log, its checkpoints and the resumption of an interrupted run, and cleanTMS, which
# removes that directory.

export runTMS, cleanTMS, SimData

"""
    flatten_phases(phases)

the phases given as nested vectors, flattened into a single list, where a phase has the position
a checkpoint records
"""
flatten_phases(p::AbstractVector) = reduce(vcat, map(flatten_phases, p); init = Any[])
flatten_phases(p) = Any[p]

"""
    check_phase(phase)

refuse an object that is not a subtype of `AbstractPhase`, or lacks one of the fields `name`,
`time_start` and `final_measurements`
"""
function check_phase(phase)
    if !(phase isa AbstractPhase)
        error("$(typeof(phase)) is not a phase, a subtype of TensorMixedStates.AbstractPhase")
    end
    for f in (:name, :time_start, :final_measurements)
        if !hasfield(typeof(phase), f)
            error("$(typeof(phase)) is not a phase, it has no field $f")
        end
    end
    return phase
end

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
    append!(lists, [ q.measurements for q in phases if hasfield(typeof(q), :measurements) ])
    return [ measurement_lists(l) for l in lists ]
end

"""
    refused_reads(phase)

the values the lists of measurements of `phase` cannot read, each with the reason given when one
does: none by default
"""
refused_reads(::AbstractPhase) = ()

"""
    check_values(phases, final_measurements)

refuse a value measured by one of the lists of the phases, or by the final measurements, that no
phase provides, see `provided_values`, or that its phase refuses, see `refused_reads`, and
return the names of those they provide
"""
function check_values(phases::Vector, final_measurements)
    for q in phases, (value, why) in refused_reads(q)
        lists = Any[q.final_measurements]
        if hasfield(typeof(q), :measurements)
            push!(lists, q.measurements)
        end
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
                if r ∉ provided && r ∉ (:step, :rng)
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
        names = hasfield(typeof(p), :on) && p.on isa Union{Symbol, Vector{Symbol}} ?
                state_names(p.on) : Symbol[]
        for n in [ names ; collect(provided_values(p)) ]
            if n ≠ :state && n in reserved_values
                error("$(nameof(typeof(p))) cannot name a state $n, a name the run reserves")
            end
        end
        for n in names
            if count(==(n), names) > 1
                error("$(nameof(typeof(p))) names the state $n twice in its on, which would act " *
                      "on it twice")
            end
        end
    end
end

"""
    check_destinations(phases, final_measurements)

refuse, before a run with a directory starts rather than when it writes there, a file of the
simulation directory as the destination of a list of measurements of the phases or of the final
measurements, or as a file a phase writes, see `phase_files` and `simulation_files`, and a text
file given two lists, as one holds one list, the lists one call sends it making one
"""
function check_destinations(phases::Vector, final_measurements)
    for p in phases, f in phase_files(p)
        if is_simulation_file(f)
            error("cannot save to $f, a file of the simulation directory: choose another name")
        end
    end
    held = Dict{String, Vector{String}}()
    for l in run_lists(phases, final_measurements), (d, m) in by_destination(Outputs(), l)
        if !(d isa Union{TextFile, JsonFile})
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
    phases::Vector{Any}
    provided::Vector{Symbol}
    function SimData(name, description, time_start, final_measurements, format, seed, phases)
        # which Julia 1.10 refuses
        if !isnothing(seed) && seed < 0
            error("the seed of a simulation is an integer of zero or above, got $seed")
        end
        ps = map(check_phase, flatten_phases(phases))
        provided = check_values(ps, final_measurements)
        check_names(ps)
        return new(name, description, time_start, final_measurements, check_format(format), seed,
                   ps, provided)
    end
end

SimData(; name::AbstractString = "simulation", description::AbstractString = "",
        time_start::Number = 0., final_measurements = [], format::AbstractString = "%.13g",
        seed::Union{Nothing, Integer} = nothing, phases) =
    SimData(name, description, time_start, final_measurements, format, seed, phases)

show(io::IO, s::SimData) =
    print(io, """
        SimData(
            name = $(repr(s.name)),
            description = $(repr(s.description)),
            time_start = $(s.time_start),
            final_measurements = $(s.final_measurements),
            format = $(repr(s.format)),
            seed = $(repr(s.seed)),
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
                    error("$(files[rel]) and $path would both be copied to prog/$rel: rename " *
                          "one of them")
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
    adapt_threading(sim, threading)

set the threading of the run before a phase: `:dense` or `:blocks` as asked, or for `:auto` the
mode suited to the system of the state, once there is one; a change is logged
"""
function adapt_threading(sim::Simulation, threading)
    mode = threading ≠ :auto ? threading :
           isnothing(sim.state) ? nothing : threading_mode(sim.state.system)
    if isnothing(mode)
        return nothing
    end
    before = threading_settings()
    set_threading(mode)
    after = threading_settings()
    if after != before
        log_message(sim, "Threading set to :$mode: BLAS threads $(after.blas), Strided threads " *
                         "$(after.strided), block sparse multithreading " *
                         (after.blocksparse ? "on" : "off"))
    end
    return nothing
end

"""
    run_phases(sim, phases, start, done, ended, threading)

run the phases from the phase `start`, its first `done` steps done, and the phase itself when
`ended`, see `Commit`: each logged, given its `time_start`, run by `run_phase` and measured by
its `final_measurements`, committed at its end, with a checkpoint when one is due, a stop asked
for, or after the last phase when the run writes checkpoints, so that a new run does nothing
"""
function run_phases(sim::Simulation, phases::Vector, start::Int, done::Int, ended::Bool, threading)
    run = getfield(sim, :run)
    c = run.checkpointer
    for i in start:length(phases)
        phase = phases[i]
        d = i == start ? done : 0
        run.phase = i
        # the lists of the phase before are measured no more
        empty!(run.measures)
        adapt_threading(sim, threading)
        if d == 0
            log_message(sim, "\n***** Starting phase \"$(phase.name)\" *****")
            if !isnothing(phase.time_start)
                sim = Simulation(sim, sim.state, phase.time_start)
            end
            run.phase_time = sim.time
            sim = with_values(sim; step = 0)
        else
            log_message(sim, "\n***** Resuming phase \"$(phase.name)\" after step $d *****")
        end
        run.step_start = time()
        td = @timed begin
            if !(i == start && ended)
                # called for the phase and the state as they are, see output
                sim = invokelatest(run_phase, phase, sim, d)
            end
            phase_output(sim, phase.final_measurements)
            sim
        end
        sim = td.value
        report_dropped(run)
        log_message(sim, "***** Ending phase \"$(phase.name)\" after " *
                         "$(round(td.time; digits = 3)) seconds, " *
                         "$(Base.format_bytes(td.bytes)) allocated, peak memory " *
                         "$(Base.format_bytes(Sys.maxrss())) *****")
        # absent before Julia 1.11; it includes the recompilation
        if haskey(td, :compile_time)
            log_message(sim, "compilation time was $(round(td.compile_time; digits = 3)) seconds")
        end
        commit!(sim, i + 1, 0)
        checkpoint_or_stop(run; force = i == length(phases) && (c.interval > 0 || c.generation ≠ 0))
    end
    return sim
end

"""
    stop_place(c, phases)

where a run stops, from its last commit
"""
function stop_place(c::Checkpointer, phases::Vector)
    k = c.last
    return k.phase > length(phases) ? "after the last phase" :
           k.step == 0 ? "before phase \"$(phases[k.phase].name)\"" :
           "after step $(k.step) of phase \"$(phases[k.phase].name)\""
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
    output_final(sim, sd)

the final measurements of the simulation, for `runTMS`, whose keyword `output` hides the function
`output` there
"""
output_final(sim::Simulation, sd::SimData) = phase_output(sim, sd.final_measurements)

"""
    summary_line(sd, run)

the line the console takes once a run returns: whether it ended or stopped, where its results
are, and how many warnings its log received
"""
function summary_line(sd::SimData, run::RunData)
    at = stop_place(run.checkpointer, sd.phases)
    state = !run.stopped ? "done" :
            isempty(run.dir) ? "stopped $at, with no directory to resume from" :
            "stopped $at: run it again to resume it"
    where = isempty(run.dir) ? "" : ", results in $(run.dir)"
    n = run.warnings
    warned = n == 0 ? "" : isempty(run.dir) ? ", $n warning$(n == 1 ? "" : "s") in its output" :
             ", $n warning$(n == 1 ? "" : "s") in its log"
    return "simulation \"$(sd.name)\" $state$where$warned"
end

"""
    error_line(run, sd, e)

the line the log and the file `error` take when the run fails with `e`: the phase it was in and
the time of its last commit
"""
function error_line(run, sd::SimData, e)
    k = isnothing(run) ? nothing : run.checkpointer.last
    where = isnothing(run) || run.phase == 0 ? "" :
            " in phase $(run.phase) \"$(sd.phases[run.phase].name)\", last committed at time " *
            time_text(run, k.time)
    return "ERROR$where: " * sprint(showerror, e)
end

"""
    runTMS(sim_data; max_time = Inf, checkpoint_interval = 0, threading = :dense, output = nothing)

run the simulation, see `SimData`, in a directory named after it, and return the `Simulation`
it ends with. A checkpoint found in the directory is resumed from, if it was written by the same
program, `prog.jl` in the directory, including the same files, copied into its directory `prog`,
given the same arguments: copying an edited program there, or setting `TMS_ACCEPT_PROGRAM=1`,
resumes it with that one. A file `restart` in the directory
starts the simulation over, and `cleanTMS` removes the directory. A destination that is a file of
the directory, as `log` or `stamp`, or a text file given two lists, is refused before the run
starts.

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
    if live
        check_destinations(sd.phases, sd.final_measurements)
    end
    dir = live ? abspath(sd.name) : ""
    file(name) = joinpath(dir, name)
    outermost = enter_run()
    # only the outermost run sets the threading of the process
    threading = outermost ? threading : nothing
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
        if live
            claim_directory(dir)
            claimed = true
            if isdir(dir)
                stale = check_unlocked(dir)
                if isfile(file("restart"))
                    check_removable(sd.name)
                    rm(dir; recursive = true)
                    stale = nothing
                end
            end
            mkpath(dir)
            src_path = Base.source_path()
            resume = has_checkpoint(dir)
            if resume && !same_program(dir, src_path) && get(ENV, "TMS_ACCEPT_PROGRAM", "") ≠ "1"
                error("the checkpoint of \"$(sd.name)\" was written by another program, or one " *
                      "given other arguments: to resume with this one, copy it onto prog.jl, and the " *
                      "files it includes into prog, or set TMS_ACCEPT_PROGRAM=1")
            end
            started = true
            write(file("running"), "$(gethostname()) $(getpid()) $(Dates.now())\n")
            # left over from the run before
            rm(file("stop"); force = true)
            rm(file("error"); force = true)
            # a resumed run continues the log, keeping the history of every run
            logio = open(file("log"), resume ? "a" : "w")
            if sd.description ≠ ""
                write(file("description"), sd.description)
            end
            # cp refuses to copy the program onto itself, when it is the copy that is run
            if !isnothing(src_path) && src_path ≠ ""
                if !(isfile(file("prog.jl")) && samefile(src_path, file("prog.jl")))
                    cp(src_path, file("prog.jl"); force = true)
                end
                copy_program_files(dir, src_path)
                # no arguments, no file, as `same_program` reads it
                if isempty(ARGS)
                    rm(file("prog_args.json"); force = true)
                else
                    write(file("prog_args.json"), JSON.json(ARGS))
                end
            end
        else
            logio = output
        end
        commit, generation, seed = resume ? load_checkpoint(dir) :
                                   (nothing, 0, something(sd.seed, rand(RandomDevice(), UInt64)))
        c = Checkpointer(dir, seed; interval = checkpoint_interval, max_time)
        outputs = Outputs(; dir, sd.format)
        run = RunData(dir, outputs, logio, output, c, sd.provided)
        if resume
            c.last = commit
            c.written = commit
            c.generation = generation
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
            Simulation(nothing, sd.time_start, Dict{Symbol, Any}(:rng => Xoshiro(seed), :step => 0),
                       run)
        end
        # each launch, so that the log of a run resumed elsewhere says where each part ran
        log_line(run, (resume ? "\n" : "") * "Run on $(gethostname()), $(Dates.now()), Julia " *
                      "$VERSION, TensorMixedStates $(pkgversion(TensorMixedStates))")
        if !isnothing(stale)
            log_line(run, "Warning: the run before, process $(stale[2]) on $(stale[1]), ended " *
                          "without removing running")
        end
        if !resume
            commit!(sim, 1, 0)
        end
        start, done, ended = resume ? (commit.phase, commit.step, commit.last) : (1, 0, false)
        sim = Logging.with_logger(SimLogger(run, Logging.current_logger())) do
            if resume
                log_line(run, "Resuming from the checkpoint: phase $(commit.phase), step " *
                              "$(commit.step), time $(time_text(run, commit.time))")
                # before anything is written: a destination opened anew would empty a file the
                # checkpoint continues
                reopen!(outputs, commit.marks)
            end
            try
                sim = run_phases(sim, sd.phases, start, done, ended, threading)
                output_final(sim, sd)
                report_dropped(run)
            catch e
                if !(e isa StopRun || (e isa InterruptException && live))
                    rethrow()
                end
                report_dropped(run)
                if e isa InterruptException
                    log_line(run, "\n***** Interrupted, writing a checkpoint *****")
                    write_checkpoint(c)
                end
                log_line(run, stop_line(c, sd.phases))
                run.stopped = true
                # a nested run passes the interrupt on, so that the run it is in stops as well
                if e isa InterruptException && !outermost
                    rethrow()
                end
                # what the run resumes from, its last commit
                k = c.last
                sim = Simulation(k.state, k.time, k.values, run)
            end
            return sim
        end
        # what a program, or an AI running it, reads without opening the log
        @info summary_line(sd, run)
        return sim
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
