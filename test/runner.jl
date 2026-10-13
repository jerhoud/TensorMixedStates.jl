# Simulations as runTMS runs them: SimData and the phases, the outputs and the log of a run, its
# checkpoints and its resumes, its lock, its program, the settings of a launch and the threading.
#
# Goes here: everything of src/runner. The tests run in temporary directories, which they make
# the working directory while they run, safe as long as the groups stay sequential.

using TensorMixedStates: AbstractPhase

"""
    request_stop()

ask the run under way to stop, as a user does by creating the file `stop` in its directory: the
directory of the working directory that holds the marker `running`, a measurement not being
told the directory of its simulation
"""
function request_stop()
    for d in readdir()
        if isfile(joinpath(d, "running"))
            touch(joinpath(d, "stop"))
        end
    end
end

"""
    stopper(n)

a measurement asking the run to stop the `n[]`-th time it is taken from now, `n[] = 0` never;
it gives 0 armed or not, so that every run writes the same columns
"""
stopper(n::Ref{Int}) = "stopper" => function (_)
    if n[] > 0
        n[] -= 1
        if n[] == 0
            request_stop()
        end
    end
    return 0.0
end

"""
    in_tmp(f)

`f()` with a temporary directory as working directory
"""
in_tmp(f) = mktempdir(dir -> cd(f, dir))

"""
    run_through(sd; kwargs...)

the simulation `sd` run again and again until it is not stopped
"""
function run_through(sd; kwargs...)
    for _ in 1:500
        sim = runTMS(sd; kwargs...)
        if !stopped(sim)
            return sim
        end
    end
    error("the simulation did not end")
end

# a phase of one's own made of steps: kicks on the first qubit, each taking 0.1 of time
Base.@kwdef struct Kicks <: AbstractPhase
    name::String = "kicks"
    time_start = nothing
    final_measurements = []
    nkicks::Int = 4
    measurements = []
end

TMS.prepare_phase(p::Kicks, _) = p
TMS.step_count(p::Kicks) = p.nkicks
TMS.run_step(::Kicks, sim, _) = Simulation(sim, apply(exp(-0.3im * X)(1), sim.state), sim.time + 0.1)

# what is not a phase, lacking its fields, and a phase that has no way to run
struct NoFields <: AbstractPhase end

struct Unrunnable <: AbstractPhase
    name::String
    time_start
    final_measurements
end

# a phase of one's own of a single call, summing what it reads after each kick, which it
# provides as a value
Base.@kwdef struct SummingKicks <: AbstractPhase
    name::String = "summing kicks"
    time_start = nothing
    final_measurements = []
    nkicks::Int = 4
    measurements = []
end

TMS.provided_values(::SummingKicks) = (:total,)

function TMS.run_phase(p::SummingKicks, sim::Simulation, done::Int)
    total = done == 0 ? 0.0 : sim.total
    for k in done + 1:p.nkicks
        st = apply(exp(-0.3im * X)(1), sim.state)
        total += real(expect(st, Z(1)))
        sim = TMS.step_done(Simulation(sim, st; total), p, k)
    end
    return sim
end

@testset "SimData" begin
    phases = [CreateState{Pure}(2, Qubit(), "Up"), Gates(gates = X(1))]
    # nested vectors are flattened
    @test length(SimData(phases = [phases[1], [[phases[2]], phases[2]]]).phases) == 3
    # the first phase creates the state, and every phase is one
    @test_throws "must create the state" SimData(phases = [Gates(gates = X(1))])
    @test_throws "is not a phase" SimData(phases = [phases; 3])
    @test_throws "is not a phase" SimData(phases = [phases; SimData(; phases)])
    @test_throws "has no field name" SimData(phases = [phases; NoFields()])
    # a value measured that no phase provides
    @test_throws "reads energy" SimData(; phases, final_measurements = "f" => :energy)
    @test_throws "reads ground" SimData(; phases, final_measurements = "f" => Fidelity(:ground))
    @test_ok SimData(phases = [phases; CopyState(:state => :ground)],
                     final_measurements = "f" => [Fidelity(:ground), :step, "r" => c -> rand(c.rng)])
    @test_ok SimData(phases = [phases; GroundState(hamiltonian = Z(1), algo = Dmrg(nsweeps = 1))],
                     final_measurements = "f" => :energy)
    @test_throws "does not write one number" SimData(; phases, format = "%s %s")
    # the reserved names
    @test_throws "a name the run reserves" CopyState(:state => :rng)
    @test_throws "under another name" CopyState(:a => :a)
    # a phase prints the values of its fields, those of a phase of one's own too
    @test occursin("nkicks = 4", repr(Kicks()))
end

@testset "A run stopped at every step resumes as the uninterrupted one" begin
    # a deadline already past stops every run after one step, or one phase without steps: the
    # outputs, the draws of the generator of the run included, come out the same, and the state
    # up to the rounding two uninterrupted runs differ by
    in_tmp() do
        H = -sum(Z(i) * Z(i + 1) for i in 1:3) - 0.8 * sum(X(i) for i in 1:4)
        lim = Limits(maxdim = 8, cutoff = 1e-14)
        draw = "draw" => c -> rand(c.rng)
        # the imaginary part of the trace is rounding, which differs from one run to another as
        # the state does
        tr = "trace" => c -> real(trace(c.state))
        phases = [CreateState{Pure}(4, Qubit(), "Up"; randomize = 4),
                  CopyState(:state => :start),
                  Evolve(evolver = [-im * H, -im * X(2)] => [t -> 1, t -> cos(3t)], duration = 0.3,
                         time_step = 0.1, limits = lim,
                         measurements = "evolve" => [X(1), Z, :step, draw, Fidelity(:start)]),
                  GroundState(hamiltonian = H, algo = Dmrg(nsweeps = 30, tol = 1e-10), limits = lim,
                              measurements = "search" => [:energy, :step, Z(1), draw],
                              final_measurements = "found" => [:energy, Fidelity(:start)]),
                  ToMixed(),
                  Evolve(evolver = -im * H + Dissipator(sqrt(0.3) * Sm)(2), duration = 0.2,
                         time_step = 0.1, algo = ApproxW(order = 2), limits = lim,
                         measurements = "open" => [X(1), Z, Purity]),
                  Gates(gates = X(1)),
                  Thermalize(hamiltonian = H, beta = 0.4, beta_step = 0.1, limits = lim,
                             measurements = "thermo" => [:beta, :log_trace, Purity, draw]),
                  CopyState(:start => :state; move = true),
                  ToMixed(),
                  SteadyState(lindbladian = -im * H + Dissipator(Sm)(1) + Dissipator(Sp)(4),
                              algo = Dmrg(nsweeps = 4), limits = lim,
                              measurements = "steady" => [:energy, Z(1), tr])]
        sd(name) = SimData(; name, seed = 1, phases,
                           final_measurements = "end" => [tr, :energy, draw, Fidelity(:start)])
        ref = runTMS(sd("ref"))
        sim = run_through(sd("chk"); max_time = -1)
        for f in ("evolve", "search", "found", "open", "thermo", "steady", "end")
            @test read("chk/$f", String) == read("ref/$f", String)
        end
        @test norm(sim.state - ref.state) < 1e-12
        @test sim.energy ≈ ref.energy
        @test !haskey(getfield(sim, :values), :start)
        # the search ended on its tolerance, before its 30 sweeps
        @test 3 < length(readlines("ref/search")) < 31
        @test count("Resuming from the checkpoint", read("chk/log", String)) > 30
    end
end

@testset "Interrupts, kills and stops" begin
    in_tmp() do
        fail_in, crash_in, stop_in = Ref(0), Ref(0), Ref(0)
        breaker = "breaker" => function (_)
            if fail_in[] > 0
                fail_in[] -= 1
                if fail_in[] == 0
                    throw(InterruptException())
                end
            end
            if crash_in[] > 0
                crash_in[] -= 1
                if crash_in[] == 0
                    error("simulated kill")
                end
            end
            return 0.0
        end
        lim = Limits(maxdim = 4, cutoff = 1e-15)
        evolve(op) = Evolve(evolver = -im * op, duration = 0.3, time_step = 0.1, limits = lim,
                            measurements = "data" => [X(1), Y(1), breaker, stopper(stop_in)])
        phases = [CreateState{Pure}(3, Qubit(), "X+"), evolve(Z(1)), Gates(gates = Z(1)),
                  evolve(Z(2))]
        runTMS(SimData(; name = "ref", phases))
        reference = read("ref/data", String)

        # an interrupt writes the last commit and stops, as a stop does
        fail_in[] = 4
        sd = SimData(; name = "int", phases)
        @test stopped(runTMS(sd))
        @test fail_in[] == 0
        @test occursin("Interrupted, writing a checkpoint", read("int/log", String))
        @test !stopped(runTMS(sd))
        @test read("int/data", String) == reference

        # a run killed writes no checkpoint: it resumes from an older one, the lines written
        # since cut back, and marks the error with its phase
        stop_in[] = 2
        sd = SimData(; name = "kill", phases)
        runTMS(sd; checkpoint_interval = 1e9)
        checkpointed = length(readlines("kill/data"))
        crash_in[] = 3
        @test_throws "simulated kill" runTMS(sd; checkpoint_interval = 1e9)
        @test crash_in[] == 0
        @test length(readlines("kill/data")) > checkpointed
        @test occursin("in phase 4 \"Time evolution\"", read("kill/error", String))
        @test occursin("ERROR in phase 4", read("kill/log", String))
        @test !isfile("kill/running")
        runTMS(sd)
        @test read("kill/data", String) == reference
        # the marker describes the last run, which succeeded
        @test !isfile("kill/error")

        # an interrupt without a directory reaches the caller
        fail_in[] = 1
        @test_throws InterruptException runTMS(sd; output = devnull)
    end
end

@testset "A stop after the last phase, and a completed simulation run again" begin
    in_tmp() do
        stop_in = Ref(0)
        phases = [CreateState{Pure}(2, Qubit(), "Up"),
                  Evolve(evolver = -im * X(1), duration = 0.3, time_step = 0.1,
                         final_measurements = "data" => [Z(1), stopper(stop_in)])]
        sd(name) = SimData(; name, phases, final_measurements = "fin" => Z(1))
        ref = runTMS(sd("ref"))
        # the final measurements of the simulation, which the stopped run leaves, are taken
        # at the time reached
        stop_in[] = 1
        @test stopped(runTMS(sd("chk")))
        @test !isfile("chk/fin")
        sim = runTMS(sd("chk"))
        @test sim.time ≈ ref.time
        @test read("chk/fin", String) == read("ref/fin", String)
        @test read("chk/data", String) == read("ref/data", String)
        # the checkpoint a stop wrote is written again after the last phase, so that the
        # simulation run again does nothing, and its files stay as they are
        again = runTMS(sd("chk"))
        @test again.time ≈ ref.time
        @test read("chk/fin", String) == read("ref/fin", String)
        @test occursin("Resuming from the checkpoint: phase 3", read("chk/log", String))
        # and so with periodic checkpoints
        runTMS(sd("per"); checkpoint_interval = 1e9)
        runTMS(sd("per"); checkpoint_interval = 1e9)
        @test read("per/fin", String) == read("ref/fin", String)
        # without them, no checkpoint is written, and the simulation runs again whole
        runTMS(sd("none"))
        @test !isfile("none/checkpoint.json")
    end
end

@testset "The log, and the files of a resume" begin
    in_tmp() do
        stop_in = Ref(2)
        phases = [CreateState{Pure}(2, Qubit(), "Up"),
                  Kicks(measurements = "data" => [Z(1), :step, stopper(stop_in)],
                        final_measurements = "final" => Z(1))]
        sd = SimData(; name = "sim", phases)
        s = runTMS(sd)
        @test stopped(s)
        # a stopped run hands back what it resumes from, and takes no final measurement
        @test s.time ≈ 0.2
        @test !isfile("sim/final")
        # a file cut by hand is continued as it is, with a warning
        write("sim/data", "")
        @test !stopped(runTMS(sd))
        log = read("sim/log", String)
        @test occursin("WARNING: $(abspath("sim/data")) is shorter than at the checkpoint", log)
        # the log keeps the history of every run
        stop = findfirst("Stopping after step 2 of phase \"kicks\"", log)
        @test !isnothing(stop)
        @test first(stop) < first(findfirst("Resuming from the checkpoint: phase 2, step 2", log))
        @test occursin("Resuming phase \"kicks\" after step 2", log)
        @test length(readlines("sim/data")) == 2
    end
end

@testset "Phases of one's own" begin
    in_tmp() do
        # made of steps, and of a single call committing each of its steps, its value carried
        stop_in = Ref(0)
        phases = [CreateState{Pure}(2, Qubit(), "Up"),
                  Kicks(measurements = "kicks" => [Z(1), stopper(stop_in)]),
                  SummingKicks(measurements = "sums" => [:total, :step, stopper(stop_in)],
                               final_measurements = "total" => :total)]
        sd(name) = SimData(; name, phases)
        ref = runTMS(sd("ref"))
        @test ref.time ≈ 0.4
        sim = run_through(sd("chk"); max_time = -1)
        for f in ("kicks", "sums", "total")
            @test read("chk/$f", String) == read("ref/$f", String)
        end
        @test sim.total ≈ ref.total
        # one that has neither interface
        @test_throws "does not know how to run it" runTMS(
            SimData(name = "none", phases = [phases[1], Unrunnable("n", nothing, [])]))
    end
end

@testset "The program that resumes a checkpoint" begin
    # the same file, byte for byte, given the same arguments, or the consent of
    # TMS_ACCEPT_PROGRAM=1; the copy of the program run from the directory runs again
    program(state) = """
        using TensorMixedStates, .Qubits
        runTMS(SimData(name = "sim", description = "$state",
                       phases = [CreateState{Pure}(2, Qubit(), "$state"), Gates(gates = X(1))]);
               max_time = -1)
        """
    function with_args(f, args)
        saved = copy(ARGS)
        empty!(ARGS)
        append!(ARGS, args)
        try
            f()
        finally
            empty!(ARGS)
            append!(ARGS, saved)
        end
    end
    in_tmp() do
        run_prog(args = String[]) = with_args(() -> include(joinpath(pwd(), "prog.jl")), args)
        write("prog.jl", program("Up"))
        run_prog()
        @test !isfile("sim/prog_args.json")
        stamp = read("sim/stamp", String)
        # refused before anything is written
        write("prog.jl", program("Dn"))
        @test_throws "written by another program" run_prog()
        @test read("sim/description", String) == "Up"
        @test read("sim/stamp", String) == stamp
        @test read("sim/prog.jl", String) == program("Up")
        @test !isfile("sim/error") && !isfile("sim/running")
        # other arguments are another simulation
        write("prog.jl", program("Up"))
        @test_throws "written by another program" run_prog(["other"])
        # the consent, for one launch, copies the program there
        write("prog.jl", program("Dn"))
        withenv(() -> run_prog(), "TMS_ACCEPT_PROGRAM" => "1")
        @test read("sim/prog.jl", String) == program("Dn")
        @test read("sim/description", String) == "Dn"
        # the copy of the program runs again
        with_args(() -> include(joinpath(pwd(), "sim", "prog.jl")), String[])
        @test !isfile("sim/error")
        # a directory that lost its program is another program
        mv("sim/prog.jl", "saved.jl")
        @test_throws "written by another program" run_prog()
        mv("saved.jl", "sim/prog.jl")
        # the arguments are kept, a run without them refused, a fresh start without them
        # removing those of the run before
        rm("sim"; recursive = true)
        run_prog(["x"])
        @test TMS.JSON.parsefile("sim/prog_args.json") == ["x"]
        @test_throws "written by another program" run_prog()
        rm("sim/checkpoint.json")
        run_prog()
        @test !isfile("sim/prog_args.json")
        # with no program file, from the REPL, nothing is compared and the checkpoint resumed
        write("sim/prog.jl", "something else")
        task_local_storage(:SOURCE_PATH, nothing) do
            runTMS(SimData(name = "sim", phases = [CreateState{Pure}(2, Qubit(), "Dn"),
                                                   Gates(gates = X(1))]); max_time = -1)
        end
        @test !isfile("sim/error")
        @test read("sim/prog.jl", String) == "something else"
    end
end

@testset "A directory another run may be using is refused" begin
    phases = [CreateState{Pure}(2, Qubit(), "Up")]
    julia = `$(Base.julia_cmd()) --startup-file=no`
    in_tmp() do
        runTMS(SimData(; name = "sim", phases))
        @test !isfile("sim/running")
        # a process of this machine still there, another machine, or a file naming none
        p = run(`$julia -e "sleep(300)"`; wait = false)
        try
            for marker in ("$(gethostname()) $(getpid(p)) 2026-10-13T12:00:00\n",
                           "elsewhere 12 2026-10-13T12:00:00\n", "")
                write("sim/running", marker)
                @test_throws "may still be going on" runTMS(SimData(; name = "sim", phases))
                @test_throws "may still be going on" cleanTMS("sim")
                @test read("sim/running", String) == marker
                @test !isfile("sim/error")
            end
            write("sim/running", "elsewhere 12 2026-10-13T12:00:00\n")
            @test_throws "started 2026-10-13T12:00:00" runTMS(SimData(; name = "sim", phases))
        finally
            kill(p)
        end
        # a process of this machine that has ended, or this very one, whose number it had
        q = run(`$julia -e "sleep(1)"`; wait = false)
        ended = getpid(q)
        wait(q)
        for pid in (ended, getpid())
            write("sim/running", "$(gethostname()) $pid 2026-10-13T12:00:00\n")
            runTMS(SimData(; name = "sim", phases))
            @test !isfile("sim/running")
            @test occursin("process $pid on $(gethostname()), ended without removing running",
                           read("sim/log", String))
        end
    end
end

@testset "restart, cleanTMS and the working directory" begin
    in_tmp() do
        here = pwd()
        phases = [CreateState{Pure}(2, Qubit(), "Up"), SaveState(file = "s.h5")]
        sd = SimData(; name = "sim", phases)
        runTMS(sd)
        # the working directory is left as it is, the files written in the directory of the
        # simulation, LoadState reading there too
        @test pwd() == here
        @test isfile("sim/s.h5") && !isfile("s.h5")
        runTMS(SimData(name = "second", phases = [LoadState(file = "../sim/s.h5"), SaveState(file = "t.h5")]))
        @test isfile("second/t.h5")
        # a file restart starts the simulation over, and goes with the directory
        touch("sim/old")
        touch("sim/restart")
        runTMS(sd)
        @test !isfile("sim/old") && !isfile("sim/restart") && isfile("sim/log")
        # cleanTMS removes it
        cleanTMS(sd)
        @test !isdir("sim")
        @test_ok cleanTMS("absent")
        # neither removes the current directory
        mkpath("sub")
        touch("sub/precious")
        cd("sub") do
            for name in (".", "..", pwd())
                @test_throws "contains the current directory" cleanTMS(name)
            end
            mkpath("x")
            touch("x/restart")
            cd("x") do
                @test_throws "contains the current directory" runTMS(SimData(name = ".", phases = phases[1:1]))
            end
        end
        @test isfile("sub/precious")
        # a directory that cannot be made leaves the caller's alone
        write("clash", "")
        write("running", "mine")
        for (name, err) in (("clash", Base.IOError), ("", "needs a name"))
            @test_throws err runTMS(SimData(; name, phases = phases[1:1]))
            @test !isfile("error")
            @test read("running", String) == "mine"
            @test pwd() == here
        end
    end
end

@testset "The files of the simulation are its own" begin
    in_tmp() do
        start = CreateState{Pure}(2, Qubit(), "Up")
        for file in ("checkpoint-1.h5", "log", "./stamp")
            @test_throws "a file of the simulation directory" runTMS(
                SimData(name = "a", phases = [start, SaveState(; file)]))
            @test_throws "a file of the simulation directory" runTMS(
                SimData(name = "b", phases = [start], final_measurements = file => Z(1)))
        end
        # without a directory, every name goes to the stream
        @test_ok runTMS(SimData(name = "c", phases = [start], final_measurements = "log" => Z(1));
                        output = devnull)
        # two names of one file are one destination
        runTMS(SimData(name = "d", phases = [start, Gates(gates = X(1),
                                                         final_measurements = ["data" => X(1), "./data" => X(1)])]))
        @test length(readlines("d/data")) == 3
    end
end

@testset "States named in the simulation" begin
    in_tmp() do
        phases = [CreateState{Pure}(2, Qubit(), "Up"),
                  CopyState(:state => :up),
                  Gates(gates = X(1)),
                  SaveState(file = "states.h5", statename = "up", on = :up),
                  LoadState(file = "states.h5", statename = "up", on = :again),
                  CopyState(:again => :state; move = true)]
        sim = runTMS(SimData(; name = "sim", phases,
                             final_measurements = "f" => [Fidelity(:up), Overlap(:up), Z(1)]))
        @test real(expect(sim.state, Z(1))) ≈ 1
        @test sim.up isa State{Pure}
        @test !haskey(getfield(sim, :values), :again)
        @test parse.(Float64, split(readlines("sim/f")[2])[2:3]) ≈ [1, 1]
        # a state the simulation does not have
        @test_throws "has no state named" runTMS(SimData(name = "x", phases = [phases[1],
            SaveState(file = "s.h5", on = :nowhere)]))
        # a value under a reserved name
        @test_throws "a name the run reserves" Simulation(sim, sim.state; step = 1)
    end
end

@testset "The settings of a launch" begin
    in_tmp() do
        phases = [CreateState{Pure}(2, Qubit(), "Up"),
                  Kicks(measurements = "data" => Z(1))]
        sd = SimData(; name = "sim", phases)
        # an environment variable overrides the keyword, and the stamp records the settings
        @test stopped(withenv(() -> runTMS(sd), "TMS_MAX_TIME" => "-1"))
        @test occursin("Max time -1.0", read("sim/stamp", String))
        @test_throws "TMS_MAX_TIME is \"soon\"" withenv(() -> runTMS(sd), "TMS_MAX_TIME" => "soon")
        @test_throws "TMS_THREADING is \"many\"" withenv(() -> runTMS(sd), "TMS_THREADING" => "many")
        @test_throws "unknown threading" runTMS(sd; threading = :many)
        withenv("TMS_THREADING" => "nothing", "TMS_CHECKPOINT_INTERVAL" => "1e9") do
            run_through(sd)
        end
        stamp = readlines("sim/stamp")
        @test "Threading nothing" in stamp && "Checkpoint interval 1.0e9" in stamp
        @test "Julia threads $(Threads.nthreads())" in stamp && "CPU threads $(Sys.CPU_THREADS)" in stamp
        @test any(l -> startswith(l, "BLAS ") && endswith(l, " threads"), stamp)
        @test any(l -> startswith(l, "Block sparse multithreading "), stamp)
        # the seed is written, and drawn when none is given
        seed = only(l for l in stamp if startswith(l, "Seed "))
        runTMS(SimData(name = "other", phases = phases[1:1]))
        @test only(l for l in readlines("other/stamp") if startswith(l, "Seed ")) ≠ seed
        runTMS(SimData(name = "fixed", phases = phases[1:1], seed = 42))
        @test "Seed 42" in readlines("fixed/stamp")
    end
end

@testset "A run without a directory, and runs within a run" begin
    in_tmp() do
        phases = [CreateState{Pure}(2, Qubit(), "Up"), Kicks(measurements = "data" => Z(1))]
        io = IOBuffer()
        sim = runTMS(SimData(; name = "sim", phases, final_measurements = Data("d") => Z(1));
                     output = io)
        text = String(take!(io))
        @test occursin("Starting phase \"kicks\"", text) && occursin("Z(1)", text)
        @test !isdir("sim")
        @test length(sim.data["d"]) == 1
        # a deadline stops it, with nothing to resume
        @test stopped(runTMS(SimData(; name = "sim", phases); output = devnull, max_time = -1))
        # a run within a run, in a phase of one's own, leaves the threading and the behaviour
        # of Ctrl-C to the outermost
        start = ThreadingState()
        inner = Ref{Any}(nothing)
        nested = "nested" => function (_)
            inner[] = runTMS(SimData(name = "inner", phases = phases[1:1]); threading = :blocks)
            return 0.0
        end
        runTMS(SimData(name = "outer", phases = phases[1:1], final_measurements = "n" => nested))
        @test inner[] isa Simulation
        @test "Threading nothing" in readlines("inner/stamp")
        @test ThreadingState() == start
        # runs on several threads at once
        sims = fetch.([ Threads.@spawn runTMS(SimData(; name = "t$k", phases); output = devnull)
                        for k in 1:4 ])
        @test all(s -> s.time ≈ 0.4, sims)
        @test TMS.active_runs.count[] == 0
    end
end

@testset "Schedules of the sweeps and of the steps across a resume" begin
    # a resumed search takes the rest of its schedules, and an evolution reads them at the step
    # it runs: a schedule started over would truncate a state the uninterrupted run kept
    in_tmp() do
        H = sum(-Z(i) * Z(i + 1) for i in 1:5) - sum(1.0 * X(i) for i in 1:6)
        ms(f) = f => [MaxLinkdim, X(1), Z(1) * Z(2), :step]
        phases = [CreateState{Pure}(6, Qubit(), "Z+"),
                  GroundState(hamiltonian = H, algo = Dmrg(nsweeps = 4, noise = [1e-2, 1e-3, 0]),
                              limits = Limits(maxdim = [2, 2, 8, 8], cutoff = 1e-14),
                              measurements = ms("search"), measurements_period = 2),
                  CreateState{Pure}(6, Qubit(), "X+"),
                  Evolve(evolver = -im * H, duration = 0.4, time_step = 0.1,
                         limits = Limits(maxdim = [2, 2, 8, 8], cutoff = 1e-14),
                         measurements = ms("evolve"))]
        runTMS(SimData(; name = "ref", phases))
        # the bond dimension follows the schedule, else the comparison below would say nothing
        @test length(unique(l -> split(l)[2], readlines("ref/evolve")[2:end])) > 1
        run_through(SimData(; name = "chk", phases); max_time = -1)
        @test read("chk/search", String) == read("ref/search", String)
        @test read("chk/evolve", String) == read("ref/evolve", String)
    end
end

"""
    quietly(f)

`f()`, what it prints thrown away: ITensors warns when block sparse multithreading is switched
on while Julia has a single thread
"""
quietly(f) = redirect_stdout(f, devnull)

@testset "ThreadingState and set_threading" begin
    now = threading_settings()
    s = ThreadingState()
    @test (s.blas, s.strided, s.blocksparse) == (now.blas, now.strided, now.blocksparse)
    # a field left out takes the value in force, and it prints as the call that builds it
    @test ThreadingState(blas = 3) == ThreadingState(3, s.strided, s.blocksparse)
    t = ThreadingState(blas = 2, strided = 1, blocksparse = true)
    @test repr(t) == "ThreadingState(blas = 2, strided = 1, blocksparse = true)"
    @test eval(Meta.parse(repr(t))) == t
    start = ThreadingState()
    try
        # every change returns the threading it replaces, which set_threading takes back
        @test set_threading(:dense) == start
        dense = ThreadingState()
        @test dense == ThreadingState(blas = start.blas, strided = 1, blocksparse = false)
        @test quietly(() -> set_threading(:blocks)) == dense
        @test ThreadingState() == ThreadingState(blas = 1, strided = 1, blocksparse = true)
        set_threading(dense)
        @test ThreadingState() == dense
        # :dense gives BLAS back the threads Julia starts it with after :blocks
        set_threading(ThreadingState(blas = 2))
        set_threading(:dense)
        @test threading_settings().blas == 2
        quietly(() -> set_threading(:blocks))
        set_threading(:dense)
        @test threading_settings().blas == TMS.default_blas_threads()
        # the mode suited to a system or a state: blocks when it conserves something and Julia
        # has threads to run them on
        set_threading(System(2, Qubit()))
        @test !ThreadingState().blocksparse
        quietly(() -> set_threading(State{Pure}(System(2, Fermion(conserve = N)), "Occ")))
        @test ThreadingState().blocksparse == (Threads.nthreads() > 1)
        @test_throws "unknown threading mode" set_threading(:blocksparse)
    finally
        quietly(() -> set_threading(start))
    end
end

@testset "The threading of a run" begin
    start = ThreadingState()
    phases = [CreateState{Pure}(2, Qubit(), "Up")]
    try
        in_tmp() do
            # a mode that does not depend on the state is set from the start, which the stamp
            # records, and the threading found is put back at the end
            runTMS(SimData(name = "dense"; phases))
            stamp = readlines("dense/stamp")
            @test "Threading :dense" in stamp && "Strided threads 1" in stamp
            @test ThreadingState() == start
            runTMS(SimData(name = "left"; phases); threading = nothing)
            @test "Threading nothing" in readlines("left/stamp")
            @test ThreadingState() == start
            # :auto chooses before each phase from the system of the state, and logs each change
            hopping = sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3)
            evolve = Evolve(evolver = -im * hopping, duration = 0.1, time_step = 0.1)
            auto = [CreateState{Mixed}(4, Fermion(conserve = N), ["Occ", "Emp", "Occ", "Emp"]),
                    evolve, Weaken(target = ()), evolve]
            quietly(() -> runTMS(SimData(name = "auto", phases = auto); threading = :auto))
            @test "Strided threads 1" in readlines("auto/stamp")
            changes = filter(l -> startswith(l, "Threading set to"), readlines("auto/log"))
            if Threads.nthreads() > 1
                @test length(changes) == 2
                @test startswith(changes[1], "Threading set to :blocks")
                @test startswith(changes[2], "Threading set to :dense")
            else
                @test isempty(changes)
            end
            @test ThreadingState() == start
        end
    finally
        quietly(() -> set_threading(start))
    end
end

@testset "An excited state of named states" begin
    # GroundState orthogonal to the states the simulation named, resumed after each sweep as the
    # uninterrupted run
    in_tmp() do
        H = sum(X(i) * X(i + 1) + Y(i) * Y(i + 1) + Z(i) * Z(i + 1) for i in 1:3) + 0.3 * sum(Z(i) for i in 1:4)
        algo = Dmrg(nsweeps = 30, tol = 1e-12)
        lim = Limits(maxdim = 16)
        phases = [CreateState{Pure}(4, Qubit(), "Up"; randomize = 4),
                  CopyState(:state => :start),
                  GroundState(hamiltonian = H, algo = algo, limits = lim),
                  CopyState(:state => :ground),
                  CopyState(:start => :state),
                  GroundState(hamiltonian = H, algo = algo, limits = lim, orthogonal_to = [:ground],
                              measurements = "excited" => [:energy, Z(1)])]
        sd(name) = SimData(; name, phases, seed = 3, final_measurements = "end" => [:energy, Z])
        ref = runTMS(sd("ref"))
        LA, sys = TMS.LinearAlgebra, ref.state.system
        E = LA.eigvals(LA.Hermitian(Matrix(TMS.system_matrix(Pure, sys, TMS.lower(sys, H)))))
        @test ref.energy ≈ E[2] atol = 1e-9
        @test abs2(inner(ref.ground, ref.state)) < 1e-12
        # the weight by default, 3 bonds of 3 terms of width 2 and 4 fields of width 0.6
        @test occursin("Searching an excited state, of weight 20.4", read("ref/log", String))
        run_through(sd("chk"); max_time = -1)
        @test read("chk/excited", String) == read("ref/excited", String)
        @test read("chk/end", String) == read("ref/end", String)
    end
end

@testset "Update" begin
    in_tmp() do
        # a value set, a state kept and removed, the time changed, each by an assignment
        phases = [CreateState{Pure}(2, Qubit(), "Up"),
                  Update(c -> c.ground = c.state; provides = [:ground]),
                  Gates(gates = X(1)),
                  Update(c -> c.x = real(expect(c.state, Z(1))); provides = [:x],
                         final_measurements = "kept" => [:x, Fidelity(:ground)]),
                  Update(c -> begin
                      c.ground = nothing
                      c.time = 5.
                  end; final_measurements = "removed" => [:x, Fidelity(:ground)])]
        sim = runTMS(SimData(; name = "sim", phases))
        @test sim.x ≈ -1
        @test !haskey(getfield(sim, :values), :ground)
        @test sim.time == 5.
        @test parse.(Float64, split(readlines("sim/kept")[2])) ≈ [0, -1, 0] atol = 1e-12
        @test split(readlines("sim/removed")[2], '\t')[3] == ""
        # a value measured that no Update provides, and the names the run reserves
        @test_throws "reads y" SimData(phases = [phases[1], Update(c -> c.y = 1)],
                                       final_measurements = "f" => :y)
        @test_throws "a name the run reserves" Update(c -> nothing; provides = [:rng])
        @test_throws "step cannot be changed" runTMS(SimData(name = "res", phases = [phases[1],
                                                                 Update(c -> c.step = 3)]))
        # a failure leaves the simulation without what the function had changed
        @test_throws "half done" runTMS(SimData(name = "fail", phases = [phases[1],
            Update(c -> begin
                c.state = apply(X(1), c.state)
                error("half done")
            end)]))
        # step by step: kicks measured at random every second step, output and log_message
        # within, resumed after each step as the uninterrupted run
        kick = exp(-0.4im * X)(1)
        steps = Update(nsteps = 6, provides = [:m], measurements = "kicks" => [:m, :step, Z(1)],
                       measurements_period = 2) do c
            c.m, _, c.state = collapse(apply(kick, c.state), Z(1); rng = c.rng)
            output(c, "inside" => [:step, Z(1)])
            log_message(c, "kicked at step $(c.step + 1)")
        end
        sd(name) = SimData(; name, seed = 5, phases = [phases[1], steps], final_measurements = "m" => :m)
        ref = runTMS(sd("ref"))
        @test length(readlines("ref/kicks")) == 4
        @test occursin("kicked at step 6", read("ref/log", String))
        run_through(sd("chk"); max_time = -1)
        for f in ("kicks", "inside", "m")
            @test read("chk/$f", String) == read("ref/$f", String)
        end
    end
end

@testset "What a phase refuses when it is written, and what the console shows" begin
    @test_throws "GroundState needs its search method: algo = Dmrg(nsweeps = 10)" GroundState(
        hamiltonian = Z(1))
    @test_throws "SteadyState needs its search method" SteadyState(lindbladian = Dissipator(Sm)(1))
    @test_throws "mpo_algo is \"naive\" or \"zipup\", not \"zipp\"" SteadyState(
        lindbladian = Dissipator(Sm)(1), algo = Dmrg(nsweeps = 1), mpo_algo = "zipp")
    in_tmp() do
        # a line once the run returns, where its results are, whether it stopped, and the
        # warnings its log received
        phases = [CreateState{Pure}(2, Qubit(), "Up"), Kicks(measurements = "data" => Z(1))]
        @test_logs (:info, "simulation \"sim\" done, results in $(abspath("sim"))") runTMS(
            SimData(; name = "sim", phases))
        @test_logs (:info, "simulation \"stop\" stopped before phase \"kicks\": run it again to " *
                           "resume it, results in $(abspath("stop"))") runTMS(
            SimData(; name = "stop", phases); max_time = -1)
        warned = Kicks(measurements = "c" => ["c" => c -> Check(1., 2.; tol = 0.5, warn = true)])
        @test_logs (:info, r"simulation \"warn\" done, results in .*warn, 4 warnings in its log$") runTMS(
            SimData(name = "warn", phases = [phases[1], warned]))
        @test count("WARNING: Check c failed", read("warn/log", String)) == 4
        @test_logs (:info, "simulation \"out\" done, 4 warnings in its output") runTMS(
            SimData(name = "out", phases = [phases[1], warned]); output = devnull)
    end
end

# a phase of one's own with a parameter, which its printing shows
struct Valued{N} <: AbstractPhase
    name::String
    time_start
    final_measurements
end

TMS.run_phase(::Valued, sim::Simulation, ::Int) = sim

@testset "Printing phases and simulations" begin
    @test startswith(repr(Valued{3}("valued", nothing, [])), "\nValued{3}(")
    @test startswith(repr(CreateState{Mixed}(2, Qubit(), "Up")), "\nCreateState{Mixed}(")
    @test occursin("phases =", repr(SimData(phases = [CreateState{Mixed}(2, Qubit(), "Up"),
                                                     Valued{3}("valued", nothing, [])])))
end

@testset "The phases that make and transform the state" begin
    in_tmp() do
        run(phases; kwargs...) = runTMS(SimData(; phases, kwargs...); output = devnull)
        # a State given to CreateState, mixed if asked, a mixed one not made pure, and a random
        # state without a system
        pure = State{Pure}(System(2, Qubit()), ["Up", "Dn"])
        @test run([CreateState(type = Mixed(), state = pure)]).state isa State{Mixed}
        @test_throws "asked for a pure state but given a mixed one" run(
            [CreateState(type = Pure(), state = mix(pure))])
        @test_throws "a system and randomize" run([CreateState(type = Pure(), randomize = 10)])
        @test_throws "cannot randomize a State into a mixed state" run(
            [CreateState(type = Mixed(), state = pure, randomize = 4)])
        # its time_start, and a random mixed state drawn in a sector
        @test run([CreateState{Pure}(2, Qubit(), "Up"; time_start = 1.)]).time == 1.
        conf = ["Occ", "Emp", "Occ", "Emp"]
        s = run([CreateState{Mixed}(4, Fermion(conserve = N), conf; randomize = 4)]).state
        @test s isa State{Mixed} && maxlinkdim(s) ≤ 4
        # ToMixed truncates a state already mixed
        @test maxlinkdim(run([CreateState{Mixed}(4, Qubit(), "Up"; randomize = 16),
                              ToMixed(limits = Limits(maxdim = 2))]).state) ≤ 2
        # a mixed GHZ state, then the partial trace keeping three of its sites
        sys = System(6, Qubit())
        ghz = mix((State{Pure}(sys, "Up") + State{Pure}(sys, "Dn")) / sqrt(2))
        check(name, f, ref) = name => c -> Check(f(c.state), ref; tol = 1e-10)
        sim = run([CreateState(type = Mixed(), state = ghz, final_measurements = "" => [
                       check("purity", purity, 1), check("x", st -> real(expect(st, prod(X(i) for i in 1:6))), 1)]),
                   PartialTrace(positions = [2, 3, 5], keep = true, final_measurements = "" => [
                       check("z", st -> real(expect1(st, Z)), zeros(3)),
                       check("zz", st -> real(expect2(st, (Z, Z))), ones(3, 3))])])
        @test length(sim.state) == 3
        # dephasing allows a strong symmetry and loss does not: weakening in between gives what
        # a weak symmetry gives from the start
        h = -sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3)
        dephasing = -im * h + sum(Dissipator(sqrt(0.3) * N)(i) for i in 1:4)
        loss = -im * h + Dissipator(sqrt(0.2) * C)(2)
        ev(evolver) = Evolve(; evolver, duration = 0.5, time_step = 0.05,
                             limits = Limits(cutoff = 1e-12, maxdim = 32))
        chain(site, middle) = run([CreateState{Mixed}(4, site, conf), ev(dephasing), middle...,
                                   ev(loss)]).state
        numbers(s) = [ real(trace(s)) ; real.(expect1(s, N)) ; expect(s, dag(C)(1) * C(3)) ]
        reference = numbers(chain(Fermion(conserve = N), []))
        s = chain(Fermion(conserve = strong(N)), [Weaken()])
        @test repr(symmetries(s.system)) == "N"
        @test numbers(s) ≈ reference
        s = chain(Fermion(conserve = strong(N)), [Weaken(target = ())])
        @test !TMS.is_charged(s.system)
        @test numbers(s) ≈ reference
        @test_throws "drop `strong`" chain(Fermion(conserve = strong(N)), [])
    end
end

@testset "A complex time, and checkpoints of another version" begin
    in_tmp() do
        # a complex time, of no imaginary part, keeps its two columns across a resume
        stop_in = Ref(0)
        phases = [CreateState{Pure}(2, Qubit(), "Up"),
                  Gates(gates = X(1), time_start = complex(0.5),
                        final_measurements = "data" => [Z(1), stopper(stop_in)]),
                  Gates(gates = X(1), final_measurements = "data" => [Z(1), stopper(stop_in)])]
        runTMS(SimData(; name = "ref", phases))
        @test startswith(readlines("ref/data")[1], "Re(time)\tIm(time)")
        stop_in[] = 1
        sd = SimData(; name = "chk", phases)
        @test stopped(runTMS(sd))
        k, _, _ = TMS.load_checkpoint(abspath("chk"))
        @test k.time === complex(0.5) && k.phase == 3
        @test !stopped(runTMS(sd))
        @test read("chk/data", String) == read("ref/data", String)
        # a checkpoint of another layout, or of another minor version of Julia, is refused
        stop_in[] = 1
        runTMS(SimData(; name = "old", phases))
        meta = "old/checkpoint.json"
        text = read(meta, String)
        write(meta, replace(text, "\"version\":5" => "\"version\":4"))
        @test_throws "has the layout 4 of another version" runTMS(SimData(; name = "old", phases))
        write(meta, replace(text, "\"julia\":\"$VERSION\"" => "\"julia\":\"1.9.0\""))
        @test_throws "written by Julia 1.9.0" runTMS(SimData(; name = "old", phases))
    end
end

@testset "A checkpoint is replaced whole" begin
    # a crash while a checkpoint is written leaves the one before it: the state file of the
    # other generation, half written, is not read, and is written again
    in_tmp() do
        stop_in = Ref(0)
        phases = [CreateState{Pure}(2, Qubit(), "Up"),
                  Kicks(measurements = "data" => [Z(1), stopper(stop_in)])]
        runTMS(SimData(; name = "ref", phases))
        stop_in[] = 2
        sd = SimData(; name = "chk", phases)
        runTMS(sd)
        g = TMS.load_checkpoint(abspath("chk"))[2]
        write(joinpath("chk", TMS.state_file(3 - g)), "left half written")
        stop_in[] = 1
        runTMS(sd)
        runTMS(sd)
        @test read("chk/data", String) == read("ref/data", String)
    end
end

@testset "Every kind of destination and of value survives a resume" begin
    # a text file is continued from where it was, a JSON Lines file too, and the rows of a Data,
    # held by the checkpoint, come back with their types, complex numbers and matrices included;
    # an interrupt in the final measurements of a phase of no steps resumes from the phase before
    in_tmp() do
        stop_in, fail_in = Ref(0), Ref(0)
        breaker = "breaker" => function (_)
            if fail_in[] > 0
                fail_in[] -= 1
                if fail_in[] == 0
                    throw(InterruptException())
                end
            end
            return 0.0
        end
        ms = [Sp(1), (X, Z), (Sp, Sm)]
        measurements = ["data" => [ms ; stopper(stop_in)], "out.json" => ms, Data("d") => ms]
        phases = [CreateState{Pure}(2, Qubit(), "X+"),
                  Evolve(evolver = -im * (Z(1) + X(1) * X(2)), duration = 0.6, time_step = 0.1,
                         limits = Limits(maxdim = 10, cutoff = 1e-15); measurements),
                  Gates(gates = Z(1)),
                  Gates(gates = X(1), final_measurements = "final" => [X(1), breaker])]
        ref = runTMS(SimData(; name = "ref", phases))
        stop_in[] = 3
        sd = SimData(; name = "chk", phases)
        @test stopped(runTMS(sd))
        fail_in[] = 1
        @test stopped(runTMS(sd))
        @test fail_in[] == 0
        sim = runTMS(sd)
        @test !stopped(sim)
        for f in ("data", "out.json", "final")
            @test read("chk/$f", String) == read("ref/$f", String)
        end
        resumed, whole = sim.data["d"], ref.data["d"]
        @test length(resumed) == length(whole)
        @test all(map(r -> map(typeof ∘ last, r), resumed) .== map(r -> map(typeof ∘ last, r), whole))
        @test all(isapprox(last(a), last(b)) for (ra, rb) in zip(resumed, whole) for (a, b) in zip(ra, rb))
    end
end

@testset "Nested phases, a dmrg resumed with a measurement period, the log and the files" begin
    in_tmp() do
        # nested vectors of phases are the flat list, resumed after each step alike
        k1 = Kicks(nkicks = 2, measurements = "data" => [Z(1), :step])
        flat = SimData(name = "flat", phases = [CreateState{Pure}(2, Qubit(), "Up"), k1, k1, k1])
        nested = SimData(name = "nested", phases = [[CreateState{Pure}(2, Qubit(), "Up")], [k1, [k1]], k1])
        runTMS(flat)
        run_through(nested; max_time = -1)
        @test read("nested/data", String) == read("flat/data", String)
        # a search resumed after each sweep logs the same sweeps, the last one included, which
        # its last step commits as the end of the phase
        h = -sum(Z(i) * Z(i + 1) for i in 1:3) - 0.8 * sum(X(i) for i in 1:4)
        phases = [CreateState{Pure}(4, Qubit(), "Up"),
                  GroundState(hamiltonian = h, algo = Dmrg(nsweeps = 5),
                              limits = Limits(maxdim = 8, cutoff = 1e-14),
                              measurements = "data" => [Z(1), :step], measurements_period = 2)]
        runTMS(SimData(; name = "gref", phases))
        run_through(SimData(; name = "gchk", phases); max_time = -1)
        @test read("gchk/data", String) == read("gref/data", String)
        steps(f) = [ parse(Int, m[1]) for m in eachmatch(r"\nstep (\d+), ", read(f, String)) ]
        energies(f) = [ parse(Float64, m[1]) for m in eachmatch(r", energy (\S+)\n", read(f, String)) ]
        @test steps("gchk/log") == steps("gref/log") == 1:5
        @test energies("gchk/log") ≈ energies("gref/log")
        # the time a resume logs is the one the checkpoint reached
        stop_in = Ref(2)
        sd = SimData(name = "time", phases = [CreateState{Pure}(2, Qubit(), "Up"),
                                              Kicks(measurements = "data" => [Z(1), stopper(stop_in)])])
        runTMS(sd)
        runTMS(sd)
        m = match(r"Resuming from the checkpoint: phase 2, step 2, time (\S+)", read("time/log", String))
        @test parse(Float64, m[1]) ≈ 0.2
        # with its checkpoints, the directory holds the files of the simulation and its outputs only
        runTMS(SimData(name = "all", phases = [CreateState{Pure}(2, Qubit(), "Up"), k1]);
               checkpoint_interval = 1e-9)
        @test issubset(setdiff(readdir("all"), ["data"]), TMS.simulation_files)
        @test "checkpoint.json" in readdir("all")
        for file in ("stop", "checkpoint.json", "./running")
            @test_throws "a file of the simulation directory" runTMS(SimData(name = "b",
                phases = [CreateState{Pure}(2, Qubit(), "Up")], final_measurements = file => Z(1)))
        end
    end
end

@testset "The warnings of a run" begin
    # a warning of the package goes to the log of the run, a warning of one's own to the logger
    # of the caller, as the line of the end does
    in_tmp() do
        warned = Kicks(nkicks = 1, measurements = "c" => [
            "c" => c -> Check(1., 2.; tol = 0.5, warn = true),
            "u" => c -> (@warn "from the user"; 0.)])
        @test_logs (:warn, "from the user") (:info, r"simulation \"w\" done") runTMS(
            SimData(name = "w", phases = [CreateState{Pure}(2, Qubit(), "Up"), warned]))
        @test occursin("WARNING: Check c failed", read("w/log", String))
        @test !occursin("from the user", read("w/log", String))
    end
end
